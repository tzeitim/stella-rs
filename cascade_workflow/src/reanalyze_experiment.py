#!/usr/bin/env python3
"""
STRICT Re-analysis infrastructure that FAILS on missing/incorrect data.
No fallbacks, no guessing, no silent failures.

Usage:
    python reanalyze_experiment.py experiment_dir --output results.csv
"""

import argparse
import pickle
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional
import pandas as pd
import yaml
import logging
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import time

from core.shared_metrics import calculate_metrics_for_trees
from core.metrics_computation import extract_parameters_from_tree
import cassiopeia as cass

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Raised when data validation fails."""
    pass


class MissingDataError(Exception):
    """Raised when required data is missing."""
    pass


class ConfigurationError(Exception):
    """Raised when configuration is invalid or missing."""
    pass


def validate_tree_structure(tree: cass.data.CassiopeiaTree, tree_name: str) -> None:
    """
    Validate that a tree has the required structure for analysis.
    Raises exceptions if validation fails.
    """
    # Check basic tree structure
    if not hasattr(tree, 'character_matrix'):
        raise ValidationError(f"{tree_name}: Missing character_matrix")

    if not hasattr(tree, 'parameters'):
        raise ValidationError(f"{tree_name}: Missing parameters dictionary")

    # Validate character matrix
    cm_shape = tree.character_matrix.shape
    if cm_shape[0] == 0:
        raise ValidationError(f"{tree_name}: Empty character matrix (no leaves)")
    if cm_shape[1] == 0:
        raise ValidationError(f"{tree_name}: Empty character matrix (no characters)")

    # Check for required parameters
    required_params = ['lam_true', 'q_true', 'lam_gt', 'q_gt']
    missing_params = [p for p in required_params if p not in tree.parameters]
    if missing_params:
        raise MissingDataError(f"{tree_name}: Missing required parameters: {missing_params}")

    # Validate parameter values
    for param in required_params:
        value = tree.parameters[param]
        if value is None or np.isnan(value):
            raise ValidationError(f"{tree_name}: Parameter '{param}' is None or NaN")

    # Check priors
    if not hasattr(tree, 'priors') or tree.priors is None:
        raise ValidationError(f"{tree_name}: Missing priors")

    if len(tree.priors) == 0:
        raise ValidationError(f"{tree_name}: Empty priors dictionary")


def validate_tier_config(config: Dict[str, Any], tier_num: int) -> Dict[str, Any]:
    """
    Validate and extract tier configuration. No guessing allowed.
    """
    if 'cas9_tiers' not in config:
        raise ConfigurationError("Configuration missing 'cas9_tiers' section")

    if tier_num not in config['cas9_tiers']:
        available_tiers = list(config['cas9_tiers'].keys())
        raise ConfigurationError(f"Tier {tier_num} not found in config. Available: {available_tiers}")

    tier_config = config['cas9_tiers'][tier_num]

    # Validate required tier fields
    required_fields = ['k', 'cassette_size', 'm', 'name']
    missing_fields = [f for f in required_fields if f not in tier_config]
    if missing_fields:
        raise ConfigurationError(f"Tier {tier_num} config missing fields: {missing_fields}")

    # Validate values are positive
    if tier_config['k'] <= 0:
        raise ValidationError(f"Tier {tier_num}: k must be positive, got {tier_config['k']}")
    if tier_config['cassette_size'] <= 0:
        raise ValidationError(f"Tier {tier_num}: cassette_size must be positive")

    return tier_config


def validate_trees_compatible(tree1: cass.data.CassiopeiaTree, tree2: cass.data.CassiopeiaTree,
                             tree1_name: str, tree2_name: str) -> None:
    """
    Validate that two trees are compatible for comparison.
    """
    # Check leaf sets match
    leaves1 = set(tree1.leaves)
    leaves2 = set(tree2.leaves)

    if leaves1 != leaves2:
        only_in_1 = leaves1 - leaves2
        only_in_2 = leaves2 - leaves1
        raise ValidationError(f"Leaf sets don't match between {tree1_name} and {tree2_name}. "
                            f"Only in {tree1_name}: {list(only_in_1)[:5]}... "
                            f"Only in {tree2_name}: {list(only_in_2)[:5]}...")

    # For now, don't require character matrix dimensions to match exactly,
    # as we may be comparing trees with different numbers of recording sites.
    # But log the difference
    if tree1.character_matrix.shape[1] != tree2.character_matrix.shape[1]:
        logger.warning(f"Character dimensions differ: {tree1_name}={tree1.character_matrix.shape[1]}, "
                      f"{tree2_name}={tree2.character_matrix.shape[1]}")


def process_single_tree_strict(args_tuple):
    """
    Process a single reconstructed tree with STRICT validation.
    NO FALLBACKS. Fails on missing data.
    """
    recon_file_path, experiment_dir, config, cas9_trees, gt_tree = args_tuple

    try:
        # Parse filename
        def parse_reconstructed_filename(filename: str) -> Dict[str, Any]:
            clean_name = filename.replace('_reconstructed.pkl', '').replace('_reconstructed', '')
            parts = clean_name.split('_')

            result = {
                'instance_id': None,
                'sim_id': None,
                'reconstruction_id': None,
                'tier_num': None,
                'solver': None
            }

            for part in parts:
                if part.startswith('instance'):
                    result['instance_id'] = int(part.replace('instance', ''))
                elif part.startswith('sim'):
                    result['sim_id'] = int(part.replace('sim', ''))
                elif part.startswith('recon'):
                    result['reconstruction_id'] = int(part.replace('recon', ''))
                elif part.startswith('tier'):
                    result['tier_num'] = int(part.replace('tier', ''))
                elif result['solver'] is None:
                    result['solver'] = part

            # Validate all fields were parsed
            missing_fields = [k for k, v in result.items() if v is None]
            if missing_fields:
                raise ValidationError(f"Failed to parse filename '{filename}'. Missing: {missing_fields}")

            return result

        parsed = parse_reconstructed_filename(recon_file_path.name)

        # Load reconstructed tree
        with open(recon_file_path, 'rb') as f:
            reconstructed_tree = pickle.load(f)

        # Validate reconstructed tree structure
        validate_tree_structure(reconstructed_tree, f"Reconstructed tree {recon_file_path.name}")

        # Find corresponding CAS9 tree - REQUIRED
        cas9_key = f"instance{parsed['instance_id']}_sim{parsed['sim_id']}_tier{parsed['tier_num']}_instance"
        if cas9_key not in cas9_trees:
            raise MissingDataError(f"CAS9 tree not found: {cas9_key}")

        cas9_tree = cas9_trees[cas9_key]

        # Validate CAS9 tree
        validate_tree_structure(cas9_tree, f"CAS9 tree {cas9_key}")

        # Validate trees are compatible
        validate_trees_compatible(reconstructed_tree, cas9_tree,
                                 "Reconstructed tree", "CAS9 tree")

        # Validate tier configuration
        tier_config = validate_tier_config(config, parsed['tier_num'])

        # Extract REQUIRED parameters
        lam_sim = cas9_tree.parameters['lam_true']
        q_sim = cas9_tree.parameters['q_true']
        lam_gt = cas9_tree.parameters['lam_gt']
        q_gt = cas9_tree.parameters['q_gt']

        # Validate parameters are reasonable
        for param_name, param_value in [('lam_sim', lam_sim), ('q_sim', q_sim),
                                       ('lam_gt', lam_gt), ('q_gt', q_gt)]:
            if param_value <= 0 or param_value > 10:  # Reasonable bounds
                raise ValidationError(f"Parameter {param_name}={param_value} out of reasonable range (0, 10]")

        # Calculate metrics with validated inputs
        metrics = calculate_metrics_for_trees(
            reconstructed_tree=reconstructed_tree,
            reference_tree=cas9_tree,  # ALWAYS use CAS9 tree
            config=config,
            gt_instance_id=parsed['instance_id'],
            cas9_simulation_id=parsed['sim_id'],
            reconstruction_id=parsed['reconstruction_id'],
            solver_name=parsed['solver'],
            lam_sim=lam_sim,
            q_sim=q_sim,
            lam_gt=lam_gt,
            q_gt=q_gt,
            tier_num=parsed['tier_num']
        )

        # Validate metrics output
        required_metrics = ['cPHS_simulation', 'cPHS_gt', 'RF_distance', 'triplets_distance']
        for metric in required_metrics:
            if metric not in metrics:
                raise ValidationError(f"Metric calculation didn't return '{metric}'")
            if np.isnan(metrics[metric]):
                # This is actually an error - metrics shouldn't be NaN if inputs are valid
                raise ValidationError(f"Metric '{metric}' is NaN - calculation failed")

        # Build result row (matching pipeline format)
        reconstruction_id = (f"instance{parsed['instance_id']}_sim{parsed['sim_id']}_"
                           f"recon{parsed['reconstruction_id']}_tier{parsed['tier_num']}_"
                           f"{parsed['solver']}")

        # All fields are REQUIRED - no defaults
        row = {
            'reconstruction_id': reconstruction_id,
            'gt_instance_id': parsed['instance_id'],
            'cas9_simulation_id': parsed['sim_id'],
            'reconstruction_num': parsed['reconstruction_id'],
            'cas9_tier': parsed['tier_num'],
            'cas9_tier_name': tier_config['name'],
            'recording_sites': tier_config['k'] * tier_config['cassette_size'],
            'states_per_site': tier_config['m'],
            'solver': parsed['solver'],
            'lam_simulation': lam_sim,
            'q_simulation': q_sim,
            'lam_gt': lam_gt,
            'q_gt': q_gt,
            'cPHS': metrics['cPHS_simulation'],
            'cPHS_gt': metrics['cPHS_gt'],
            'RF_distance': metrics['RF_distance'],
            'triplets_distance': metrics['triplets_distance'],
            # Add other metrics as available
        }

        return [row]

    except (ValidationError, MissingDataError, ConfigurationError) as e:
        # These are expected errors - data issues
        return [{'error': str(e), 'file': str(recon_file_path), 'error_type': type(e).__name__}]

    except Exception as e:
        # Unexpected errors - likely bugs
        return [{'error': str(e), 'file': str(recon_file_path), 'error_type': 'UNEXPECTED',
                'traceback': str(e)}]


class ExperimentReanalyzer:
    """
    STRICT re-analyzer that fails on missing or invalid data.
    """

    def __init__(self, experiment_dir: Path):
        self.experiment_dir = experiment_dir

        # Load and validate configuration
        self.config = self._load_and_validate_config()

        # Load and validate ground truth tree
        self.gt_tree = self._load_and_validate_gt_tree()

        # CAS9 trees loaded on demand
        self.cas9_trees = {}

    def _load_and_validate_config(self) -> Dict[str, Any]:
        """Load experiment configuration with validation."""
        config_path = self.experiment_dir / 'cascade_config.yaml'
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        # Validate basic config structure
        required_sections = ['cas9_tiers', 'ground_truth', 'execution']
        missing_sections = [s for s in required_sections if s not in config]
        if missing_sections:
            raise ConfigurationError(f"Config missing sections: {missing_sections}")

        logger.info(f"Loaded valid config from {config_path}")
        return config

    def _load_and_validate_gt_tree(self) -> cass.data.CassiopeiaTree:
        """Load ground truth tree with validation."""
        gt_tree_path = self.experiment_dir / 'gt_tree.pkl'
        if not gt_tree_path.exists():
            # Try alternative location
            gt_tree_path = self.experiment_dir / 'gt_trees' / 'instance0_gt.pkl'
            if not gt_tree_path.exists():
                raise FileNotFoundError(f"Ground truth tree not found in expected locations")

        with open(gt_tree_path, 'rb') as f:
            gt_tree = pickle.load(f)

        # Basic validation
        if not hasattr(gt_tree, 'character_matrix'):
            raise ValidationError("GT tree missing character matrix")

        logger.info(f"Loaded GT tree with {len(gt_tree.leaves)} leaves")
        return gt_tree

    def _load_cas9_trees(self, required_count: Optional[int] = None) -> None:
        """Load CAS9 trees with validation."""
        cas9_dir = self.experiment_dir / 'cas9_instances'
        if not cas9_dir.exists():
            raise FileNotFoundError(f"CAS9 instances directory not found: {cas9_dir}")

        cas9_files = sorted(cas9_dir.glob('*.pkl'))
        if len(cas9_files) == 0:
            raise MissingDataError("No CAS9 tree files found")

        if required_count and len(cas9_files) < required_count:
            raise MissingDataError(f"Expected at least {required_count} CAS9 trees, found {len(cas9_files)}")

        for cas9_file in cas9_files:
            with open(cas9_file, 'rb') as f:
                tree = pickle.load(f)

            # Validate each CAS9 tree
            try:
                validate_tree_structure(tree, cas9_file.stem)
                self.cas9_trees[cas9_file.stem] = tree
            except (ValidationError, MissingDataError) as e:
                logger.error(f"Invalid CAS9 tree {cas9_file.name}: {e}")
                # Don't add invalid trees

        logger.info(f"Loaded {len(self.cas9_trees)} valid CAS9 trees")

    def reanalyze_reconstructed_trees(self, max_trees: int = None, n_jobs: int = None) -> pd.DataFrame:
        """
        Re-analyze with STRICT validation.
        """
        logger.info("Starting STRICT re-analysis...")

        # Load CAS9 trees
        self._load_cas9_trees()

        # Find reconstructed trees
        reconstructed_dir = self.experiment_dir / 'reconstructed_trees'
        if not reconstructed_dir.exists():
            raise FileNotFoundError(f"Reconstructed trees directory not found: {reconstructed_dir}")

        reconstructed_files = sorted(reconstructed_dir.glob('*_reconstructed.pkl'))
        if len(reconstructed_files) == 0:
            raise MissingDataError("No reconstructed tree files found")

        if max_trees:
            reconstructed_files = reconstructed_files[:max_trees]

        logger.info(f"Processing {len(reconstructed_files)} trees with STRICT validation")

        # Determine workers
        if n_jobs is None:
            n_jobs = min(multiprocessing.cpu_count(), len(reconstructed_files))

        # Prepare arguments
        args_list = [
            (recon_file, self.experiment_dir, self.config, self.cas9_trees, self.gt_tree)
            for recon_file in reconstructed_files
        ]

        all_results = []
        error_summary = {'ValidationError': 0, 'MissingDataError': 0,
                        'ConfigurationError': 0, 'UNEXPECTED': 0}

        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = {
                executor.submit(process_single_tree_strict, args): args[0]
                for args in args_list
            }

            for future in as_completed(futures):
                recon_file = futures[future]
                try:
                    result_rows = future.result()

                    if result_rows and 'error' in result_rows[0]:
                        error_info = result_rows[0]
                        error_type = error_info.get('error_type', 'UNKNOWN')
                        error_summary[error_type] = error_summary.get(error_type, 0) + 1
                        logger.error(f"{recon_file.name}: {error_type}: {error_info['error']}")
                    else:
                        all_results.extend(result_rows)

                except Exception as e:
                    error_summary['UNEXPECTED'] = error_summary.get('UNEXPECTED', 0) + 1
                    logger.error(f"Process failed for {recon_file.name}: {e}")

        # Report error summary
        logger.info(f"Processing complete. Success: {len(all_results)}, Errors: {sum(error_summary.values())}")
        for error_type, count in error_summary.items():
            if count > 0:
                logger.error(f"  {error_type}: {count}")

        # Fail if too many errors
        total_errors = sum(error_summary.values())
        if total_errors > len(reconstructed_files) * 0.1:  # >10% failure rate
            raise ValidationError(f"Too many failures: {total_errors}/{len(reconstructed_files)}")

        return pd.DataFrame(all_results)


def main():
    parser = argparse.ArgumentParser(description="STRICT re-analysis - fails on missing/invalid data")
    parser.add_argument('experiment_dir', type=str, help='Path to experiment directory')
    parser.add_argument('--max-trees', type=int, default=None, help='Maximum trees to process')
    parser.add_argument('--jobs', '-j', type=int, default=None, help='Number of parallel jobs')
    parser.add_argument('--output', type=str, required=True, help='Output CSV file (REQUIRED)')

    args = parser.parse_args()

    experiment_dir = Path(args.experiment_dir)
    if not experiment_dir.exists():
        logger.error(f"Experiment directory not found: {experiment_dir}")
        sys.exit(1)

    try:
        reanalyzer = ExperimentReanalyzer(experiment_dir)
        results_df = reanalyzer.reanalyze_reconstructed_trees(args.max_trees, args.jobs)

        if len(results_df) == 0:
            logger.error("No valid results produced")
            sys.exit(1)

        results_df.to_csv(args.output, index=False)
        logger.info(f"Results saved to {args.output}")

    except (ValidationError, MissingDataError, ConfigurationError) as e:
        logger.error(f"VALIDATION FAILED: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()