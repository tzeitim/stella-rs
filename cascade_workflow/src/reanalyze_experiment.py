#!/usr/bin/env python3
"""
Re-analysis infrastructure using the main source code.
This script generates experiment-specific re-analysis using the canonical
metrics computation functions from the cascade workflow.

Usage:
    python reanalyze_experiment.py experiment_dir --mode [reconstructed|cas9|both]
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


class ValidationError(Exception):
    """Raised when data validation fails"""
    pass


class MissingDataError(Exception):
    """Raised when required data is missing"""
    pass


class ConfigurationError(Exception):
    """Raised when configuration is invalid or missing"""
    pass

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def process_single_tree(args_tuple):
    """
    Process a single reconstructed tree. This function is designed to be used with multiprocessing.

    Args:
        args_tuple: Tuple containing (recon_file_path, experiment_dir, config, cas9_trees, gt_tree)

    Returns:
        List of dictionaries containing metrics for this tree
    """
    recon_file_path, experiment_dir, config, cas9_trees, gt_tree = args_tuple

    try:
        # Parse filename (copied from main class)
        def parse_reconstructed_filename(filename: str) -> Dict[str, Any]:
            clean_name = filename.replace('_reconstructed.pkl', '').replace('_reconstructed', '')
            parts = clean_name.split('_')

            result = {
                'instance_id': 0,
                'sim_id': 0,
                'reconstruction_id': 0,
                'tier_num': None,
                'solver': None
            }

            solver_found = False
            for part in parts:
                if part.startswith('instance'):
                    try:
                        result['instance_id'] = int(part.replace('instance', ''))
                    except:
                        pass
                elif part.startswith('sim'):
                    try:
                        result['sim_id'] = int(part.replace('sim', ''))
                    except:
                        pass
                elif part.startswith('recon'):
                    try:
                        result['reconstruction_id'] = int(part.replace('recon', ''))
                    except:
                        pass
                elif part.startswith('tier'):
                    try:
                        result['tier_num'] = int(part.replace('tier', ''))
                    except:
                        pass
                elif not solver_found and not any(part.startswith(x) for x in ['instance', 'sim', 'recon', 'tier']):
                    result['solver'] = part
                    solver_found = True

            return result

        def get_tier_info(tier_num: int, config: Dict[str, Any]) -> Dict[str, Any]:
            tier_config = config.get('cas9_tiers', {}).get(tier_num, {})
            return {
                'cas9_tier': tier_num,
                'cas9_tier_name': tier_config.get('name', f'Tier {tier_num}'),
                'recording_sites': tier_config.get('k', 10) * tier_config.get('cassette_size', 10),
                'states_per_site': tier_config.get('m', 20),
                'gt_tree_size': config.get('ground_truth', {}).get('tree_config', {}).get('N', 1000),
                'sampled_tree_size': config.get('ground_truth', {}).get('tree_config', {}).get('n', 150),
                'run_name': config.get('execution', {}).get('run_name', 'unknown'),
                'experiment_id': config.get('execution', {}).get('run_name', 'unknown')
            }

        # Parse filename
        parsed = parse_reconstructed_filename(recon_file_path.name)
        if parsed['tier_num'] is None or parsed['solver'] is None:
            return []

        # Load reconstructed tree (handle both old and new formats)
        with open(recon_file_path, 'rb') as f:
            loaded_data = pickle.load(f)

        # Handle new format (tree + metadata) vs old format (just tree)
        if isinstance(loaded_data, dict) and 'tree' in loaded_data and 'metadata' in loaded_data:
            # New format: extract tree and use stored reconstruction_id
            reconstructed_tree = loaded_data['tree']
            stored_metadata = loaded_data['metadata']
            # Override parsed reconstruction_id with stored one if available
            if 'base_reconstruction_id' in stored_metadata:
                base_id = stored_metadata['base_reconstruction_id']
                logger.info(f"Using stored reconstruction_id: {base_id}")
        else:
            # Old format: just the tree
            reconstructed_tree = loaded_data

        # Find corresponding CAS9 tree
        cas9_key = f"instance{parsed['instance_id']}_sim{parsed['sim_id']}_tier{parsed['tier_num']}_instance"
        cas9_tree = cas9_trees.get(cas9_key)
        if cas9_tree is None:
            raise MissingDataError(f"Required CAS9 tree missing: {cas9_key}")

        # Get tier information
        tier_info = get_tier_info(parsed['tier_num'], config)

        # Extract parameters - REQUIRED for PHS calculation
        if not hasattr(cas9_tree, 'parameters'):
            raise MissingDataError(f"CAS9 tree missing parameters attribute: {cas9_key}")

        lam_sim = cas9_tree.parameters.get('lam_true')
        q_sim = cas9_tree.parameters.get('q_true')
        lam_gt = cas9_tree.parameters.get('lam_gt')
        q_gt = cas9_tree.parameters.get('q_gt')

        if lam_sim is None or q_sim is None:
            raise MissingDataError(f"Critical simulation parameters missing for {cas9_key}: lam_sim={lam_sim}, q_sim={q_sim}")

        # Use the same metrics calculation as the pipeline
        metrics = calculate_metrics_for_trees(
            reconstructed_tree=reconstructed_tree,
            reference_tree=cas9_tree if cas9_tree else gt_tree,
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

        # Create reconstruction ID to match pipeline format
        reconstruction_id = (f"instance{parsed['instance_id']}_sim{parsed['sim_id']}_"
                           f"recon{parsed['reconstruction_id']}_tier{parsed['tier_num']}_"
                           f"{parsed['solver']}")

        # Extract tier info
        tier_info = get_tier_info(parsed['tier_num'], config)

        # Create row matching exact pipeline format
        row = {
            # Basic identification (columns 0-3)
            'reconstruction_id': reconstruction_id,
            'gt_instance_id': parsed['instance_id'],
            'cas9_simulation_id': parsed['sim_id'],
            'reconstruction_num': parsed['reconstruction_id'],

            # Tier information (columns 4-7)
            'cas9_tier': parsed['tier_num'],
            'cas9_tier_name': tier_info['cas9_tier_name'],
            'recording_sites': tier_info['recording_sites'],
            'states_per_site': tier_info['states_per_site'],

            # Solver and timing (columns 8-9)
            'solver': parsed['solver'],
            'computation_time_seconds': metrics.get('computation_time_ms', 0) / 1000.0,

            # Experiment metadata (columns 10-13)
            'run_name': tier_info['run_name'],
            'experiment_id': tier_info['experiment_id'],
            'gt_tree_size': tier_info['gt_tree_size'],
            'sampled_tree_size': tier_info['sampled_tree_size'],

            # Simulation parameters (columns 14-17)
            'lam_simulation': metrics.get('lam_from_simulation', lam_sim),
            'q_simulation': metrics.get('q_from_simulation', q_sim),
            'proportion_mutated_simulation': metrics.get('proportion_mutated_simulation', np.nan),
            'q_simulation_source': metrics.get('q_simulation_source', 'simulation'),

            # Ground truth parameters (columns 18-20)
            'lam_gt': lam_gt,
            'q_gt': q_gt,
            'proportion_mutated_gt': metrics.get('proportion_mutated_gt', np.nan),

            # PHS source information (columns 21-22)
            'phs_lam_source': 'simulation',
            'phs_q_source': 'simulation',

            # Distance metrics (columns 23-24)
            'triplets_distance': metrics.get('triplets_distance', np.nan),
            'RF_distance': metrics.get('RF_distance', np.nan),

            # Parsimony metrics (columns 25-28)
            'parsimony_total_mutations': np.nan,
            'parsimony_computation_time_ms': np.nan,
            'parsimony_method_used': np.nan,
            'parsimony_internal_states_inferred': np.nan,

            # PHS scores (columns 29-30)
            'cPHS': metrics.get('cPHS_simulation', np.nan),
            'cPHS_gt': metrics.get('cPHS_gt', np.nan),

            # Likelihood timing (column 31)
            'likelihood_computation_time_ms': np.nan,

            # Worker information (columns 32-35)
            'worker_cas9_instance_path': '',
            'worker_solver': parsed['solver'],
            'worker_tier': parsed['tier_num'],
            'worker_tier_name': tier_info['cas9_tier_name'],

            # Final scores (columns 36-39)
            'parsimony_score': np.nan,
            'log_likelihood': np.nan,
            'log_likelihood_simulation': np.nan,
            'log_likelihood_gt': np.nan
        }

        # Handle parsimony score if available
        if isinstance(metrics.get('parsimony_score'), dict):
            parsimony_dict = metrics['parsimony_score']
            row.update({
                'parsimony_total_mutations': parsimony_dict.get('total_mutations', np.nan),
                'parsimony_computation_time_ms': parsimony_dict.get('computation_time_ms', np.nan),
                'parsimony_method_used': parsimony_dict.get('method_used', ''),
                'parsimony_internal_states_inferred': parsimony_dict.get('internal_states_inferred', np.nan),
                'parsimony_score': parsimony_dict.get('parsimony_score', np.nan)
            })

        # Handle likelihood scores if available
        if isinstance(metrics.get('likelihood_score'), dict):
            likelihood_dict = metrics['likelihood_score']
            row.update({
                'log_likelihood': likelihood_dict.get('log_likelihood', np.nan)
            })
        if isinstance(metrics.get('likelihood_score_simulation'), dict):
            likelihood_sim_dict = metrics['likelihood_score_simulation']
            row.update({
                'log_likelihood_simulation': likelihood_sim_dict.get('log_likelihood', np.nan)
            })
        if isinstance(metrics.get('likelihood_score_gt'), dict):
            likelihood_gt_dict = metrics['likelihood_score_gt']
            row.update({
                'log_likelihood_gt': likelihood_gt_dict.get('log_likelihood', np.nan)
            })

        return [row]

    except (ValidationError, MissingDataError, ConfigurationError) as e:
        # Re-raise validation errors to fail fast
        raise
    except Exception as e:
        # Return error info for unexpected errors
        return [{'error': str(e), 'file': str(recon_file_path)}]


class ExperimentReanalyzer:
    """
    Re-analyzes experiment results using the main source code.
    Ensures consistency with the original pipeline implementation.
    """

    def __init__(self, experiment_dir: Path, config: Optional[Dict[str, Any]] = None):
        self.experiment_dir = experiment_dir
        self.config = config or {}
        self.gt_tree = None
        self.cas9_trees = {}

        # Load experiment configuration
        self._load_experiment_config()

        # Load ground truth tree
        self._load_ground_truth_tree()

    def _load_experiment_config(self):
        """Load experiment configuration from cascade_config.yaml."""
        config_path = self.experiment_dir / 'cascade_config.yaml'
        if config_path.exists():
            with open(config_path, 'r') as f:
                self.config = yaml.safe_load(f)
            logger.info(f"Loaded experiment config: {config_path}")
        else:
            logger.warning(f"No config found at {config_path}")

    def _load_ground_truth_tree(self):
        """Load the ground truth tree."""
        gt_tree_path = self.experiment_dir / 'gt_tree.pkl'
        if gt_tree_path.exists():
            with open(gt_tree_path, 'rb') as f:
                self.gt_tree = pickle.load(f)
            logger.info("Loaded ground truth tree")
        else:
            raise FileNotFoundError(f"Ground truth tree not found: {gt_tree_path}")

    def _load_cas9_trees(self, max_trees: int = None):
        """Load CAS9 trees for parameter extraction."""
        cas9_dir = self.experiment_dir / 'cas9_instances'
        if not cas9_dir.exists():
            logger.warning("CAS9 instances directory not found")
            return

        cas9_files = sorted(cas9_dir.glob('*.pkl'))
        if max_trees:
            cas9_files = cas9_files[:max_trees]

        for cas9_file in cas9_files:
            with open(cas9_file, 'rb') as f:
                self.cas9_trees[cas9_file.stem] = pickle.load(f)

        logger.info(f"Loaded {len(self.cas9_trees)} CAS9 trees")

    def _parse_reconstructed_filename(self, filename: str) -> Dict[str, Any]:
        """Parse reconstructed tree filename to extract metadata."""
        clean_name = filename.replace('_reconstructed.pkl', '').replace('_reconstructed', '')
        parts = clean_name.split('_')

        result = {
            'instance_id': 0,
            'sim_id': 0,
            'reconstruction_id': 0,
            'tier_num': None,
            'solver': None
        }

        solver_found = False
        for part in parts:
            if part.startswith('instance'):
                try:
                    result['instance_id'] = int(part.replace('instance', ''))
                except:
                    pass
            elif part.startswith('sim'):
                try:
                    result['sim_id'] = int(part.replace('sim', ''))
                except:
                    pass
            elif part.startswith('recon'):
                try:
                    result['reconstruction_id'] = int(part.replace('recon', ''))
                except:
                    pass
            elif part.startswith('tier'):
                try:
                    result['tier_num'] = int(part.replace('tier', ''))
                except:
                    pass
            elif not solver_found and not any(part.startswith(x) for x in ['instance', 'sim', 'recon', 'tier']):
                result['solver'] = part
                solver_found = True

        return result

    def _get_tier_info(self, tier_num: int) -> Dict[str, Any]:
        """Get tier configuration information."""
        tier_config = self.config.get('cas9_tiers', {}).get(tier_num, {})

        return {
            'cas9_tier': tier_num,
            'cas9_tier_name': tier_config.get('name', f'Tier {tier_num}'),
            'recording_sites': tier_config.get('k', 10) * tier_config.get('cassette_size', 10),
            'states_per_site': tier_config.get('m', 20),
            'gt_tree_size': self.config.get('ground_truth', {}).get('tree_config', {}).get('N', 1000),
            'sampled_tree_size': self.config.get('ground_truth', {}).get('tree_config', {}).get('n', 150),
            'run_name': self.config.get('execution', {}).get('run_name', 'unknown'),
            'experiment_id': self.config.get('execution', {}).get('run_name', 'unknown')
        }

    def reanalyze_reconstructed_trees(self, max_trees: int = None, n_jobs: int = None) -> pd.DataFrame:
        """
        Re-analyze reconstructed trees using the main source metrics computation.

        Args:
            max_trees: Maximum number of trees to process
            n_jobs: Number of parallel processes (None for auto-detect)
        """
        logger.info("Re-analyzing reconstructed trees...")

        # Load CAS9 trees for parameter extraction
        self._load_cas9_trees()

        # Find reconstructed trees
        reconstructed_dir = self.experiment_dir / 'reconstructed_trees'
        if not reconstructed_dir.exists():
            logger.error(f"Reconstructed trees directory not found: {reconstructed_dir}")
            return pd.DataFrame()

        reconstructed_files = sorted(reconstructed_dir.glob('*_reconstructed.pkl'))
        if max_trees:
            reconstructed_files = reconstructed_files[:max_trees]

        logger.info(f"Processing {len(reconstructed_files)} reconstructed trees")

        # Determine number of workers
        if n_jobs is None:
            n_jobs = min(multiprocessing.cpu_count(), len(reconstructed_files))

        logger.info(f"Using {n_jobs} parallel processes")

        # Prepare arguments for parallel processing
        args_list = [
            (recon_file, self.experiment_dir, self.config, self.cas9_trees, self.gt_tree)
            for recon_file in reconstructed_files
        ]

        all_results = []
        completed_count = 0
        error_count = 0
        start_time = time.time()

        # Process trees in parallel
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            # Submit all jobs
            future_to_file = {
                executor.submit(process_single_tree, args): args[0]
                for args in args_list
            }

            # Collect results as they complete
            for future in as_completed(future_to_file):
                recon_file = future_to_file[future]
                try:
                    result_rows = future.result()

                    # Check for errors in individual processing
                    if result_rows and 'error' in result_rows[0]:
                        logger.error(f"Error processing {recon_file.name}: {result_rows[0]['error']}")
                        error_count += 1
                    else:
                        all_results.extend(result_rows)
                        completed_count += 1

                        # Log progress
                        if completed_count % 10 == 0:
                            elapsed = time.time() - start_time
                            rate = completed_count / elapsed
                            remaining = len(reconstructed_files) - completed_count - error_count
                            eta = remaining / rate if rate > 0 else 0
                            logger.info(f"Completed {completed_count}/{len(reconstructed_files)} trees "
                                      f"({rate:.1f} trees/sec, ETA: {eta:.0f}s)")

                except (ValidationError, MissingDataError, ConfigurationError) as e:
                    # Validation errors should stop processing
                    logger.error(f"Validation failed for {recon_file.name}: {e}")
                    raise
                except Exception as e:
                    logger.error(f"Failed to process {recon_file.name}: {e}")
                    error_count += 1

        elapsed_total = time.time() - start_time
        logger.info(f"Parallel processing completed: {completed_count} trees processed, "
                   f"{error_count} errors in {elapsed_total:.1f}s")

        # Convert to DataFrame and add relative parsimony
        df = pd.DataFrame(all_results)
        if not df.empty:
            df = self._add_relative_parsimony(df)

        return df

    def reanalyze_reconstructed_trees_sequential(self, max_trees: int = None) -> pd.DataFrame:
        """
        Sequential version of reanalyze_reconstructed_trees for comparison/debugging.
        """
        logger.info("Re-analyzing reconstructed trees (sequential)...")

        # Load CAS9 trees for parameter extraction
        self._load_cas9_trees()

        # Find reconstructed trees
        reconstructed_dir = self.experiment_dir / 'reconstructed_trees'
        if not reconstructed_dir.exists():
            logger.error(f"Reconstructed trees directory not found: {reconstructed_dir}")
            return pd.DataFrame()

        reconstructed_files = sorted(reconstructed_dir.glob('*_reconstructed.pkl'))
        if max_trees:
            reconstructed_files = reconstructed_files[:max_trees]

        logger.info(f"Processing {len(reconstructed_files)} reconstructed trees")

        all_results = []

        for recon_file in reconstructed_files:
            logger.info(f"Processing: {recon_file.name}")

            # Parse filename
            parsed = self._parse_reconstructed_filename(recon_file.name)
            if parsed['tier_num'] is None or parsed['solver'] is None:
                logger.warning(f"Skipping {recon_file.name}: couldn't parse tier/solver")
                continue

            # Load reconstructed tree
            with open(recon_file, 'rb') as f:
                reconstructed_tree = pickle.load(f)

            # Find corresponding CAS9 tree
            cas9_key = f"instance{parsed['instance_id']}_sim{parsed['sim_id']}_tier{parsed['tier_num']}_instance"
            cas9_tree = self.cas9_trees.get(cas9_key)

            if cas9_tree is None:
                raise MissingDataError(f"Required CAS9 tree missing: {cas9_key}")

            # Get tier information
            tier_info = self._get_tier_info(parsed['tier_num'])

            # Compute all metrics using the shared pipeline implementation
            try:
                # Extract parameters - REQUIRED for PHS calculation
                if not hasattr(cas9_tree, 'parameters'):
                    raise MissingDataError(f"CAS9 tree missing parameters attribute: {cas9_key}")

                lam_sim = cas9_tree.parameters.get('lam_true')
                q_sim = cas9_tree.parameters.get('q_true')
                lam_gt = cas9_tree.parameters.get('lam_gt')
                q_gt = cas9_tree.parameters.get('q_gt')

                if lam_sim is None or q_sim is None:
                    raise MissingDataError(f"Critical simulation parameters missing for {cas9_key}: lam_sim={lam_sim}, q_sim={q_sim}")

                # Use the same metrics calculation as the pipeline
                metrics = calculate_metrics_for_trees(
                    reconstructed_tree=reconstructed_tree,
                    reference_tree=cas9_tree if cas9_tree else self.gt_tree,
                    config=self.config,
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

                # Create reconstruction ID to match pipeline format
                reconstruction_id = (f"instance{parsed['instance_id']}_sim{parsed['sim_id']}_"
                                   f"recon{parsed['reconstruction_id']}_tier{parsed['tier_num']}_"
                                   f"{parsed['solver']}")

                # Get tier information
                tier_info = self._get_tier_info(parsed['tier_num'])

                # Create row matching exact pipeline format (same as parallel version)
                row = {
                    # Basic identification (columns 0-3)
                    'reconstruction_id': reconstruction_id,
                    'gt_instance_id': parsed['instance_id'],
                    'cas9_simulation_id': parsed['sim_id'],
                    'reconstruction_num': parsed['reconstruction_id'],

                    # Tier information (columns 4-7)
                    'cas9_tier': parsed['tier_num'],
                    'cas9_tier_name': tier_info['cas9_tier_name'],
                    'recording_sites': tier_info['recording_sites'],
                    'states_per_site': tier_info['states_per_site'],

                    # Solver and timing (columns 8-9)
                    'solver': parsed['solver'],
                    'computation_time_seconds': metrics.get('computation_time_ms', 0) / 1000.0,

                    # Experiment metadata (columns 10-13)
                    'run_name': tier_info['run_name'],
                    'experiment_id': tier_info['experiment_id'],
                    'gt_tree_size': tier_info['gt_tree_size'],
                    'sampled_tree_size': tier_info['sampled_tree_size'],

                    # Simulation parameters (columns 14-17)
                    'lam_simulation': metrics.get('lam_from_simulation', lam_sim),
                    'q_simulation': metrics.get('q_from_simulation', q_sim),
                    'proportion_mutated_simulation': metrics.get('proportion_mutated_simulation', np.nan),
                    'q_simulation_source': metrics.get('q_simulation_source', 'simulation'),

                    # Ground truth parameters (columns 18-20)
                    'lam_gt': lam_gt,
                    'q_gt': q_gt,
                    'proportion_mutated_gt': metrics.get('proportion_mutated_gt', np.nan),

                    # PHS source information (columns 21-22)
                    'phs_lam_source': 'simulation',
                    'phs_q_source': 'simulation',

                    # Distance metrics (columns 23-24)
                    'triplets_distance': metrics.get('triplets_distance', np.nan),
                    'RF_distance': metrics.get('RF_distance', np.nan),

                    # Parsimony metrics (columns 25-28)
                    'parsimony_total_mutations': np.nan,
                    'parsimony_computation_time_ms': np.nan,
                    'parsimony_method_used': np.nan,
                    'parsimony_internal_states_inferred': np.nan,

                    # PHS scores (columns 29-30)
                    'cPHS': metrics.get('cPHS_simulation', np.nan),
                    'cPHS_gt': metrics.get('cPHS_gt', np.nan),

                    # Likelihood timing (column 31)
                    'likelihood_computation_time_ms': np.nan,

                    # Worker information (columns 32-35)
                    'worker_cas9_instance_path': '',
                    'worker_solver': parsed['solver'],
                    'worker_tier': parsed['tier_num'],
                    'worker_tier_name': tier_info['cas9_tier_name'],

                    # Final scores (columns 36-39)
                    'parsimony_score': np.nan,
                    'log_likelihood': np.nan,
                    'log_likelihood_simulation': np.nan,
                    'log_likelihood_gt': np.nan
                }

                # Handle parsimony score if available
                if isinstance(metrics.get('parsimony_score'), dict):
                    parsimony_dict = metrics['parsimony_score']
                    row.update({
                        'parsimony_total_mutations': parsimony_dict.get('total_mutations', np.nan),
                        'parsimony_computation_time_ms': parsimony_dict.get('computation_time_ms', np.nan),
                        'parsimony_method_used': parsimony_dict.get('method_used', ''),
                        'parsimony_internal_states_inferred': parsimony_dict.get('internal_states_inferred', np.nan),
                        'parsimony_score': parsimony_dict.get('parsimony_score', np.nan)
                    })

                # Handle likelihood scores if available
                if isinstance(metrics.get('likelihood_score'), dict):
                    likelihood_dict = metrics['likelihood_score']
                    row.update({
                        'log_likelihood': likelihood_dict.get('log_likelihood', np.nan)
                    })
                if isinstance(metrics.get('likelihood_score_simulation'), dict):
                    likelihood_sim_dict = metrics['likelihood_score_simulation']
                    row.update({
                        'log_likelihood_simulation': likelihood_sim_dict.get('log_likelihood', np.nan)
                    })
                if isinstance(metrics.get('likelihood_score_gt'), dict):
                    likelihood_gt_dict = metrics['likelihood_score_gt']
                    row.update({
                        'log_likelihood_gt': likelihood_gt_dict.get('log_likelihood', np.nan)
                    })

                rows = [row]

                all_results.extend(rows)

                rf_val = metrics.get('RF_distance', 'N/A')
                cphs_val = metrics.get('cPHS_simulation', 'N/A')
                rf_str = f"{rf_val:.3f}" if isinstance(rf_val, (int, float)) else str(rf_val)
                cphs_str = f"{cphs_val:.6f}" if isinstance(cphs_val, (int, float)) else str(cphs_val)
                logger.debug(f"Computed metrics: RF={rf_str}, cPHS={cphs_str}")

            except (ValidationError, MissingDataError, ConfigurationError) as e:
                # Validation errors should stop processing
                logger.error(f"Validation failed for {recon_file.name}: {e}")
                raise
            except Exception as e:
                logger.error(f"Failed to compute metrics for {recon_file.name}: {e}")

        # Convert to DataFrame and add relative parsimony
        df = pd.DataFrame(all_results)
        if not df.empty:
            df = self._add_relative_parsimony(df)

        return df

    def _add_relative_parsimony(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add relative parsimony scores (solver parsimony / minimum solver parsimony for each tree instance)."""
        if df.empty or 'parsimony_score' not in df.columns:
            logger.warning("Cannot compute relative parsimony: empty DataFrame or missing parsimony_score column")
            df['relative_parsimony'] = np.nan
            return df

        logger.info("Computing relative parsimony scores")

        # Group by tree instance (gt_instance_id, cas9_simulation_id, cas9_tier)
        # These columns identify unique tree instances across different solvers
        grouping_cols = ['gt_instance_id', 'cas9_simulation_id', 'cas9_tier']

        # Check if all grouping columns exist
        missing_cols = [col for col in grouping_cols if col not in df.columns]
        if missing_cols:
            logger.warning(f"Missing grouping columns for relative parsimony: {missing_cols}")
            df['relative_parsimony'] = np.nan
            return df

        # Compute minimum parsimony for each tree instance
        min_parsimony = df.groupby(grouping_cols)['parsimony_score'].transform('min')

        # Compute relative parsimony (solver parsimony / best solver parsimony)
        df['relative_parsimony'] = df['parsimony_score'] / min_parsimony

        # Log statistics
        relative_stats = df.groupby('solver')['relative_parsimony'].agg(['mean', 'std', 'min', 'max', 'count'])
        logger.info(f"Relative parsimony statistics by solver:\n{relative_stats}")

        return df

    def reanalyze_cas9_trees(self, max_trees: int = None) -> pd.DataFrame:
        """
        Re-analyze CAS9 trees using the main source metrics computation.
        """
        logger.info("Re-analyzing CAS9 trees...")

        # Load CAS9 trees
        self._load_cas9_trees(max_trees)

        all_results = []

        for cas9_key, cas9_tree in self.cas9_trees.items():
            logger.info(f"Processing: {cas9_key}")

            # Parse CAS9 filename: instance{id}_sim{sim}_tier{tier}_instance
            parts = cas9_key.split('_')
            instance_id = 0
            sim_id = 0
            tier_num = None

            for part in parts:
                if part.startswith('instance') and parts.index(part) == 0:
                    try:
                        instance_id = int(part.replace('instance', ''))
                    except:
                        pass
                elif part.startswith('sim'):
                    try:
                        sim_id = int(part.replace('sim', ''))
                    except:
                        pass
                elif part.startswith('tier'):
                    try:
                        tier_num = int(part.replace('tier', ''))
                    except:
                        pass

            if tier_num is None:
                logger.warning(f"Skipping {cas9_key}: couldn't parse tier")
                continue

            # Get tier information
            tier_info = self._get_tier_info(tier_num)

            # For CAS9 trees, create basic row (simplified for now)
            try:
                params = extract_parameters_from_tree(cas9_tree, "simulation")

                # Create reconstruction ID
                reconstruction_id = f"instance{instance_id}_sim{sim_id}_tier{tier_num}_cas9"

                # Create basic row for CAS9 tree
                row = {
                    'reconstruction_id': reconstruction_id,
                    'solver': 'cas9_recording',
                    'tier': tier_num,
                    'phs_lam_source': 'simulation',
                    'cPHS': np.nan,  # Not applicable without reconstruction
                    'RF_distance': np.nan,
                    'triplets_distance': np.nan,
                    **params
                }

                rows = [row]

                all_results.extend(rows)

                logger.debug(f"Added CAS9 tree: {reconstruction_id}")

            except Exception as e:
                logger.error(f"Failed to compute metrics for {cas9_key}: {e}")

        return pd.DataFrame(all_results)


def main():
    parser = argparse.ArgumentParser(
        description='Re-analyze experiment using main source code'
    )
    parser.add_argument('experiment_dir', type=str,
                        help='Path to experiment directory')
    parser.add_argument('--mode', choices=['reconstructed', 'cas9', 'both'],
                        default='both',
                        help='What to re-analyze')
    parser.add_argument('--output', type=str, default=None,
                        help='Output CSV file')
    parser.add_argument('--max-trees', type=int, default=None,
                        help='Maximum trees to process')
    parser.add_argument('--jobs', '-j', type=int, default=None,
                        help='Number of parallel jobs (default: auto-detect CPU count)')
    parser.add_argument('--sequential', action='store_true',
                        help='Force sequential processing (for debugging)')

    args = parser.parse_args()

    experiment_dir = Path(args.experiment_dir)
    if not experiment_dir.exists():
        logger.error(f"Experiment directory not found: {experiment_dir}")
        sys.exit(1)

    # Initialize reanalyzer
    reanalyzer = ExperimentReanalyzer(experiment_dir)

    results_df = pd.DataFrame()

    # Re-analyze based on mode
    if args.mode in ['reconstructed', 'both']:
        if args.sequential:
            logger.info("Using sequential processing")
            recon_results = reanalyzer.reanalyze_reconstructed_trees_sequential(args.max_trees)
        else:
            logger.info("Using parallel processing")
            recon_results = reanalyzer.reanalyze_reconstructed_trees(args.max_trees, args.jobs)
        results_df = pd.concat([results_df, recon_results], ignore_index=True)

    if args.mode in ['cas9', 'both']:
        cas9_results = reanalyzer.reanalyze_cas9_trees(args.max_trees)
        results_df = pd.concat([results_df, cas9_results], ignore_index=True)

    if not results_df.empty:
        logger.info(f"Re-analysis complete: {len(results_df)} rows generated")

        # Print summary
        if 'solver' in results_df.columns:
            logger.info("Results by solver:")
            for solver in sorted(results_df['solver'].unique()):
                solver_data = results_df[results_df['solver'] == solver]
                sim_data = solver_data[solver_data['phs_lam_source'] == 'simulation']
                if not sim_data.empty:
                    rf_mean = sim_data['RF_distance'].mean()
                    cphs_mean = sim_data['cPHS'].mean()
                    logger.info(f"  {solver}: {len(sim_data)} trees, "
                              f"RF={rf_mean:.3f}, cPHS={cphs_mean:.6f}")

        # Save results
        if args.output:
            results_df.to_csv(args.output, index=False)
            logger.info(f"Results saved to: {args.output}")
        else:
            # Default output name
            output_file = experiment_dir / f"reanalysis_results_{args.mode}.csv"
            results_df.to_csv(output_file, index=False)
            logger.info(f"Results saved to: {output_file}")

    else:
        logger.warning("No results generated")


if __name__ == "__main__":
    main()