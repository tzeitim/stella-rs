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


from core.shared_metrics import calculate_metrics_for_trees
import cassiopeia as cass

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


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

    def reanalyze_reconstructed_trees(self, max_trees: int = None) -> pd.DataFrame:
        """
        Re-analyze reconstructed trees using the main source metrics computation.
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
                logger.warning(f"No matching CAS9 tree found for {cas9_key}")

            # Get tier information
            tier_info = self._get_tier_info(parsed['tier_num'])

            # Compute all metrics using the shared pipeline implementation
            try:
                # Extract parameters based on mode
                lam_sim = lam_gt = q_sim = q_gt = None
                if cas9_tree and hasattr(cas9_tree, 'parameters'):
                    lam_sim = cas9_tree.parameters.get('lam_true')
                    q_sim = cas9_tree.parameters.get('q_true')
                    lam_gt = cas9_tree.parameters.get('lam_gt')
                    q_gt = cas9_tree.parameters.get('q_gt')

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
                    q_gt=q_gt
                )

                # Create reconstruction ID
                reconstruction_id = (f"instance{parsed['instance_id']}_sim{parsed['sim_id']}_"
                                   f"recon{parsed['reconstruction_id']}_tier{parsed['tier_num']}_"
                                   f"{parsed['solver']}")

                # Flatten the metrics for DataFrame creation
                row = {
                    'reconstruction_id': reconstruction_id,
                    'solver': parsed['solver'],
                    'tier': parsed['tier_num'],
                    **metrics
                }

                # Handle parsimony score special case
                if isinstance(metrics.get('parsimony_score'), dict):
                    row.update(metrics['parsimony_score'])
                    del row['parsimony_score']

                rows = [row]

                all_results.extend(rows)

                logger.debug(f"Computed metrics: RF={metrics.get('RF_distance', 'N/A'):.3f}, "
                           f"cPHS={metrics.get('phs_sim', 'N/A'):.6f}")

            except Exception as e:
                logger.error(f"Failed to compute metrics for {recon_file.name}: {e}")

        return pd.DataFrame(all_results)

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

            # Extract parameters and compute basic metrics for CAS9 tree
            try:
                params = extract_parameters_from_tree(cas9_tree, "simulation")

                # For CAS9 trees, we don't compute RF/triplets (no reconstruction)
                # but we can compute PHS using the tree itself
                from core.metrics_computation import calculate_phs_scores

                phs_scores = calculate_phs_scores(
                    cas9_tree,
                    params['lam_from_simulation'],
                    params['q_from_simulation'],
                    params.get('lam_gt'),
                    params.get('q_gt')
                )

                # Combine all metrics
                metrics = {**params, **phs_scores}
                metrics.update({
                    'RF_distance': np.nan,  # Not applicable for CAS9 trees
                    'triplets_distance': np.nan,  # Not applicable for CAS9 trees
                })

                # Create reconstruction ID
                reconstruction_id = f"instance{instance_id}_sim{sim_id}_tier{tier_num}_cas9"

                # Create standardized rows
                rows = create_metrics_rows(
                    base_metrics=metrics,
                    reconstruction_id=reconstruction_id,
                    solver='cas9_recording',
                    tier_info=tier_info
                )

                all_results.extend(rows)

                logger.debug(f"Computed CAS9 metrics: cPHS={metrics.get('phs_sim', 'N/A'):.6f}")

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
        recon_results = reanalyzer.reanalyze_reconstructed_trees(args.max_trees)
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