#!/usr/bin/env python3
"""
Local Test Pipeline for Fast Metrics Testing

This module provides a local, non-LSF version of the cascade pipeline
for rapid iteration and testing of metrics computation fixes.

Key features:
- Smart caching of GT trees and CAS9 instances
- Local function calls instead of job submissions
- Reuses production code paths for realistic testing
- Fast execution (1-2 minutes vs hours)

Usage:
    python local_test_pipeline.py --config test_config_mini.yaml
"""

import pickle
import hashlib
import logging
import time
import yaml
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Tuple
import sys
import os

# Add the parent directories to Python path for imports
current_dir = Path(__file__).parent
src_dir = current_dir.parent
cascade_dir = src_dir.parent
sys.path.extend([str(src_dir), str(cascade_dir)])

# Import production workers
from workers.master_gt_worker import generate_ground_truth_tree
from workers.cas9_recording_worker import Cas9RecordingWorker
from workers.reconstruction_worker import reconstruct_and_calculate_metrics

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class LocalTestPipeline:
    """
    Local test pipeline that bypasses LSF and runs everything locally.
    Includes smart caching to avoid regenerating expensive data.
    """

    def __init__(self, config_path: str, cache_dir: str = None):
        """
        Initialize the local test pipeline.

        Args:
            config_path: Path to test configuration YAML file
            cache_dir: Directory for caching generated data (default: ./cache)
        """
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.cache_dir = Path(cache_dir or current_dir / "cache")
        self.results_dir = current_dir / "results"

        # Ensure directories exist
        self.cache_dir.mkdir(exist_ok=True)
        self.results_dir.mkdir(exist_ok=True)

        logger.info(f"Initialized LocalTestPipeline with config: {config_path}")
        logger.info(f"Cache dir: {self.cache_dir}")
        logger.info(f"Results dir: {self.results_dir}")

    def _load_config(self) -> Dict[str, Any]:
        """Load and validate test configuration."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_path}")

        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)

        logger.info(f"Loaded config with {config.get('execution', {}).get('num_gt_instances', 0)} GT instances")
        return config

    def _config_hash(self, config_section: Dict[str, Any]) -> str:
        """Generate a hash for configuration caching."""
        config_str = yaml.dump(config_section, sort_keys=True)
        return hashlib.md5(config_str.encode()).hexdigest()[:8]

    def get_or_generate_gt_trees(self) -> List[Any]:
        """
        Get GT trees from cache or generate new ones.
        Uses production GTWorker for realistic generation.
        """
        gt_config = self.config['ground_truth']
        execution_config = self.config['execution']

        cache_key = self._config_hash({
            'ground_truth': gt_config,
            'random_seed': execution_config['random_seed']
        })
        cache_file = self.cache_dir / f"gt_trees_{cache_key}.pkl"

        if cache_file.exists():
            logger.info(f"📁 Reusing cached GT trees: {cache_file.name}")
            with open(cache_file, 'rb') as f:
                return pickle.load(f)

        logger.info("🌳 Generating new GT trees...")
        start_time = time.time()

        gt_trees = []
        num_instances = execution_config['num_gt_instances']
        seeds = gt_config.get('instance_seeds', [])

        # Ensure we have enough seeds
        if len(seeds) < num_instances:
            base_seed = execution_config['random_seed']
            seeds.extend([base_seed + i for i in range(len(seeds), num_instances)])

        for i in range(num_instances):
            logger.info(f"  Generating GT instance {i} with seed {seeds[i]}")

            # Set the seed for this instance
            import numpy as np
            np.random.seed(seeds[i])

            # Use production generate_ground_truth_tree function
            gt_tree = generate_ground_truth_tree(self.config)

            gt_trees.append({
                'instance_id': i,
                'seed': seeds[i],
                'tree': gt_tree
            })

        # Cache the results
        with open(cache_file, 'wb') as f:
            pickle.dump(gt_trees, f)

        generation_time = time.time() - start_time
        logger.info(f"✅ Generated {len(gt_trees)} GT trees in {generation_time:.1f}s")

        return gt_trees

    def get_or_generate_cas9_instances(self, gt_trees: List[Any]) -> List[Dict[str, Any]]:
        """
        Get CAS9 instances from cache or generate new ones.
        Uses production Cas9RecordingWorker.
        """
        cas9_config = {
            'cas9_tiers': self.config['cas9_tiers'],
            'cas9_simulations_per_gt': self.config['cas9_simulations_per_gt'],
            'execution': self.config['execution']
        }

        cache_key = self._config_hash(cas9_config)
        cache_file = self.cache_dir / f"cas9_instances_{cache_key}.pkl"

        if cache_file.exists():
            logger.info(f"📁 Reusing cached CAS9 instances: {cache_file.name}")
            with open(cache_file, 'rb') as f:
                return pickle.load(f)

        logger.info("🧬 Generating new CAS9 instances...")
        start_time = time.time()

        cas9_instances = []
        simulations_per_gt = self.config['cas9_simulations_per_gt']

        for gt_data in gt_trees:
            gt_tree = gt_data['tree']
            instance_id = gt_data['instance_id']

            for sim_id in range(simulations_per_gt):
                for tier_num, tier_config in self.config['cas9_tiers'].items():
                    logger.info(f"  CAS9 instance {instance_id}, sim {sim_id}, tier {tier_num}")

                    # Copy config to expected location for Cas9RecordingWorker
                    config_copy = self.results_dir / "cascade_config.yaml"
                    import shutil
                    shutil.copy2(self.config_path, config_copy)

                    # Use production Cas9RecordingWorker (need to create a temp file)
                    # Save GT tree temporarily
                    temp_gt_file = self.cache_dir / f"temp_gt_{instance_id}.pkl"
                    with open(temp_gt_file, 'wb') as f:
                        pickle.dump(gt_tree, f)

                    cas9_worker = Cas9RecordingWorker(
                        gt_tree_path=str(temp_gt_file),
                        tier=tier_num,
                        instance_id=instance_id,
                        cas9_simulation_id=sim_id,
                        output_dir=self.results_dir / "cas9_instances",
                        shared_dir=self.results_dir
                    )

                    # Import and create Cas9TierConfig object
                    sys.path.append(str(src_dir.parent))
                    from config_loader import Cas9TierConfig

                    # Convert dict config to tier config object
                    tier_config_obj = Cas9TierConfig(**tier_config)

                    # Generate mutation rates for the tier
                    tier_config_obj.mutation_rates = tier_config_obj.generate_mutation_rates()

                    # Apply CAS9 recording using the static method
                    from workers.cas9_recording_worker import apply_cas9_recording_to_tree
                    cas9_tree = apply_cas9_recording_to_tree(gt_tree, tier_config_obj, tier_num, self.config)

                    cas9_instances.append({
                        'gt_instance_id': instance_id,
                        'cas9_simulation_id': sim_id,
                        'tier_num': tier_num,
                        'tier_config': tier_config_obj,  # Save the constructed object, not the dict
                        'cas9_tree': cas9_tree,
                        'gt_tree': gt_tree
                    })

        # Cache the results
        with open(cache_file, 'wb') as f:
            pickle.dump(cas9_instances, f)

        generation_time = time.time() - start_time
        logger.info(f"✅ Generated {len(cas9_instances)} CAS9 instances in {generation_time:.1f}s")

        return cas9_instances

    def run_reconstructions(self, cas9_instances: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Run reconstructions on CAS9 instances using production reconstruction worker.
        This is where we test our metrics fixes!
        """
        logger.info("🔧 Running reconstructions...")
        start_time = time.time()

        results = []
        solvers = [s for s in self.config['solvers']['enabled']]

        total_reconstructions = len(cas9_instances) * len(solvers)
        current_recon = 0

        for cas9_data in cas9_instances:
            cas9_tree = cas9_data['cas9_tree']
            gt_instance_id = cas9_data['gt_instance_id']
            cas9_simulation_id = cas9_data['cas9_simulation_id']
            tier_num = cas9_data['tier_num']
            tier_config = cas9_data['tier_config']

            for solver_name in solvers:
                current_recon += 1
                logger.info(f"  [{current_recon}/{total_reconstructions}] "
                          f"GT{gt_instance_id}_sim{cas9_simulation_id}_tier{tier_num}_{solver_name}")

                try:
                    # Use production reconstruction worker
                    result = reconstruct_and_calculate_metrics(
                        cas9_tree=cas9_tree,
                        solver_name=solver_name,
                        tier_num=tier_num,
                        tier_config=tier_config,
                        gt_instance_id=gt_instance_id,
                        cas9_simulation_id=cas9_simulation_id,
                        reconstruction_id=0,  # Always 0 for testing
                        config=self.config
                    )

                    # Add metadata for analysis
                    result.update({
                        'test_run': True,
                        'config_path': str(self.config_path),
                        'timestamp': time.time()
                    })

                    # Apply production dual-row logic: create separate rows for each parameter source
                    # 1. Simulation parameters result (always created)
                    result_simulation = result.copy()
                    result_simulation['phs_lam_source'] = 'simulation'
                    result_simulation['phs_q_source'] = 'simulation'
                    results.append(result_simulation)

                    # 2. Ground truth parameters result (only if GT parameters are available)
                    if result.get('cPHS_gt') is not None and not np.isnan(result.get('cPHS_gt', np.nan)):
                        result_gt = result.copy()
                        result_gt['phs_lam_source'] = 'ground_truth'
                        result_gt['phs_q_source'] = 'ground_truth'
                        results.append(result_gt)
                        logger.info(f"    📈 Added ground truth row (cPHS_gt={result['cPHS_gt']:.6f})")

                    # Log key metrics to verify our fixes
                    likelihood_sim = result.get('likelihood_score_simulation', {})
                    likelihood_gt = result.get('likelihood_score_gt', {})

                    if isinstance(likelihood_sim, dict):
                        ll_sim = likelihood_sim.get('log_likelihood', 'N/A')
                        ll_gt = likelihood_gt.get('log_likelihood', 'N/A') if isinstance(likelihood_gt, dict) else 'N/A'
                        logger.info(f"    📊 Likelihood: sim={ll_sim}, gt={ll_gt}")

                except Exception as e:
                    logger.error(f"    ❌ Reconstruction failed: {e}")
                    results.append({
                        'gt_instance_id': gt_instance_id,
                        'cas9_simulation_id': cas9_simulation_id,
                        'tier_num': tier_num,
                        'solver': solver_name,
                        'error': str(e),
                        'test_run': True,
                        'timestamp': time.time()
                    })

        reconstruction_time = time.time() - start_time
        successful_results = len([r for r in results if 'error' not in r])
        logger.info(f"✅ Completed {successful_results}/{len(results)} reconstructions in {reconstruction_time:.1f}s")

        return results

    def save_results(self, results: List[Dict[str, Any]]) -> str:
        """Save results to CSV file with timestamp."""
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        filename = self.results_dir / f"test_results_{timestamp}.csv"

        # Convert results to flat format for CSV
        flattened_results = []
        for result in results:
            flat_result = self._flatten_result_dict(result)
            flattened_results.append(flat_result)

        import pandas as pd
        df = pd.DataFrame(flattened_results)
        df.to_csv(filename, index=False)

        logger.info(f"💾 Results saved to: {filename}")
        return str(filename)

    def _flatten_result_dict(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Flatten nested result dictionary for CSV export."""
        flat = {}
        for key, value in result.items():
            if isinstance(value, dict) and value:  # Only flatten non-empty dictionaries
                # Flatten nested dictionary
                for sub_key, sub_value in value.items():
                    flat[f"{key}_{sub_key}"] = sub_value
                # Don't include the original nested key to avoid empty columns
            elif not isinstance(value, dict):  # Only include non-dict values
                flat[key] = value
            # Skip empty dictionaries and nested dict keys entirely
        return flat

    def run_full_pipeline(self) -> str:
        """
        Run the complete test pipeline: GT generation → CAS9 recording → Reconstruction.
        Returns path to results file.
        """
        logger.info("🚀 Starting Local Test Pipeline")
        pipeline_start = time.time()

        # Phase 1: GT Trees
        gt_trees = self.get_or_generate_gt_trees()

        # Phase 2: CAS9 Instances
        cas9_instances = self.get_or_generate_cas9_instances(gt_trees)

        # Phase 3: Reconstructions
        results = self.run_reconstructions(cas9_instances)

        # Phase 4: Save Results
        results_file = self.save_results(results)

        total_time = time.time() - pipeline_start
        logger.info(f"🎉 Pipeline completed in {total_time:.1f}s")
        logger.info(f"📈 Total tests: {len(results)}")
        logger.info(f"📄 Results: {results_file}")

        return results_file


def main():
    """Command-line interface for the local test pipeline."""
    import argparse

    parser = argparse.ArgumentParser(description="Local Test Pipeline for Fast Metrics Testing")
    parser.add_argument("--config", required=True, help="Path to test configuration YAML")
    parser.add_argument("--cache-dir", help="Directory for caching (default: ./cache)")

    args = parser.parse_args()

    # Run the pipeline
    pipeline = LocalTestPipeline(args.config, args.cache_dir)
    results_file = pipeline.run_full_pipeline()

    print(f"\n✅ Test pipeline completed!")
    print(f"Results saved to: {results_file}")


if __name__ == "__main__":
    main()