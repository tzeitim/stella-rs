#!/usr/bin/env python3
"""
Master GT Worker - Level 1 of LSF Cascading Job System

Generates a ground truth tree and submits Cas9 recording jobs for all tiers.
This is the entry point that kicks off the entire analysis chain.
"""

import argparse
import json
import logging
import os
import pickle
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any
import numpy as np
import cassiopeia as cass
from cassiopeia.simulator import BirthDeathFitnessSimulator, UniformLeafSubsampler, Cas9LineageTracingDataSimulator
import yaml

# Add utils directory to path for job throttling
utils_dir = Path(__file__).parent.parent / "utils"
sys.path.append(str(utils_dir))
from job_throttling import create_throttler_from_config

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_config(config_path: str = "cascade_config.yaml") -> Dict[str, Any]:
    """Load configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.warning(f"Could not load config {config_path}: {e}")
        return {}

# CAS9 Tier Configuration
class Cas9SimulationTier:
    def __init__(self, name: str, recording_sites: int, states_per_site: int, 
                 mutation_rate: float, description: str = ""):
        self.name = name
        self.recording_sites = recording_sites
        self.states_per_site = states_per_site
        self.mutation_rate = mutation_rate
        self.description = description
    
    @property
    def total_sites(self) -> int:
        return self.recording_sites
    
    def get_effective_mutation_rate(self) -> float:
        return self.mutation_rate

CAS9_TIERS = {
    1: Cas9SimulationTier(
        name="Ultra High Fidelity",
        recording_sites=5000,
        states_per_site=50,
        mutation_rate=0.5,
        description="Ultra high fidelity recording with maximum sites and states"
    ),
    2: Cas9SimulationTier(
        name="High Fidelity", 
        recording_sites=1000,
        states_per_site=30,
        mutation_rate=0.4,
        description="High fidelity recording for optimal balance"
    ),
    3: Cas9SimulationTier(
        name="Medium Fidelity",
        recording_sites=200,
        states_per_site=20,
        mutation_rate=0.3,
        description="Medium fidelity for standard analysis"
    ),
    4: Cas9SimulationTier(
        name="Low Fidelity",
        recording_sites=15,
        states_per_site=10,
        mutation_rate=0.2,
        description="Low fidelity for constrained scenarios"
    )
}

def generate_state_priors(k, m, exp):
    """Generate state priors (q_i) for character states

    IMPORTANT: State 0 is reserved for unedited state, so mutation states start from 1
    """
    state_priors = np.array([np.random.exponential(exp) for _ in range(m)])
    state_priors /= np.sum(state_priors)
    # Fix: states should be 1, 2, 3, ..., m (not 0, 1, 2, ..., m-1)
    return {i+1: state_priors[i] for i in range(m)}

def generate_ground_truth_tree(config: Dict[str, Any]):
    """Generate a ground truth tree using Cassiopeia simulation"""
    logger.info("Simulating ground truth tree with birth-death model...")
    
    # Get tree configuration from config file
    gt_config = config.get('ground_truth', {}).get('tree_config', {})
    
    # Tree configuration with config values, falling back to defaults
    tree_config = {
        'N': gt_config.get('N', 1e3),  # Number of cells in original tree
        'n': gt_config.get('n', 1e2),  # Number of cells after subsampling  
        'fitness': {
            'birth_waiting_distribution': lambda scale: np.random.exponential(1/scale),
            'initial_birth_scale': gt_config.get('fitness', {}).get('initial_birth_scale', 2),
            'death_waiting_distribution': lambda: np.inf,  # Cells don't die
            'mutation_distribution': lambda: 1 if np.random.uniform() < 0.5 else 0,
            'fitness_distribution': lambda: np.random.normal(0.5, gt_config.get('fitness', {}).get('fitness_std', 0.25)),
            'fitness_base': gt_config.get('fitness', {}).get('fitness_base', 1.1)
        },
        'k': 50,  # Sequence length (can be made configurable later)
        'm': 50,  # Number of unique mutations (can be made configurable later)
        'rho': 0.5,  # Character mutation probability (can be made configurable later)
        'state_priors_exponents': 1e-5  # (can be made configurable later)
    }
    
    logger.info(f"Tree config - N: {tree_config['N']}, n: {tree_config['n']}")
    logger.info(f"Fitness config - base: {tree_config['fitness']['fitness_base']}, birth_scale: {tree_config['fitness']['initial_birth_scale']}")
    
    # Simulate the original tree topology
    topology_simulator = BirthDeathFitnessSimulator(
        **tree_config['fitness'], 
        num_extant=int(tree_config['N'])
    )
    original_topology = topology_simulator.simulate_tree()
    
    # Subsample the tree topology to get the GT tree
    leaf_subsampler = UniformLeafSubsampler(number_of_leaves=int(tree_config['n']))
    gt_tree = leaf_subsampler.subsample_leaves(original_topology)
    gt_tree.scale_to_unit_length()
    
    # Generate lineage tracing data sequences for the GT tree
    lam = -np.log(1.0 - tree_config["rho"])
    k, m, exp = tree_config['k'], tree_config['m'], tree_config['state_priors_exponents']
    priors = generate_state_priors(k, m, exp)
    
    lt_simulator = Cas9LineageTracingDataSimulator(
        number_of_cassettes=k,
        size_of_cassette=1,
        number_of_states=m,
        mutation_rate=lam,
        state_priors=priors
    )
    lt_simulator.overlay_data(gt_tree)
    
    # Update tree parameters
    gt_tree.priors = {i: priors for i in range(k)}
    gt_tree.parameters["stochastic_missing_rate"] = 0
    gt_tree.parameters["heritable_missing_rate"] = 0

    # Calculate GT lam and q for propagation to downstream workers
    cm_gt = gt_tree.character_matrix.to_numpy()
    proportion_mutated_gt = np.sum(cm_gt > 0) / np.sum(cm_gt >= 0)
    lam_gt = -np.log(1.0 - proportion_mutated_gt)

    # Calculate q from GT priors
    priors_array = np.array(list(priors.values()))
    q_gt = np.sum(priors_array ** 2)

    # Store these for downstream workers
    gt_tree.parameters["lam_gt"] = lam_gt
    gt_tree.parameters["q_gt"] = q_gt
    gt_tree.parameters["proportion_mutated_gt"] = proportion_mutated_gt

    logger.info(f"Generated GT tree with {len(gt_tree.leaves)} leaves")
    logger.info(f"GT parameters: lam={lam_gt:.4f}, q={q_gt:.4f}, prop_mutated={proportion_mutated_gt:.4f}")
    return gt_tree


class MasterGTWorker:
    """Master worker that generates GT tree and submits Cas9 recording jobs."""
    
    def __init__(self, output_dir: str, shared_dir: str):
        self.output_dir = Path(output_dir)
        self.shared_dir = Path(shared_dir)
        self.gt_tree_path = None

        # Load configuration
        config_path = self.shared_dir / "cascade_config.yaml"
        self.config = load_config(str(config_path))
        logger.info(f"Loaded config from {config_path}")

        # Initialize job throttler with dynamic config support
        self.throttler = create_throttler_from_config(self.config, self.shared_dir)

        # Ensure directories exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.shared_dir / "status").mkdir(parents=True, exist_ok=True)
        (self.shared_dir / "jobs").mkdir(parents=True, exist_ok=True)
        
    def generate_gt_trees(self) -> list:
        """Generate multiple ground truth trees based on num_gt_instances."""
        num_instances = self.config.get('execution', {}).get('num_gt_instances', 1)
        logger.info(f"Generating {num_instances} ground truth tree instances...")
        
        gt_tree_paths = []
        
        for instance_id in range(num_instances):
            logger.info(f"Generating GT instance {instance_id + 1}/{num_instances}...")
            
            try:
                gt_tree = generate_ground_truth_tree(self.config)
                
                # Save GT tree to configured directory
                config_shared_dir = self.config.get('output', {}).get('shared_dir', str(self.shared_dir))
                gt_output_dir = Path(config_shared_dir) / "gt_trees"
                gt_output_dir.mkdir(parents=True, exist_ok=True)
                
                # Create instance-specific filename
                gt_tree_path = gt_output_dir / f"gt_tree_instance_{instance_id}.pkl"
                
                # Save the GT tree
                with open(gt_tree_path, 'wb') as f:
                    pickle.dump(gt_tree, f)
                    
                logger.info(f"GT tree instance {instance_id} saved to: {gt_tree_path}")
                gt_tree_paths.append(gt_tree_path)
                
            except Exception as e:
                logger.error(f"Failed to generate GT instance {instance_id}: {e}")
                raise
        
        # Also save a copy of the first instance as gt_tree.pkl for compatibility
        if gt_tree_paths:
            config_shared_dir = self.config.get('output', {}).get('shared_dir', str(self.shared_dir))
            compatibility_path = Path(config_shared_dir) / "gt_tree.pkl" 
            with open(gt_tree_paths[0], 'rb') as src:
                with open(compatibility_path, 'wb') as dst:
                    dst.write(src.read())
            logger.info(f"Compatibility GT tree saved to: {compatibility_path}")
        
        # Update status
        # self.update_status("level1", "gt_generation", "completed")
        return gt_tree_paths
    
    def submit_cas9_jobs(self, gt_tree_paths: list) -> None:
        """Submit LSF jobs for all Cas9 tiers, GT instances, and simulations."""
        num_instances = len(gt_tree_paths)
        
        # Get cas9_simulations_per_gt from config
        cas9_simulations_per_gt = self.config.get('execution', {}).get('cas9_simulations_per_gt', 1)
        
        logger.info(f"Submitting Cas9 recording jobs for {num_instances} instances with {cas9_simulations_per_gt} simulations each...")
        
        submitted_jobs = []
        
        # Get tiers from config, default to all 4 if not specified
        config_tiers = self.config.get('cas9_tiers', {})
        if not config_tiers:
            logger.warning("No cas9_tiers found in config, using default tiers 1-4")
            tier_nums = [1, 2, 3, 4]
        else:
            tier_nums = list(config_tiers.keys())
            logger.info(f"Using tiers from config: {tier_nums}")
        
        # Submit jobs for each GT instance, each simulation, and each tier
        for instance_id, gt_tree_path in enumerate(gt_tree_paths):
            for cas9_simulation_id in range(cas9_simulations_per_gt):
                for tier_num in tier_nums:
                    try:
                        job_id = self.submit_cas9_job(tier_num, instance_id, cas9_simulation_id, str(gt_tree_path))
                        submitted_jobs.append({
                            'instance': instance_id,
                            'cas9_simulation_id': cas9_simulation_id,
                            'tier': tier_num,
                            'job_id': job_id,
                            'status': 'submitted',
                            'gt_tree_path': str(gt_tree_path)
                        })
                        logger.info(f"Submitted Cas9 Instance {instance_id} Sim {cas9_simulation_id} Tier {tier_num} job: {job_id}")
                        
                    except Exception as e:
                        logger.error(f"Failed to submit Cas9 Instance {instance_id} Sim {cas9_simulation_id} Tier {tier_num} job: {e}")
                        submitted_jobs.append({
                            'instance': instance_id,
                            'cas9_simulation_id': cas9_simulation_id,
                            'tier': tier_num,
                            'job_id': None,
                            'status': 'failed',
                            'error': str(e),
                            'gt_tree_path': str(gt_tree_path)
                        })
        
        # Update status with submitted jobs
        # self.update_status("level2", "cas9_jobs", "submitted", {
        #     'num_instances': num_instances,
        #     'submitted_jobs': submitted_jobs
        # })
        
    def submit_cas9_job(self, tier: int, instance_id: int, cas9_simulation_id: int, gt_tree_path: str) -> str:
        """Submit LSF job for a specific Cas9 tier, instance, and simulation with throttling."""

        def _do_submit():
            """Internal function to perform the actual submission"""
            job_script = self.shared_dir / "jobs" / f"cas9_tier{tier}_job.lsf"

            # Use the actual shared directory that was passed to this worker
            # Don't use config's shared_dir as it may not match the actual directory structure
            output_base = self.shared_dir

            # Build bsub command with instance and simulation-specific naming
            cmd = [
                'bsub',
                '-J', f'cas9_instance{instance_id}_sim{cas9_simulation_id}_tier{tier}_analysis',
                '-oo', f"{self.shared_dir.resolve()}/logs/cas9_instance{instance_id}_sim{cas9_simulation_id}_tier{tier}_%J.out",
                '-eo', f"{self.shared_dir.resolve()}/logs/cas9_instance{instance_id}_sim{cas9_simulation_id}_tier{tier}_%J.err",
                '-W', '0:30',  # 30 minutes
                '-n', '15', '-R', 'span[hosts=1]',
                '-R', 'rusage[mem=1.5GB]',
                'python', str(Path(__file__).parent / 'cas9_recording_worker.py'),
                '--gt_tree_path', gt_tree_path,
                '--tier', str(tier),
                '--instance', str(instance_id),
                '--cas9_simulation_id', str(cas9_simulation_id),
                '--output_dir', str(output_base / "cas9_instances"),
                '--shared_dir', str(output_base)
            ]

            # Submit job
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                raise RuntimeError(f"bsub failed: {result.stderr}")

            # Extract job ID from bsub output
            job_id = result.stdout.strip().split('<')[1].split('>')[0]
            logger.info(f"Submitted CAS9 job {job_id} for instance {instance_id}, sim {cas9_simulation_id}, tier {tier}")
            return job_id

        # Use throttling to submit the job
        return self.throttler.submit_with_throttling(_do_submit, 'cas9')
    
    def update_status(self, level: str, component: str, status: str, details: Any = None) -> None:
        """Update job status in shared status file."""
        status_file = self.shared_dir / "status" / f"{level}_status.json"
        
        # Load existing status
        if status_file.exists():
            with open(status_file, 'r') as f:
                status_data = json.load(f)
        else:
            status_data = {}
        
        # Update status
        status_data[component] = {
            'status': status,
            'timestamp': str(Path(__file__).stat().st_mtime),
            'details': details
        }
        
        # Save updated status
        with open(status_file, 'w') as f:
            json.dump(status_data, f, indent=2)
    
    def run(self) -> None:
        """Execute the master GT workflow."""
        logger.info("=== Master GT Worker Starting ===")
        
        try:
            # Step 1: Generate GT trees (multiple instances)
            gt_tree_paths = self.generate_gt_trees()
            
            # Step 2: Submit Cas9 recording jobs for all instances
            self.submit_cas9_jobs(gt_tree_paths)
            
            logger.info("=== Master GT Worker Completed Successfully ===")
            
        except Exception as e:
            logger.error(f"Master GT Worker failed: {e}")
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Master GT Worker - Generate GT tree and submit Cas9 jobs")
    parser.add_argument('--output_dir', required=True, help='Directory to save GT tree')
    parser.add_argument('--shared_dir', required=True, 
                       help='Shared directory for job coordination')
    
    args = parser.parse_args()
    
    worker = MasterGTWorker(args.output_dir, args.shared_dir)
    worker.run()


if __name__ == "__main__":
    main()