#!/usr/bin/env python3
"""
Cas9 Recording Worker - Level 2 of LSF Cascading Job System

Applies Cas9 recording to a GT tree for a specific tier and submits reconstruction jobs
for all solvers. This is the middle layer that bridges GT generation and reconstruction.
"""

import argparse
import json
import logging
import os
import pickle
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any, List
import numpy as np
import cassiopeia as cass
from dataclasses import dataclass
from copy import deepcopy
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
@dataclass
class Cas9SimulationTier:
    """Configuration for a Cas9 simulation tier."""
    
    name: str
    description: str
    k: int  # number of cassettes (integrations)
    cassette_size: int  # number of sites per cassette
    m: int  # number of unique mutations per site
    mutation_rates: list  # mutation rates for different sites/cassettes
    state_priors_exponent: float
    
    @property
    def total_sites(self) -> int:
        """Total recording sites across all cassettes."""
        return self.k * self.cassette_size
    
    def get_effective_mutation_rate(self) -> float:
        return self.mutation_rates[0] if self.mutation_rates else 0.1
#   161         number_of_characters = size_of_cassette * number_of_cassettes                                                                                     
CAS9_TIERS = {
    1: Cas9SimulationTier(
        name="Tier 1 - Ultra High Fidelity",
        description="Near-perfect recording: 100 integrations × 50 sites with time-scaled mutation rates",
        k=100,  # 100 integrations
        cassette_size=50,  # 50 sites per integration
        m=50,  # 50 unique mutations per site
        # Time-scaled mutation rates: early sites mutate frequently, later sites rarely
        mutation_rates=list(np.hstack([2.0 * np.exp(-0.1 * i) for i in range(50)]* 100 )),  # Exponential decay
        state_priors_exponent=1e-6  # Very uniform priors
    ),
    
    2: Cas9SimulationTier(
        name="Tier 2 - High Fidelity", 
        description="Good recording: 50 integrations × 20 sites with varied mutation rates",
        k=50,   # 50 integrations
        cassette_size=20,  # 20 sites per integration  
        m=30,   # 30 unique mutations per site
        # Varied mutation rates across integrations
        mutation_rates=list(np.hstack([1.5 * (1.0 + 0.5 * np.sin(0.2 * i)) for i in range(20)] * 50 )),
        state_priors_exponent=1e-5
    ),
    
    3: Cas9SimulationTier(
        name="Tier 3 - Medium Fidelity",
        description="Moderate recording: 20 integrations × 10 sites with some variation",
        k=20,   # 20 integrations
        cassette_size=10,  # 10 sites per integration
        m=20,   # 20 unique mutations per site
        # Some variation in mutation rates
        mutation_rates=list(np.hstack([1.0 + 0.3 * (i % 3 - 1) for i in range(10)] * 20)),  # 0.7, 1.0, 1.3 pattern
        state_priors_exponent=1e-4
    ),
    
    4: Cas9SimulationTier(
        name="Tier 4 - Low Fidelity",
        description="Poor recording: 5 integrations × 3 sites with uniform mutation rates",
        k=5,    # 5 integrations
        cassette_size=3,   # 3 sites per integration
        m=10,   # 10 unique mutations per site
        # Uniform mutation rates (worst case)
        mutation_rates=list(np.hstack([0.8 for _ in range(3)] * 5)),  # All the same
        state_priors_exponent=1e-3
    )
}

# Solver definitions
# SOLVERS defined in solver_config.py

def apply_cas9_recording_to_tree(gt_tree, tier: Cas9SimulationTier, tier_number: int, config: dict = None):
    """Apply Cas9 recording simulation to ground truth tree for specified tier"""
    logger.info(f"Applying Cas9 recording for {tier.name} tier")
    
    # Create a copy of the tree to avoid modifying original
    tree_with_cas9 = deepcopy(gt_tree)
    
    # Generate new Cas9 recording with tier-specific parameters
    from cassiopeia.simulator import Cas9LineageTracingDataSimulator
    
    # Generate state priors for this tier
    def generate_state_priors(k, m, exp=1e-5):
        """Generate state priors (q_i) for character states

        IMPORTANT: State 0 is reserved for unedited state, so mutation states start from 1
        """
        state_priors = np.array([np.random.exponential(exp) for _ in range(m)])
        state_priors /= np.sum(state_priors)
        # Fix: states should be 1, 2, 3, ..., m (not 0, 1, 2, ..., m-1)
        return {i+1: state_priors[i] for i in range(m)}
    
    # Set recording parameters based on tier
    k = tier.total_sites
    m = tier.m
    #lam = [-np.log(1.0 - i) for i in tier.mutation_rates]
    print(f"first 10 {tier.mutation_rates[1:10]=}")
    print(f"{len(tier.mutation_rates)=}")
    print(f"{tier.cassette_size=}")

    lam = tier.mutation_rates
    priors = generate_state_priors(k, m)
    
    # Get tier-specific missing data parameters from config
    tier_config_dict = config.get('cas9_tiers', {}).get(tier_number, {}) if config else {}
    heritable_silencing_rate = tier_config_dict.get('heritable_silencing_rate', 0)
    stochastic_silencing_rate = tier_config_dict.get('stochastic_silencing_rate', 0)

    # Warn if non-zero missing rates are specified (should be 0 to match simulation_phs.py golden standard)
    if heritable_silencing_rate != 0 or stochastic_silencing_rate != 0:
        logger.warning(f"WARNING: Non-zero missing data rates specified for Tier {tier_number}. "
                      f"heritable_silencing_rate={heritable_silencing_rate}, "
                      f"stochastic_silencing_rate={stochastic_silencing_rate}. "
                      f"PHS not supported with missing data")

    logger.info(f"Missing data rates for Tier {tier_number}: "
                f"heritable={heritable_silencing_rate}, stochastic={stochastic_silencing_rate}")

    # Apply Cas9 lineage tracing simulation with tier-specific missing data controls
    lt_simulator = Cas9LineageTracingDataSimulator(
        number_of_cassettes=k,
        size_of_cassette=1,
        number_of_states=m,
        mutation_rate=lam,
        state_priors=priors,
        heritable_silencing_rate=heritable_silencing_rate,
        stochastic_silencing_rate=stochastic_silencing_rate,
        heritable_missing_data_state=-1,  # Standard missing data encoding
        stochastic_missing_data_state=-1  # Standard missing data encoding
    )
    
    # Overlay the Cas9 data on the tree
    lt_simulator.overlay_data(tree_with_cas9)
    
    # Update tree parameters with actual missing data rates
    tree_with_cas9.priors = {i: priors for i in range(k)}
    tree_with_cas9.parameters["stochastic_missing_rate"] = stochastic_silencing_rate
    tree_with_cas9.parameters["heritable_missing_rate"] = heritable_silencing_rate

    # Propagate GT parameters from original tree for downstream comparison
    if hasattr(gt_tree, 'parameters'):
        for param in ["lam_gt", "q_gt", "proportion_mutated_gt"]:
            if param in gt_tree.parameters:
                tree_with_cas9.parameters[param] = gt_tree.parameters[param]

    logger.info(f"Cas9 recording applied: {k} sites, {m} states per site")
    return tree_with_cas9


class Cas9RecordingWorker:
    """Worker that applies Cas9 recording and submits reconstruction jobs."""
    
    def __init__(self, gt_tree_path: str, tier: int, instance_id: int, output_dir: str, shared_dir: str, 
                 cas9_simulation_id: int = 0):
        self.gt_tree_path = Path(gt_tree_path)
        self.tier = tier
        self.instance_id = instance_id
        self.cas9_simulation_id = cas9_simulation_id
        self.output_dir = Path(output_dir)
        self.shared_dir = Path(shared_dir)
        self.cas9_instance_path = None
        
        # Load configuration
        config_path = self.shared_dir / "cascade_config.yaml"
        self.config = load_config(str(config_path))
        logger.info(f"Loaded config from {config_path}")

        # Initialize job throttler with dynamic config support
        self.throttler = create_throttler_from_config(self.config, self.shared_dir)

        # Ensure directories exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.shared_dir / "logs").mkdir(parents=True, exist_ok=True)

        # Validate tier against config
        config_tiers = list(self.config.get('cas9_tiers', {}).keys())
        if tier not in config_tiers:
            raise ValueError(f"Invalid tier: {tier}. Must be one of {config_tiers}")
    
    def load_gt_tree(self):
        """Load the ground truth tree."""
        logger.info(f"Loading GT tree from: {self.gt_tree_path}")
        
        if not self.gt_tree_path.exists():
            raise FileNotFoundError(f"GT tree file not found: {self.gt_tree_path}")
        
        with open(self.gt_tree_path, 'rb') as f:
            gt_tree = pickle.load(f)
        
        return gt_tree
    
    def apply_cas9_recording(self):
        """Apply Cas9 recording to the GT tree for this tier."""
        logger.info(f"Applying Cas9 recording for Tier {self.tier}...")
        
        try:
            # Load GT tree
            gt_tree = self.load_gt_tree()
            
            # Get tier configuration from config
            tier_config_dict = self.config.get('cas9_tiers', {}).get(self.tier, {})
            if not tier_config_dict:
                raise ValueError(f"Tier {self.tier} configuration not found")

            # Create tier config object from config data
            from config_loader import Cas9TierConfig
            tier_config = Cas9TierConfig(**tier_config_dict)

            logger.info(f"Tier config: {tier_config.name}")
            logger.info(f"Recording sites: {tier_config.k * tier_config.cassette_size}")
            
            # Generate mutation rates for the tier
            tier_config.mutation_rates = tier_config.generate_mutation_rates()

            # Apply Cas9 recording with config for missing data parameters
            cas9_tree = apply_cas9_recording_to_tree(gt_tree, tier_config, self.tier, self.config)
            
            # Save Cas9 instance to configured directory with proper indexing
            config_shared_dir = self.config.get('output', {}).get('shared_dir', str(self.shared_dir))
            cas9_instances_dir = Path(config_shared_dir) / "cas9_instances"
            cas9_instances_dir.mkdir(parents=True, exist_ok=True)
            self.cas9_instance_path = cas9_instances_dir / f"instance{self.instance_id}_sim{self.cas9_simulation_id}_tier{self.tier}_instance.pkl"
            with open(self.cas9_instance_path, 'wb') as f:
                pickle.dump(cas9_tree, f)
                
            logger.info(f"Cas9 instance saved to: {self.cas9_instance_path}")
            
            # Update status
            # self.update_status("level2", f"instance{self.instance_id}_tier{self.tier}_recording", "completed", {
            #     'tier_name': tier_config.name,
            #     'recording_sites': tier_config.k * tier_config.cassette_size,
            #     'cas9_instance_path': str(self.cas9_instance_path)
            # })
            
            return cas9_tree
            
        except Exception as e:
            logger.error(f"Failed to apply Cas9 recording for Tier {self.tier}: {e}")
            # self.update_status("level2", f"instance{self.instance_id}_tier{self.tier}_recording", "failed", str(e))
            raise
    
    def submit_reconstruction_jobs(self) -> None:
        """Submit LSF jobs for all solvers for this tier."""
        logger.info(f"Submitting reconstruction jobs for Tier {self.tier} across all solvers...")
        
        submitted_jobs = []
        
        # Get enabled solvers from config
        config_solvers = self.config.get('solvers', {}).get('enabled', ['nj', 'greedy'])
        logger.info(f"Using solvers from config: {config_solvers}")
        
        for solver in config_solvers:
            try:
                job_id = self.submit_reconstruction_job(solver)
                submitted_jobs.append({
                    'solver': solver,
                    'job_id': job_id,
                    'status': 'submitted'
                })
                logger.info(f"Submitted Tier {self.tier} {solver} reconstruction job: {job_id}")
                
            except Exception as e:
                logger.error(f"Failed to submit Tier {self.tier} {solver} job: {e}")
                submitted_jobs.append({
                    'solver': solver,
                    'job_id': None,
                    'status': 'failed',
                    'error': str(e)
                })
        
        # Update status with submitted jobs
        # self.update_status("level3", f"instance{self.instance_id}_tier{self.tier}_reconstructions", "submitted", {
        #     'cas9_instance_path': str(self.cas9_instance_path),
        #     'tier': self.tier,
        #     'submitted_jobs': submitted_jobs
        # })
    
    def submit_reconstruction_job(self, solver: str) -> str:
        """Submit LSF job for a specific solver reconstruction with throttling."""

        def _do_submit(recon_id: int):
            """Internal function to perform the actual submission"""
            # Get configured output directory or use default
            config_shared_dir = self.config.get('output', {}).get('shared_dir', str(self.shared_dir))
            output_base = Path(config_shared_dir)

            # Get LSF configuration for reconstruction jobs
            lsf_config = self.config.get('lsf', {})
            reconstruction_queue = lsf_config.get('queues', {}).get('reconstruction', 'long')
            recon_resources = lsf_config.get('resources', {}).get('reconstruction', {})

            cores = recon_resources.get('cores', 65)
            memory_gb = recon_resources.get('memory_gb', 1.5)

            # Build bsub command with config-based resources (let queue handle time limits)
            cmd = [
                'bsub',
                '-J', f'reconstruct_instance{self.instance_id}_sim{self.cas9_simulation_id}_recon{recon_id}_tier{self.tier}_{solver}',
                '-oo', f"{self.shared_dir.resolve()}/logs/reconstruct_instance{self.instance_id}_sim{self.cas9_simulation_id}_recon{recon_id}_tier{self.tier}_{solver}_%J.out",
                '-eo', f"{self.shared_dir.resolve()}/logs/reconstruct_instance{self.instance_id}_sim{self.cas9_simulation_id}_recon{recon_id}_tier{self.tier}_{solver}_%J.err",
                '-q', reconstruction_queue,                # Use configured queue (let queue set time limits)
                '-n', str(cores), '-R', 'span[hosts=1]',  # Use configured cores
                '-R', f'rusage[mem={memory_gb}GB]',        # Use configured memory
                'python', str(Path(__file__).parent / 'reconstruction_worker.py'),
                '--cas9_instance_path', str(self.cas9_instance_path),
                '--solver', solver,
                '--tier', str(self.tier),
                '--output_dir', str(output_base / "results"),
                '--shared_dir', str(output_base),
                '--gt_instance_id', str(self.instance_id),
                '--cas9_simulation_id', str(self.cas9_simulation_id),
                '--reconstruction_id', str(recon_id)
            ]

            # Submit job
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                raise RuntimeError(f"bsub failed for reconstruction {recon_id}: {result.stderr}")

            # Extract job ID from bsub output
            job_id = result.stdout.strip().split('<')[1].split('>')[0]
            logger.info(f"Submitted reconstruction job {job_id} for solver {solver}, recon {recon_id}")
            return job_id

        # Get reconstructions_per_solver from config
        reconstructions_per_solver = self.config.get('solvers', {}).get('reconstructions_per_solver', 1)

        # Submit multiple reconstruction jobs if reconstructions_per_solver > 1
        submitted_job_ids = []

        for recon_id in range(reconstructions_per_solver):
            # Use throttling to submit the job
            job_id = self.throttler.submit_with_throttling(_do_submit, 'reconstruction', recon_id)
            if job_id:
                submitted_job_ids.append(job_id)

        return ','.join(submitted_job_ids)  # Return all job IDs
    
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
        """Execute the Cas9 recording workflow."""
        logger.info(f"=== Cas9 Recording Worker Starting - Tier {self.tier} ===")
        
        try:
            # Step 1: Apply Cas9 recording
            self.apply_cas9_recording()
            
            # Step 2: Submit reconstruction jobs
            self.submit_reconstruction_jobs()
            
            logger.info(f"=== Cas9 Recording Worker Completed Successfully - Tier {self.tier} ===")
            
        except Exception as e:
            logger.error(f"Cas9 Recording Worker failed for Tier {self.tier}: {e}")
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Cas9 Recording Worker - Apply Cas9 recording and submit reconstruction jobs")
    parser.add_argument('--gt_tree_path', required=True, help='Path to ground truth tree')
    parser.add_argument('--tier', type=int, required=True, help='Cas9 recording tier (1-4)')
    parser.add_argument('--instance', type=int, required=True, help='GT instance ID')
    parser.add_argument('--cas9_simulation_id', type=int, default=0, help='Cas9 simulation ID for this GT instance')
    parser.add_argument('--output_dir', required=True, help='Directory to save Cas9 instances')
    parser.add_argument('--shared_dir', required=True, 
                       help='Shared directory for job coordination')
    
    args = parser.parse_args()
    
    worker = Cas9RecordingWorker(
        args.gt_tree_path, 
        args.tier, 
        args.instance,
        args.output_dir, 
        args.shared_dir,
        args.cas9_simulation_id
    )
    worker.run()


if __name__ == "__main__":
    main()
