#!/usr/bin/env python3
"""
Configuration loader for the Cascading LSF Job System

Loads and validates configuration from YAML files, with support for:
- Multiple GT instances with different seeds
- Flexible solver selection
- Parameterized resource allocation
- Test mode for development
"""

import yaml
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)

@dataclass 
class LSFResource:
    cores: int
    memory_gb: float
    time_limit: str
    queue: str = "short"

@dataclass
class Cas9TierConfig:
    name: str
    description: str
    k: int  # number of cassettes
    cassette_size: int  # sites per cassette  
    m: int  # mutations per site
    mutation_rate_pattern: str
    base_mutation_rate: float
    state_priors_exponent: float
    
    @property
    def total_sites(self) -> int:
        return self.k * self.cassette_size
    
    def generate_mutation_rates(self) -> List[float]:
        """Generate mutation rates based on the specified pattern."""
        total_sites = self.total_sites
        
        if self.mutation_rate_pattern == "exponential_decay":
            # Exponential decay pattern
            rates = [self.base_mutation_rate * np.exp(-0.1 * i) for i in range(total_sites)]
        elif self.mutation_rate_pattern == "sine_wave":
            # Sinusoidal variation
            rates = [self.base_mutation_rate * (1.0 + 0.5 * np.sin(0.2 * i)) for i in range(total_sites)]
        elif self.mutation_rate_pattern == "linear":
            # Linear variation pattern
            rates = [self.base_mutation_rate + 0.3 * (i % 3 - 1) for i in range(total_sites)]
        elif self.mutation_rate_pattern == "uniform":
            # All sites have same rate
            rates = [self.base_mutation_rate] * total_sites
        else:
            raise ValueError(f"Unknown mutation rate pattern: {self.mutation_rate_pattern}")
            
        return rates

class CascadeConfig:
    """Main configuration class for the cascading job system."""
    
    def __init__(self, config_path: str = "cascade_config.yaml"):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self._setup_seeds()
        self._validate_config()
    
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from YAML file."""
        try:
            with open(self.config_path, 'r') as f:
                config = yaml.safe_load(f)
            logger.info(f"Loaded configuration from {self.config_path}")
            return config
        except Exception as e:
            logger.error(f"Failed to load config from {self.config_path}: {e}")
            raise
    
    def _setup_seeds(self):
        """Setup seeds for reproducible random generation."""
        master_seed = self.config.get('execution', {}).get('random_seed', 42)
        num_instances = self.num_gt_instances
        
        # Generate instance seeds if not provided
        instance_seeds = self.config.get('ground_truth', {}).get('instance_seeds', [])
        if len(instance_seeds) < num_instances:
            np.random.seed(master_seed)
            missing_seeds = num_instances - len(instance_seeds)
            additional_seeds = np.random.randint(1, 10000, missing_seeds).tolist()
            instance_seeds.extend(additional_seeds)
            
        self.instance_seeds = instance_seeds[:num_instances]
        logger.info(f"Using seeds: {self.instance_seeds}")
    
    def _validate_config(self):
        """Validate configuration parameters."""
        # Check that all required sections exist
        required_sections = ['execution', 'cas9_tiers', 'solvers', 'lsf']
        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Missing required config section: {section}")
        
        # Validate solver selection
        enabled_solvers = self.enabled_solvers
        if not enabled_solvers:
            raise ValueError("At least one solver must be enabled")
        
        # Validate tier configuration
        if not self.cas9_tiers:
            raise ValueError("At least one Cas9 tier must be configured")
        
        logger.info("Configuration validation passed")
    
    @property
    def num_gt_instances(self) -> int:
        """Number of ground truth instances to generate."""
        return self.config['execution'].get('num_gt_instances', 1)
    
    @property
    def enabled_solvers(self) -> List[str]:
        """List of enabled solvers."""
        return self.config['solvers'].get('enabled', ['nj', 'greedy'])
    
    @property 
    def reconstructions_per_solver(self) -> int:
        """Number of reconstructions per solver per Cas9 simulation."""
        return self.config['solvers'].get('reconstructions_per_solver', 1)
    
    @property
    def cas9_simulations_per_gt(self) -> int:
        """Number of Cas9 simulations to generate per ground truth."""
        return self.config['execution'].get('cas9_simulations_per_gt', 1)
    
    @property
    def cas9_tiers(self) -> Dict[int, Cas9TierConfig]:
        """Dictionary of Cas9 tier configurations."""
        tiers = {}
        for tier_id, tier_data in self.config['cas9_tiers'].items():
            tiers[int(tier_id)] = Cas9TierConfig(**tier_data)
        return tiers
    
    @property
    def lsf_resources(self) -> Dict[str, LSFResource]:
        """LSF resource configurations for different job types."""
        resources = {}
        lsf_config = self.config['lsf']
        
        for job_type, resource_data in lsf_config['resources'].items():
            queue = lsf_config['queues'].get(job_type, 'short')
            resources[job_type] = LSFResource(
                cores=resource_data['cores'],
                memory_gb=resource_data['memory_gb'],
                time_limit=resource_data['time_limit'],
                queue=queue
            )
        return resources
    
    @property
    def shared_dir(self) -> str:
        """Shared directory path."""
        val =  self.config['output'].get('shared_dir', None)
        if val is None:
            raise ValueError("a shared_dir config in the output key must be specified")
        return val
    
    @property
    def log_dir(self) -> str:
        """Log directory path.""" 
        return self.config['lsf'].get('log_dir', f"{self.shared_dir}/logs")
    
    @property
    def test_mode(self) -> bool:
        """Whether to run in test mode with smaller parameters."""
        return self.config.get('debug', {}).get('test_mode', False)
    
    def get_test_params(self) -> Dict[str, Any]:
        """Get test mode parameters."""
        return self.config.get('debug', {}).get('test_params', {})
    
    def total_expected_jobs(self) -> int:
        """Calculate total number of expected jobs with multiple Cas9 simulations per GT."""
        num_gt = self.num_gt_instances
        num_cas9_per_gt = self.cas9_simulations_per_gt
        num_tiers = len(self.cas9_tiers)
        num_solvers = len(self.enabled_solvers)
        reconstructions_per = self.reconstructions_per_solver
        
        # GT jobs + Cas9 jobs + Reconstruction jobs
        gt_jobs = num_gt
        cas9_jobs = num_gt * num_cas9_per_gt * num_tiers
        reconstruction_jobs = num_gt * num_cas9_per_gt * num_tiers * num_solvers * reconstructions_per
        
        return gt_jobs + cas9_jobs + reconstruction_jobs
    
    def get_instance_config(self, instance_id: int) -> Dict[str, Any]:
        """Get configuration for a specific GT instance."""
        if instance_id >= self.num_gt_instances:
            raise ValueError(f"Instance {instance_id} exceeds configured instances ({self.num_gt_instances})")
            
        return {
            'instance_id': instance_id,
            'seed': self.instance_seeds[instance_id],
            'gt_path': f"{self.shared_dir}/gt_trees/gt_instance_{instance_id:03d}.pkl",
            'cas9_dir': f"{self.shared_dir}/cas9_instances/instance_{instance_id:03d}",
            'results_dir': f"{self.shared_dir}/results/instance_{instance_id:03d}"
        }

def load_config(config_path: str = "cascade_config.yaml") -> CascadeConfig:
    """Load cascade configuration from file."""
    return CascadeConfig(config_path)

if __name__ == "__main__":
    # Test configuration loading
    config = load_config()
    
    print("Configuration Summary:")
    print(f"GT Instances: {config.num_gt_instances}")
    print(f"Instance Seeds: {config.instance_seeds}")
    print(f"Enabled Solvers: {config.enabled_solvers}")
    print(f"Cas9 Tiers: {list(config.cas9_tiers.keys())}")
    print(f"Total Expected Jobs: {config.total_expected_jobs()}")
    print(f"Test Mode: {config.test_mode}")
    
    # Show tier details
    for tier_id, tier_config in config.cas9_tiers.items():
        print(f"\nTier {tier_id}: {tier_config.name}")
        print(f"  Sites: {tier_config.total_sites} ({tier_config.k}×{tier_config.cassette_size})")
        print(f"  Pattern: {tier_config.mutation_rate_pattern}")
        rates = tier_config.generate_mutation_rates()
        print(f"  Rate range: {min(rates):.3f} - {max(rates):.3f}")
