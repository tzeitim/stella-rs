#!/usr/bin/env python3
"""
LSF Coordinator for Multi-Tier Cas9 Reconstruction Analysis

This script orchestrates the distributed execution of the multi-tier Cas9 analysis
across an LSF cluster. It generates job parameters, submits jobs with proper 
dependencies, monitors progress, and aggregates results.
"""

import os
import json
import subprocess
import time
import argparse
import pandas as pd
import yaml
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional
import logging

# Import our tier configurations
from reconstruction_metrics_table_multi_tier import CAS9_TIERS


@dataclass
class JobParameters:
    """Parameters for a single distributed job."""
    job_id: int
    job_type: str  # 'gt_generation', 'cas9_recording', 'reconstruction'
    gt_id: Optional[int] = None
    cas9_tier: Optional[int] = None
    solver: Optional[str] = None
    output_dir: str = ""
    depends_on: List[int] = None
    
    def __post_init__(self):
        if self.depends_on is None:
            self.depends_on = []


class LSFCoordinator:
    """Coordinates distributed execution of multi-tier Cas9 analysis."""

    def __init__(self,
                 n_gt_trees: int = 50,
                 base_dir: str = "/shared/cas9_analysis",
                 queue: str = "normal",
                 debug: bool = False,
                 config_file: str = None):

        self.n_gt_trees = n_gt_trees
        self.base_dir = Path(base_dir)
        self.queue = queue
        self.debug = debug

        # Setup basic logging first
        self.logger = logging.getLogger(__name__)

        # Load configuration from YAML file
        self.config = self.load_config(config_file)

        # Available solvers from config or default
        self.solvers = self.config.get('solvers', {}).get('enabled', [
            'nj', 'maxcut', 'maxcut_greedy', 'greedy',
            'smj', 'spectral', 'spectral_greedy'
        ])

        # Setup directories
        self.setup_directories()

        # Setup full logging after directories are created
        self.setup_logging()

        # Job tracking
        self.submitted_jobs = {}  # job_id -> LSF job ID
        self.job_parameters = {}  # job_id -> JobParameters

    def load_config(self, config_file: str = None) -> Dict:
        """Load configuration from YAML file."""
        if config_file is None:
            # Try to find config file in common locations
            potential_configs = [
                "cascade_config.yaml",
                "config.yaml",
                "ricosino/cascade_config.yaml",
                "cascade_workflow/configs/cascade_config.yaml"
            ]

            for config_path in potential_configs:
                if os.path.exists(config_path):
                    config_file = config_path
                    break

            if config_file is None:
                self.logger.warning("No config file found, using defaults")
                return {}

        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
                self.logger.info(f"Loaded configuration from: {config_file}")
                return config
        except Exception as e:
            self.logger.warning(f"Could not load config {config_file}: {e}, using defaults")
            return {}

    def setup_directories(self):
        """Create necessary directory structure."""
        directories = [
            "parameters",
            "gt_trees", 
            "cas9_instances",
            "results",
            "aggregated",
            "logs",
            "scripts"
        ]
        
        for dir_name in directories:
            dir_path = self.base_dir / dir_name
            dir_path.mkdir(parents=True, exist_ok=True)
            
        print(f"Created directory structure in {self.base_dir}")
    
    def setup_logging(self):
        """Setup logging for the coordinator."""
        log_file = self.base_dir / "logs" / "coordinator.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def generate_job_parameters(self) -> Dict[int, JobParameters]:
        """Generate parameters for all jobs in the distributed workflow."""
        
        job_id = 1
        job_params = {}
        
        # Phase 1: GT Tree Generation Jobs
        self.logger.info(f"Generating parameters for {self.n_gt_trees} GT tree generation jobs...")
        gt_job_ids = []
        
        for gt_id in range(1, self.n_gt_trees + 1):
            params = JobParameters(
                job_id=job_id,
                job_type='gt_generation',
                gt_id=gt_id,
                output_dir=str(self.base_dir / "gt_trees"),
            )
            job_params[job_id] = params
            gt_job_ids.append(job_id)
            job_id += 1
        
        # Phase 2: Cas9 Recording Jobs (4 tiers × N GT trees)
        self.logger.info(f"Generating parameters for {len(CAS9_TIERS) * self.n_gt_trees} Cas9 recording jobs...")
        cas9_job_ids = []
        
        for gt_id in range(1, self.n_gt_trees + 1):
            for tier_num in CAS9_TIERS.keys():
                params = JobParameters(
                    job_id=job_id,
                    job_type='cas9_recording',
                    gt_id=gt_id,
                    cas9_tier=tier_num,
                    output_dir=str(self.base_dir / "cas9_instances"),
                    depends_on=gt_job_ids  # Wait for all GT trees to be generated
                )
                job_params[job_id] = params
                cas9_job_ids.append(job_id)
                job_id += 1
        
        # Phase 3: Reconstruction Jobs (7 solvers × 4 tiers × N GT trees)
        self.logger.info(f"Generating parameters for {len(self.solvers) * len(CAS9_TIERS) * self.n_gt_trees} reconstruction jobs...")
        
        for gt_id in range(1, self.n_gt_trees + 1):
            for tier_num in CAS9_TIERS.keys():
                for solver in self.solvers:
                    params = JobParameters(
                        job_id=job_id,
                        job_type='reconstruction',
                        gt_id=gt_id,
                        cas9_tier=tier_num,
                        solver=solver,
                        output_dir=str(self.base_dir / "results"),
                        depends_on=cas9_job_ids  # Wait for all Cas9 recordings
                    )
                    job_params[job_id] = params
                    job_id += 1
        
        self.job_parameters = job_params
        
        # Save parameters to disk
        params_file = self.base_dir / "parameters" / "job_parameters.json"
        with open(params_file, 'w') as f:
            json.dump({k: asdict(v) for k, v in job_params.items()}, f, indent=2)
        
        total_jobs = len(job_params)
        self.logger.info(f"Generated {total_jobs} total job parameters")
        self.logger.info(f"- GT generation: {len(gt_job_ids)} jobs")
        self.logger.info(f"- Cas9 recording: {len(cas9_job_ids)} jobs")  
        self.logger.info(f"- Reconstruction: {total_jobs - len(gt_job_ids) - len(cas9_job_ids)} jobs")
        
        return job_params
    
    def save_cas9_tier_config(self):
        """Save Cas9 tier configurations to shared location."""
        config_file = self.base_dir / "parameters" / "cas9_tiers.json"
        
        # Convert tier configurations to JSON-serializable format
        tiers_json = {}
        for tier_num, tier in CAS9_TIERS.items():
            tiers_json[tier_num] = {
                'name': tier.name,
                'description': tier.description,
                'k': tier.k,
                'cassette_size': tier.cassette_size,
                'm': tier.m,
                'mutation_rates': tier.mutation_rates,
                'state_priors_exponent': tier.state_priors_exponent
            }
        
        with open(config_file, 'w') as f:
            json.dump(tiers_json, f, indent=2)
            
        self.logger.info(f"Saved Cas9 tier configurations to {config_file}")
    
    def create_lsf_script(self, job_params: JobParameters) -> str:
        """Create LSF submission script for a job."""

        script_name = f"{job_params.job_type}_{job_params.job_id}.lsf"
        script_path = self.base_dir / "scripts" / script_name

        # Get LSF resources from config with fallback to defaults
        lsf_resources = self.config.get('lsf', {}).get('resources', {})

        # Job-specific resource requirements
        if job_params.job_type == 'gt_generation':
            resources = lsf_resources.get('master_gt', {})
            walltime = resources.get('time_limit', '0:30')
            memory_gb = resources.get('memory_gb', 4.0)
            cores = resources.get('cores', 1)
            python_script = "generate_gt_tree.py"
            script_args = f"--gt_id {job_params.gt_id} --output_dir {job_params.output_dir}"

        elif job_params.job_type == 'cas9_recording':
            resources = lsf_resources.get('cas9_recording', {})
            walltime = resources.get('time_limit', '0:20')
            memory_gb = resources.get('memory_gb', 6.0)
            cores = resources.get('cores', 1)
            python_script = "apply_cas9_recording.py"
            script_args = f"--gt_id {job_params.gt_id} --cas9_tier {job_params.cas9_tier} --gt_dir {self.base_dir}/gt_trees --output_dir {job_params.output_dir}"

        elif job_params.job_type == 'reconstruction':
            resources = lsf_resources.get('reconstruction', {})
            walltime = resources.get('time_limit', '1:00')
            memory_gb = resources.get('memory_gb', 8.0)
            cores = resources.get('cores', 1)
            python_script = "reconstruct_and_analyze.py"
            script_args = f"--gt_id {job_params.gt_id} --cas9_tier {job_params.cas9_tier} --solver {job_params.solver} --cas9_dir {self.base_dir}/cas9_instances --output_dir {job_params.output_dir}"

        # Warn if memory per core is higher than 1.5GB (recommended value)
        if memory_gb > 1.5:
            self.logger.warning(
                f"WARNING: {job_params.job_type} job configured with {memory_gb}GB per core "
                f"(higher than recommended 1.5GB). Total memory will be {cores} × {memory_gb}GB = {cores * memory_gb}GB"
            )
        
        # Create dependency string
        depend_str = ""
        if job_params.depends_on:
            # For LSF, we'll use job names instead of IDs for dependencies
            depend_str = f"#BSUB -w \"done({job_params.job_type}_{job_params.job_id})\""
        
        # Create LSF script content
        script_content = f"""#!/bin/bash
#BSUB -J {job_params.job_type}_{job_params.job_id}
#BSUB -q {self.queue}
#BSUB -W {walltime}
#BSUB -R rusage[mem={memory_gb}GB]
#BSUB -n {cores} -R "span[hosts=1]"
#BSUB -o {self.base_dir}/logs/{job_params.job_type}_{job_params.job_id}.out
#BSUB -e {self.base_dir}/logs/{job_params.job_type}_{job_params.job_id}.err
{depend_str}

# Load necessary modules (adjust for your cluster)
# module load python/3.8
# module load stellars

# Change to working directory
cd {os.getcwd()}

# Run the Python script
python {python_script} {script_args}

echo "Job {job_params.job_id} ({job_params.job_type}) completed at $(date)"
"""
        
        # Write script to file
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make executable
        os.chmod(script_path, 0o755)
        
        return str(script_path)
    
    def submit_job_batch(self, job_type: str, dry_run: bool = False) -> Dict[int, str]:
        """Submit a batch of jobs of the same type."""
        
        jobs_to_submit = [params for params in self.job_parameters.values() 
                         if params.job_type == job_type]
        
        self.logger.info(f"Submitting {len(jobs_to_submit)} {job_type} jobs...")
        
        submitted = {}
        for job_params in jobs_to_submit:
            script_path = self.create_lsf_script(job_params)
            
            if dry_run:
                self.logger.info(f"DRY RUN: Would submit {script_path}")
                submitted[job_params.job_id] = f"DRY_RUN_{job_params.job_id}"
            else:
                try:
                    # Submit to LSF
                    result = subprocess.run(
                        ['bsub', '<', script_path], 
                        shell=True, 
                        capture_output=True, 
                        text=True
                    )
                    
                    if result.returncode == 0:
                        # Parse LSF job ID from output
                        lsf_job_id = result.stdout.strip().split()[1].strip('<>')
                        submitted[job_params.job_id] = lsf_job_id
                        self.logger.info(f"Submitted job {job_params.job_id} as LSF job {lsf_job_id}")
                    else:
                        self.logger.error(f"Failed to submit job {job_params.job_id}: {result.stderr}")
                        
                except Exception as e:
                    self.logger.error(f"Error submitting job {job_params.job_id}: {e}")
        
        return submitted
    
    def submit_all_jobs(self, dry_run: bool = False):
        """Submit all jobs in proper dependency order."""
        
        if not self.job_parameters:
            self.logger.error("No job parameters generated. Call generate_job_parameters() first.")
            return
        
        self.logger.info("Starting distributed job submission...")
        
        # Submit in phases to respect dependencies
        phases = ['gt_generation', 'cas9_recording', 'reconstruction']
        
        for phase in phases:
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"PHASE: {phase.upper()}")
            self.logger.info(f"{'='*60}")
            
            submitted = self.submit_job_batch(phase, dry_run=dry_run)
            self.submitted_jobs.update(submitted)
            
            # Wait a bit between phases to avoid overwhelming the scheduler
            if not dry_run and phase != phases[-1]:
                self.logger.info("Waiting 30 seconds before next phase...")
                time.sleep(30)
        
        self.logger.info(f"\nSubmitted {len(self.submitted_jobs)} total jobs")
        
        # Save submitted job tracking
        tracking_file = self.base_dir / "parameters" / "submitted_jobs.json"
        with open(tracking_file, 'w') as f:
            json.dump(self.submitted_jobs, f, indent=2)
    
    def monitor_progress(self) -> Dict[str, int]:
        """Monitor progress of submitted jobs."""
        
        if not self.submitted_jobs:
            self.logger.warning("No submitted jobs to monitor")
            return {}
        
        # Query LSF for job status
        try:
            result = subprocess.run(['bjobs', '-a'], capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.error(f"Failed to query job status: {result.stderr}")
                return {}
            
            # Parse job status
            status_counts = {'RUN': 0, 'PEND': 0, 'DONE': 0, 'EXIT': 0}
            
            for line in result.stdout.split('\n')[1:]:  # Skip header
                if line.strip():
                    parts = line.split()
                    if len(parts) >= 3:
                        status = parts[2]
                        if status in status_counts:
                            status_counts[status] += 1
            
            self.logger.info(f"Job Status: {status_counts}")
            return status_counts
            
        except Exception as e:
            self.logger.error(f"Error monitoring jobs: {e}")
            return {}
    
    def wait_for_completion(self, check_interval: int = 300):
        """Wait for all jobs to complete, checking periodically."""
        
        self.logger.info(f"Monitoring job completion (checking every {check_interval} seconds)...")
        
        while True:
            status = self.monitor_progress()
            
            if not status:
                self.logger.warning("Could not get job status, waiting...")
                time.sleep(check_interval)
                continue
            
            running_jobs = status.get('RUN', 0) + status.get('PEND', 0)
            
            if running_jobs == 0:
                self.logger.info("All jobs completed!")
                break
            
            self.logger.info(f"{running_jobs} jobs still running/pending...")
            time.sleep(check_interval)
    
    def aggregate_results(self) -> pd.DataFrame:
        """Aggregate results from all completed jobs."""
        
        results_dir = self.base_dir / "results"
        result_files = list(results_dir.glob("metrics_*.json"))
        
        self.logger.info(f"Found {len(result_files)} result files to aggregate")
        
        all_results = []
        
        for result_file in result_files:
            try:
                with open(result_file, 'r') as f:
                    result_data = json.load(f)
                all_results.append(result_data)
            except Exception as e:
                self.logger.warning(f"Could not load {result_file}: {e}")
        
        if not all_results:
            self.logger.error("No valid results found!")
            return pd.DataFrame()
        
        # Create DataFrame
        results_df = pd.DataFrame(all_results)
        
        # Save aggregated results
        output_file = self.base_dir / "aggregated" / "all_results.csv"
        results_df.to_csv(output_file, index=False)
        
        self.logger.info(f"Aggregated {len(results_df)} results to {output_file}")
        
        return results_df


def main():
    """Main function for LSF coordinator."""

    parser = argparse.ArgumentParser(description='LSF Coordinator for Multi-Tier Cas9 Analysis')
    parser.add_argument('--n_gt_trees', type=int, default=50, help='Number of GT trees to generate')
    parser.add_argument('--base_dir', type=str, default='/shared/cas9_analysis', help='Base directory for analysis')
    parser.add_argument('--queue', type=str, default='normal', help='LSF queue name')
    parser.add_argument('--config', type=str, help='Path to YAML configuration file')
    parser.add_argument('--dry_run', action='store_true', help='Generate scripts but do not submit jobs')
    parser.add_argument('--monitor_only', action='store_true', help='Only monitor existing jobs')
    parser.add_argument('--aggregate_only', action='store_true', help='Only aggregate existing results')

    args = parser.parse_args()

    # Create coordinator
    coordinator = LSFCoordinator(
        n_gt_trees=args.n_gt_trees,
        base_dir=args.base_dir,
        queue=args.queue,
        debug=args.dry_run,
        config_file=args.config
    )
    
    if args.monitor_only:
        # Only monitor existing jobs
        coordinator.wait_for_completion()
        
    elif args.aggregate_only:
        # Only aggregate existing results
        results_df = coordinator.aggregate_results()
        print(f"Aggregated {len(results_df)} results")
        
    else:
        # Full workflow
        print("="*80)
        print("LSF DISTRIBUTED MULTI-TIER CAS9 ANALYSIS")
        print("="*80)
        print(f"GT Trees: {args.n_gt_trees}")
        print(f"Total Jobs: ~{args.n_gt_trees * (1 + 4 + 28)} ({args.n_gt_trees} GT + {args.n_gt_trees*4} Cas9 + {args.n_gt_trees*28} Reconstruction)")
        print(f"Base Directory: {args.base_dir}")
        print(f"LSF Queue: {args.queue}")
        print(f"Dry Run: {args.dry_run}")
        
        # Generate parameters
        coordinator.generate_job_parameters()
        coordinator.save_cas9_tier_config()
        
        # Submit jobs
        coordinator.submit_all_jobs(dry_run=args.dry_run)
        
        if not args.dry_run:
            # Monitor completion
            coordinator.wait_for_completion()
            
            # Aggregate results
            results_df = coordinator.aggregate_results()
            print(f"\nFinal Results: {len(results_df)} reconstructions completed")


if __name__ == "__main__":
    main()