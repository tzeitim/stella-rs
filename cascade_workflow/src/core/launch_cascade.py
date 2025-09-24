#!/usr/bin/env python3
"""
Launch Cascade - Entry point for LSF Cascading Job System

Configurable script to launch multiple GT instances with flexible parameters.
Supports reproducible seeds, multiple instances, and parameterized execution.
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path
import importlib.util

# Import solver validation functions and config loader
sys.path.append(str(Path(__file__).parent.parent / "utils"))
from config_loader import load_config
from solver_config import validate_requested_solvers, get_enabled_solvers, get_supported_solvers_info

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def setup_shared_directory(shared_dir: Path) -> None:
    """Set up the shared directory structure."""
    logger.info(f"Setting up shared directory: {shared_dir}")
    
    directories = [
        "jobs",
        "logs", 
        "status",
        "gt_trees",
        "cas9_instances",
        "results"
    ]
    
    for directory in directories:
        dir_path = shared_dir / directory
        dir_path.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created directory: {dir_path}")


def copy_worker_scripts(shared_dir: Path, source_dir: Path, conda_environment: str = "cas11") -> None:
    """Copy worker scripts and job templates to shared directory."""
    logger.info("Copying worker scripts to shared directory...")
    
    # Determine if we're in the organized structure or original
    workers_dir = source_dir.parent / "workers"
    utils_dir = source_dir.parent / "utils"
    
    # If not in organized structure, use source_dir directly
    if not workers_dir.exists():
        workers_dir = source_dir
        utils_dir = source_dir
    
    # Worker scripts mapping
    script_locations = {
        "master_gt_worker.py": workers_dir,
        "cas9_recording_worker.py": workers_dir,
        "reconstruction_worker.py": workers_dir,
        "reconstruction_metrics_table_multi_tier.py": workers_dir,
        "job_monitor.py": workers_dir,
        "solver_config.py": utils_dir,
        "partitioned_results_writer.py": utils_dir,
        "config_loader.py": utils_dir,
        "job_throttling.py": utils_dir,
        "update_throttling.py": utils_dir
    }
    
    for script, location in script_locations.items():
        source = location / script
        dest = shared_dir / script
        
        # Fallback to source_dir if not found
        if not source.exists():
            source = source_dir / script
        
        if source.exists():
            subprocess.run(['cp', str(source), str(dest)], check=True)
            logger.info(f"Copied {script} to shared directory")
        else:
            logger.warning(f"Script not found: {source}")
    
    # Copy LSF job templates with path substitution
    copy_job_templates_with_substitution(shared_dir, source_dir, conda_environment)


def substitute_placeholders(template_content: str, shared_dir: Path, conda_environment: str = "cas11") -> str:
    """Replace placeholders in job templates with actual paths."""
    import os

    # Get absolute path for shared directory
    shared_dir_abs = shared_dir.resolve()

    # Get conda prefix dynamically
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if not conda_prefix:
        # Fallback: try to detect common conda installations
        home_dir = os.path.expanduser('~')
        for conda_path in ['miniforge3', 'miniconda3', 'anaconda3']:
            potential_path = os.path.join(home_dir, conda_path)
            if os.path.exists(potential_path):
                conda_prefix = potential_path
                break

        if not conda_prefix:
            raise RuntimeError("Could not detect conda installation. Set CONDA_PREFIX environment variable.")

    logger.info(f"Using conda prefix: {conda_prefix}")
    logger.info(f"Using conda environment: {conda_environment}")

    # Replace placeholders with actual paths
    substituted = template_content.replace('{SHARED_DIR}', str(shared_dir_abs))
    substituted = substituted.replace('{CONDA_PREFIX}', conda_prefix)
    substituted = substituted.replace('{CONDA_ENVIRONMENT}', conda_environment)

    return substituted


def copy_job_templates_with_substitution(shared_dir: Path, source_dir: Path, conda_environment: str = "cas11") -> None:
    """Copy job templates and substitute placeholders with actual paths."""
    logger.info("Copying LSF job templates...")
    
    # Check for organized structure first
    jobs_source_dir = source_dir.parent.parent / "lsf_templates"
    
    # Fallback to original structure
    if not jobs_source_dir.exists():
        jobs_source_dir = source_dir / "jobs"
    
    if not jobs_source_dir.exists():
        logger.warning(f"Jobs directory not found: {jobs_source_dir}")
        return
    
    jobs_dest_dir = shared_dir / "jobs" 
    jobs_dest_dir.mkdir(exist_ok=True)
    
    job_templates = [
        "master_gt_job.lsf",
        "cas9_recording_job_template.lsf",
        "reconstruction_job_template.lsf"
    ]
    
    for template in job_templates:
        source_path = jobs_source_dir / template
        dest_path = jobs_dest_dir / template
        
        if source_path.exists():
            # Read template content
            with open(source_path, 'r') as f:
                template_content = f.read()
            
            # Substitute placeholders
            substituted_content = substitute_placeholders(template_content, shared_dir, conda_environment)
            
            # Write to destination
            with open(dest_path, 'w') as f:
                f.write(substituted_content)
            
            # Make executable
            dest_path.chmod(0o755)
            logger.info(f"Copied {template} to shared directory with path substitution")
        else:
            logger.warning(f"Job template not found: {source_path}")


def submit_master_job(shared_dir: Path, config_path: str) -> str:
    """Submit the master GT job to LSF."""
    logger.info("Submitting master GT job to LSF...")
    
    job_script = shared_dir / "jobs" / "master_gt_job.lsf"
    
    # Make sure job script exists
    if not job_script.exists():
        raise FileNotFoundError(f"Master job script not found: {job_script}")
    
    # Copy and update config file to shared directory with absolute paths
    config_dest = shared_dir / "cascade_config.yaml"
    
    # Load the config as a structured object and update shared_dir to absolute path
    import yaml
    with open(config_path, 'r') as f:
        config_data = yaml.safe_load(f)
    
    # Update the shared_dir to use absolute path to prevent nested directory issues
    if 'output' not in config_data:
        config_data['output'] = {}
    config_data['output']['shared_dir'] = str(shared_dir.resolve())
    
    # Write the updated config
    with open(config_dest, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)
    
    logger.info(f"Copied and updated configuration to: {config_dest}")
    
    # Submit job with config path as environment variable
    env_vars = f"CONFIG_PATH={config_dest}"
    cmd = ['bsub', '-env', env_vars, str(job_script)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Failed to submit master job: {result.stderr}")
    
    # Extract job ID
    job_id = result.stdout.strip().split('<')[1].split('>')[0]
    logger.info(f"Master GT job submitted with ID: {job_id}")
    
    return job_id


def main():
    parser = argparse.ArgumentParser(description="Launch Cascade - Start the multi-tier Cas9 analysis")
    parser.add_argument('--shared_dir', 
                       help='Shared directory for distributed computing (optional - uses config if not specified)')
    parser.add_argument('--config', default='cascade_config.yaml',
                       help='Configuration file path')
    parser.add_argument('--setup_only', action='store_true',
                       help='Only set up directories, do not submit jobs')
    parser.add_argument('--local_test', action='store_true',
                       help='Run locally for testing (not on LSF)')
    
    args = parser.parse_args()
    
    # Load configuration
    config_obj = load_config(args.config)
    logger.info(f"Loaded configuration from: {args.config}")
    logger.info(f"GT instances to generate: {config_obj.config['execution']['num_gt_instances']}")
    logger.info(f"Enabled solvers: {config_obj.config['solvers']['enabled']}")

    # Validate that all requested solvers are supported
    requested_solvers = config_obj.config['solvers']['enabled']
    is_valid, unsupported_solvers = validate_requested_solvers(requested_solvers)

    if not is_valid:
        supported_solvers = get_enabled_solvers()
        logger.error("=" * 60)
        logger.error("❌ SOLVER VALIDATION FAILED!")
        logger.error(f"Unsupported solvers in config: {unsupported_solvers}")
        logger.error(f"Supported solvers in workflow: {supported_solvers}")
        logger.error("=" * 60)
        logger.error("Please update your configuration to only include supported solvers.")
        logger.error("Alternatively, enable the required solvers in solver_config.py")
        sys.exit(1)

    logger.info("✓ Solver validation passed")

    # Calculate expected total jobs with cas9_simulations_per_gt
    num_gt_instances = config_obj.config['execution']['num_gt_instances']
    cas9_simulations_per_gt = config_obj.config['execution'].get('cas9_simulations_per_gt', 1)
    num_tiers = len(config_obj.config.get('cas9_tiers', {1: {}, 2: {}, 3: {}, 4: {}}))
    num_solvers = len(config_obj.config['solvers']['enabled'])
    reconstructions_per_solver = config_obj.config['solvers'].get('reconstructions_per_solver', 1)
    
    cas9_jobs = num_gt_instances * cas9_simulations_per_gt * num_tiers
    reconstruction_jobs = cas9_jobs * num_solvers * reconstructions_per_solver
    total_jobs = 1 + cas9_jobs + reconstruction_jobs
    
    logger.info(f"Cas9 simulations per GT: {cas9_simulations_per_gt}")
    logger.info(f"Expected total jobs: 1 + {cas9_jobs} + {reconstruction_jobs} = {total_jobs}")

    # Extract conda environment from config
    conda_environment = config_obj.config['execution'].get('conda_environment', 'cas11')
    logger.info(f"Using conda environment: {conda_environment}")

    # Use config's shared_dir if command line argument not provided
    config_shared_dir = config_obj.config.get('output', {}).get('shared_dir')
    if args.shared_dir:
        shared_dir = Path(args.shared_dir)
        logger.info(f"Using shared directory from command line: {shared_dir}")
    elif config_shared_dir:
        shared_dir = Path(config_shared_dir)
        logger.info(f"Using shared directory from config: {shared_dir}")
    else:
        raise ValueError("No shared directory specified. Provide via --shared_dir or in config file under output.shared_dir")
    source_dir = Path(__file__).parent
    
    try:
        # Step 1: Set up shared directory
        setup_shared_directory(shared_dir)
        
        # Step 2: Copy worker scripts
        copy_worker_scripts(shared_dir, source_dir, conda_environment)
        
        if args.setup_only:
            logger.info("Setup complete. Use --submit to launch the job cascade.")
            return
        
        if args.local_test:
            logger.info("Running master GT worker locally for testing...")
            # Find master worker script
            master_worker = source_dir.parent / "workers" / "master_gt_worker.py"
            if not master_worker.exists():
                master_worker = source_dir / "master_gt_worker.py"
            
            # Run master worker directly
            cmd = [
                'python', str(master_worker),
                '--output_dir', str(shared_dir / 'gt_trees'),
                '--shared_dir', str(shared_dir)
            ]
            subprocess.run(cmd, check=True)
            logger.info("Local test completed successfully!")
        else:
            # Step 3: Submit master job
            job_id = submit_master_job(shared_dir, args.config)
            
            logger.info("=" * 60)
            logger.info("🚀 CASCADE LAUNCHED SUCCESSFULLY!")
            logger.info(f"📋 Master Job ID: {job_id}")
            logger.info(f"📁 Shared Directory: {shared_dir}")
            logger.info("=" * 60)
            # Determine correct monitor script path - check both locations
            workers_dir = Path(__file__).parent.parent / "workers"
            if workers_dir.exists():
                # We're in the organized structure
                monitor_script = workers_dir / "job_monitor.py"
            else:
                # We're running from original location
                monitor_script = Path(__file__).parent / "job_monitor.py"
            
            logger.info("\nTo monitor progress:")
            logger.info(f"  python {monitor_script} --shared_dir {shared_dir}")
            logger.info("\nTo monitor continuously:")
            logger.info(f"  python {monitor_script} --shared_dir {shared_dir} --continuous")

            # Determine aggregation script path
            utils_dir = Path(__file__).parent.parent / "utils"
            if utils_dir.exists():
                aggregate_script = utils_dir / "aggregate_results.py"
            else:
                aggregate_script = Path(__file__).parent / "aggregate_results.py"

            logger.info("\nTo aggregate results when complete:")
            logger.info(f"  python {aggregate_script} --results_dir {shared_dir}/partitioned_results --output {shared_dir}/consolidated_results.parquet")

            # Also show the dedicated partitioned results reader
            partitioned_reader = utils_dir / "read_partitioned_results.py"
            if partitioned_reader.exists():
                logger.info(f"\nTo read/analyze partitioned results:")
                logger.info(f"  python {partitioned_reader} --results_dir {shared_dir}/partitioned_results --summary")

            # Determine diagram script path - check root directory
            root_dir = Path(__file__).parent.parent.parent.parent
            diagram_script = root_dir / "generate_workflow_diagram.py"

            logger.info("\nTo generate workflow diagram:")
            logger.info(f"  python {diagram_script} {args.config}")


            logger.info("\nExpected job flow:")
            logger.info(f"  1. Master GT job (generates {num_gt_instances} GT trees, submits {cas9_jobs} Cas9 jobs)")
            logger.info(f"  2. {cas9_jobs} Cas9 recording jobs (each submits {num_solvers * reconstructions_per_solver} reconstruction jobs)")
            logger.info(f"  3. {reconstruction_jobs} reconstruction jobs (final analysis)")
            logger.info(f"  Total: {total_jobs} jobs")
        
    except Exception as e:
        logger.error(f"Launch failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
