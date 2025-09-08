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
from config_loader import load_config

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


def copy_worker_scripts(shared_dir: Path, source_dir: Path) -> None:
    """Copy worker scripts and job templates to shared directory."""
    logger.info("Copying worker scripts to shared directory...")
    
    scripts = [
        "master_gt_worker.py",
        "cas9_recording_worker.py", 
        "reconstruction_worker.py",
        "reconstruction_metrics_table_multi_tier.py",
        "job_monitor.py"
    ]
    
    for script in scripts:
        source = source_dir / script
        dest = shared_dir / script
        
        if source.exists():
            subprocess.run(['cp', str(source), str(dest)], check=True)
            logger.info(f"Copied {script} to shared directory")
        else:
            logger.warning(f"Script not found: {source}")
    
    # Copy LSF job templates with path substitution
    copy_job_templates_with_substitution(shared_dir, source_dir)


def substitute_placeholders(template_content: str, shared_dir: Path) -> str:
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
    
    # Replace placeholders with actual paths
    substituted = template_content.replace('{SHARED_DIR}', str(shared_dir_abs))
    substituted = substituted.replace('{CONDA_PREFIX}', conda_prefix)
    
    return substituted


def copy_job_templates_with_substitution(shared_dir: Path, source_dir: Path) -> None:
    """Copy job templates and substitute placeholders with actual paths."""
    logger.info("Copying LSF job templates...")
    
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
            substituted_content = substitute_placeholders(template_content, shared_dir)
            
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
    
    # Read the original config and convert relative shared_dir to absolute
    with open(config_path, 'r') as f:
        config_content = f.read()
    
    # Update the shared_dir to use absolute path to prevent nested directory issues
    import re
    config_content = re.sub(
        r'shared_dir:\s*["\']?shared_medium["\']?',
        f'shared_dir: "{shared_dir.resolve()}"',
        config_content
    )
    
    # Write the updated config
    with open(config_dest, 'w') as f:
        f.write(config_content)
    
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
    logger.info(f"Expected total jobs: 1 + 4 + {4 * len(config_obj.config['solvers']['enabled']) * config_obj.config['execution']['num_gt_instances']}")
    
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
        copy_worker_scripts(shared_dir, source_dir)
        
        if args.setup_only:
            logger.info("Setup complete. Use --submit to launch the job cascade.")
            return
        
        if args.local_test:
            logger.info("Running master GT worker locally for testing...")
            # Run master worker directly
            cmd = [
                'python', str(source_dir / 'master_gt_worker.py'),
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
            logger.info("\nTo monitor progress:")
            logger.info(f"  python {source_dir / 'job_monitor.py'} --shared_dir {shared_dir}")
            logger.info("\nTo monitor continuously:")
            logger.info(f"  python {source_dir / 'job_monitor.py'} --shared_dir {shared_dir} --continuous")
            logger.info("\nExpected job flow:")
            logger.info("  1. Master GT job (generates tree, submits 4 Cas9 jobs)")
            logger.info("  2. 4 Cas9 recording jobs (each submits 7 reconstruction jobs)")
            logger.info("  3. 28 reconstruction jobs (final analysis)")
            logger.info("  Total: 33 jobs")
        
    except Exception as e:
        logger.error(f"Launch failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
