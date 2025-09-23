#!/usr/bin/env python3
"""
Dynamic Throttling Configuration Utility

Allows updating job throttling limits on-the-fly for running cascade workflows.
"""

import argparse
import sys
import yaml
from pathlib import Path
from typing import Dict, Any, Optional

def load_throttling_config(config_file: Path) -> Dict[str, Any]:
    """Load current throttling configuration"""
    if not config_file.exists():
        return {
            'throttling': {
                'max_concurrent_cas9_jobs': 20,
                'max_concurrent_reconstruction_jobs': 100,
                'job_submission_delay': 2.0,
                'batch_size': 10,
                'check_interval': 30
            },
            'metadata': {
                'description': 'Dynamic throttling configuration - edit this file to adjust job limits in real-time'
            }
        }

    try:
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"Error loading config file: {e}")
        sys.exit(1)

def save_throttling_config(config_file: Path, config: Dict[str, Any]):
    """Save throttling configuration"""
    try:
        with open(config_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        print(f"Successfully updated throttling config: {config_file}")
    except Exception as e:
        print(f"Error saving config file: {e}")
        sys.exit(1)

def show_current_config(config: Dict[str, Any]):
    """Display current throttling configuration"""
    throttling = config.get('throttling', {})

    print("Current Throttling Configuration:")
    print("=" * 40)
    print(f"Max CAS9 Jobs:           {throttling.get('max_concurrent_cas9_jobs', 'N/A')}")
    print(f"Max Reconstruction Jobs: {throttling.get('max_concurrent_reconstruction_jobs', 'N/A')}")
    print(f"Submission Delay:        {throttling.get('job_submission_delay', 'N/A')}s")
    print(f"Batch Size:              {throttling.get('batch_size', 'N/A')}")
    print(f"Check Interval:          {throttling.get('check_interval', 'N/A')}s")
    print()

def update_config(config: Dict[str, Any], args: argparse.Namespace) -> bool:
    """Update configuration with provided arguments"""
    throttling = config.setdefault('throttling', {})
    metadata = config.setdefault('metadata', {})

    updated = False

    if args.cas9_jobs is not None:
        old_val = throttling.get('max_concurrent_cas9_jobs', 'N/A')
        throttling['max_concurrent_cas9_jobs'] = args.cas9_jobs
        print(f"CAS9 Jobs:           {old_val} → {args.cas9_jobs}")
        updated = True

    if args.reconstruction_jobs is not None:
        old_val = throttling.get('max_concurrent_reconstruction_jobs', 'N/A')
        throttling['max_concurrent_reconstruction_jobs'] = args.reconstruction_jobs
        print(f"Reconstruction Jobs: {old_val} → {args.reconstruction_jobs}")
        updated = True

    if args.delay is not None:
        old_val = throttling.get('job_submission_delay', 'N/A')
        throttling['job_submission_delay'] = args.delay
        print(f"Submission Delay:    {old_val}s → {args.delay}s")
        updated = True

    if args.batch_size is not None:
        old_val = throttling.get('batch_size', 'N/A')
        throttling['batch_size'] = args.batch_size
        print(f"Batch Size:          {old_val} → {args.batch_size}")
        updated = True

    if args.check_interval is not None:
        old_val = throttling.get('check_interval', 'N/A')
        throttling['check_interval'] = args.check_interval
        print(f"Check Interval:      {old_val}s → {args.check_interval}s")
        updated = True

    if updated:
        # Update metadata
        import time
        metadata['last_updated'] = time.strftime('%Y-%m-%d %H:%M:%S')
        metadata['updated_by'] = f"update_throttling.py (PID: {os.getpid()})"

    return updated

def reset_to_original_config(config: Dict[str, Any], shared_dir: Path) -> bool:
    """Reset throttling configuration to original LSF settings from cascade config"""
    try:
        # Load the original cascade config
        cascade_config_file = shared_dir / "cascade_config.yaml"

        if not cascade_config_file.exists():
            print(f"Error: Cannot find original config file: {cascade_config_file}")
            return False

        with open(cascade_config_file, 'r') as f:
            cascade_config = yaml.safe_load(f)

        # Extract original LSF throttling settings
        original_throttling = cascade_config.get('lsf', {}).get('throttling', {})

        if not original_throttling:
            print("Warning: No original throttling configuration found in cascade config")
            print("Using default values instead")
            original_throttling = {
                'max_concurrent_cas9_jobs': 50,
                'max_concurrent_reconstruction_jobs': 200,
                'job_submission_delay': 1.0,
                'batch_size': 20,
                'check_interval': 30
            }

        throttling = config.setdefault('throttling', {})

        print("Resetting to original LSF configuration:")
        for key, value in original_throttling.items():
            old_val = throttling.get(key, 'N/A')
            throttling[key] = value
            if key.endswith('_delay') or key.endswith('_interval'):
                print(f"  {key}: {old_val} → {value}s")
            else:
                print(f"  {key}: {old_val} → {value}")

        # Update metadata
        import time
        import os
        metadata = config.setdefault('metadata', {})
        metadata['last_updated'] = time.strftime('%Y-%m-%d %H:%M:%S')
        metadata['reset_to_original'] = True
        metadata['updated_by'] = f"update_throttling.py (PID: {os.getpid()})"

        return True

    except Exception as e:
        print(f"Error resetting to original configuration: {e}")
        return False

def apply_preset(config: Dict[str, Any], preset: str) -> bool:
    """Apply a predefined throttling preset"""
    throttling = config.setdefault('throttling', {})

    presets = {
        'conservative': {
            'max_concurrent_cas9_jobs': 5,
            'max_concurrent_reconstruction_jobs': 25,
            'job_submission_delay': 10.0,
            'batch_size': 3,
            'check_interval': 60
        },
        'moderate': {
            'max_concurrent_cas9_jobs': 15,
            'max_concurrent_reconstruction_jobs': 75,
            'job_submission_delay': 3.0,
            'batch_size': 8,
            'check_interval': 45
        },
        'aggressive': {
            'max_concurrent_cas9_jobs': 40,
            'max_concurrent_reconstruction_jobs': 200,
            'job_submission_delay': 1.0,
            'batch_size': 15,
            'check_interval': 30
        }
    }

    if preset not in presets:
        print(f"Unknown preset: {preset}")
        print(f"Available presets: {', '.join(presets.keys())}")
        return False

    preset_config = presets[preset]

    print(f"Applying '{preset}' preset:")
    for key, value in preset_config.items():
        old_val = throttling.get(key, 'N/A')
        throttling[key] = value
        if key.endswith('_delay') or key.endswith('_interval'):
            print(f"  {key}: {old_val} → {value}s")
        else:
            print(f"  {key}: {old_val} → {value}")

    # Update metadata
    import time
    import os
    metadata = config.setdefault('metadata', {})
    metadata['last_updated'] = time.strftime('%Y-%m-%d %H:%M:%S')
    metadata['preset_applied'] = preset
    metadata['updated_by'] = f"update_throttling.py (PID: {os.getpid()})"

    return True

def main():
    parser = argparse.ArgumentParser(
        description="Update dynamic throttling configuration for cascade workflows",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Show current configuration
  python update_throttling.py /path/to/shared_dir

  # Reduce CAS9 job limit to 5
  python update_throttling.py /path/to/shared_dir --cas9-jobs 5

  # Increase reconstruction jobs and reduce delay
  python update_throttling.py /path/to/shared_dir --reconstruction-jobs 150 --delay 1.0

  # Apply conservative preset
  python update_throttling.py /path/to/shared_dir --preset conservative

  # Reset to original LSF configuration
  python update_throttling.py /path/to/shared_dir --reset

  # Multiple updates
  python update_throttling.py /path/to/shared_dir \\
    --cas9-jobs 10 --reconstruction-jobs 50 --batch-size 5
        """
    )

    parser.add_argument('shared_dir',
                       help='Path to the shared directory containing throttling_config.yaml')

    parser.add_argument('--cas9-jobs', type=int,
                       help='Maximum concurrent CAS9 recording jobs')

    parser.add_argument('--reconstruction-jobs', type=int,
                       help='Maximum concurrent reconstruction jobs')

    parser.add_argument('--delay', type=float,
                       help='Delay between job submissions (seconds)')

    parser.add_argument('--batch-size', type=int,
                       help='Number of jobs to submit per batch')

    parser.add_argument('--check-interval', type=int,
                       help='Interval between job count checks (seconds)')

    parser.add_argument('--preset', choices=['conservative', 'moderate', 'aggressive'],
                       help='Apply a predefined configuration preset')

    parser.add_argument('--reset', action='store_true',
                       help='Reset to original LSF throttling configuration from cascade config')

    parser.add_argument('--show-only', action='store_true',
                       help='Only show current configuration, do not update')

    args = parser.parse_args()

    # Validate shared directory
    shared_dir = Path(args.shared_dir)
    if not shared_dir.exists():
        print(f"Error: Shared directory does not exist: {shared_dir}")
        sys.exit(1)

    config_file = shared_dir / "throttling_config.yaml"

    # Load current configuration
    config = load_throttling_config(config_file)

    # Show current configuration
    show_current_config(config)

    if args.show_only:
        return

    updated = False

    # Reset to original configuration if specified
    if args.reset:
        updated = reset_to_original_config(config, shared_dir)

    # Apply preset if specified
    if args.preset:
        if updated:
            print("\nAdditional preset application:")
        preset_updated = apply_preset(config, args.preset)
        updated = updated or preset_updated

    # Apply individual parameter updates
    if any([args.cas9_jobs is not None, args.reconstruction_jobs is not None,
            args.delay is not None, args.batch_size is not None,
            args.check_interval is not None]):
        if updated:
            print("\nAdditional updates:")
        individual_updated = update_config(config, args)
        updated = updated or individual_updated

    # Save configuration if anything was updated
    if updated:
        print()
        save_throttling_config(config_file, config)
        print("\nConfiguration will be automatically picked up within 30 seconds.")
        print("Watch the workflow logs for confirmation of the update.")
    else:
        if not args.show_only:
            print("No updates specified. Use --help for usage information.")

if __name__ == "__main__":
    import os
    main()