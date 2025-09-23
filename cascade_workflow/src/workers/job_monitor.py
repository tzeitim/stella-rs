#!/usr/bin/env python3
"""
Job Monitor - Status Tracking for LSF Cascading Job System

Monitors the status of all job levels and provides comprehensive status reporting
for the multi-tier Cas9 analysis distributed computing workflow.
"""

import argparse
import json
import logging
import signal
import subprocess
import sys
import time
import yaml
import polars as pl 
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.warning(f"Could not load config {config_path}: {e}")
        return {}


class JobMonitor:
    """Monitor for tracking LSF cascading job system status."""
    
    def __init__(self, shared_dir: str):
        self.shared_dir = Path(shared_dir)
        self.status_dir = self.shared_dir / "status"
        self.logs_dir = self.shared_dir / "logs"
        self.results_dir = self.shared_dir / "results"
        
        # Check if shared directory exists
        if not self.shared_dir.exists():
            raise FileNotFoundError(f"Shared directory does not exist: {self.shared_dir}")
        
        # Check if config exists in shared directory
        primary_config = self.shared_dir / "cascade_config.yaml"
        if not primary_config.exists():
            raise FileNotFoundError(f"Configuration file not found in shared directory: {primary_config}")
        
        # Only create status directory after validation
        self.status_dir.mkdir(parents=True, exist_ok=True)
        
        # Load configuration from shared directory only
        config_paths = [primary_config]
        
        self.config = {}
        for config_path in config_paths:
            self.config = load_config(str(config_path))
            if self.config:
                logger.info(f"Loaded configuration from: {config_path}")
                logger.info(f"Config keys: {list(self.config.keys())}")
                logger.info(f"Solvers section: {self.config.get('solvers')}")
                break
        else:
            logger.warning("No configuration file found, using defaults")
        
        # Build expected job structure from configuration
        self.expected_structure = self._build_expected_structure()
    
    def _build_expected_structure(self) -> Dict[str, List[str]]:
        """Build expected job structure from configuration."""
        structure = {
            'level1': [],
            'level2': [],
            'level3': []
        }
        
        # Get number of GT instances and Cas9 simulations per GT
        num_gt_instances = self.config.get('execution', {}).get('num_gt_instances', 1)
        cas9_simulations_per_gt = self.config.get('execution', {}).get('cas9_simulations_per_gt', 1)
        
        # Build level 1 (GT generation) components - one per instance
        for instance in range(num_gt_instances):
            structure['level1'].append(f'gt_instance_{instance}')
        
        # Get active tiers and solvers from config
        active_tiers = []
        if self.config.get('cas9_tiers'):
            active_tiers = list(self.config['cas9_tiers'].keys())
        else:
            # Default to all tiers if not specified
            active_tiers = [1, 2, 3, 4]
        
        active_solvers = []
        if self.config.get('solvers', {}).get('enabled'):
            active_solvers = self.config['solvers']['enabled']
        else:
            # Default solvers if not specified
            active_solvers = ['nj', 'maxcut', 'greedy', 'spectral', 'smj', 'dmj', 'ilp']
        
        logger.info(f"GT instances: {num_gt_instances}")
        logger.info(f"Cas9 simulations per GT: {cas9_simulations_per_gt}")
        logger.info(f"Active tiers: {active_tiers}")  
        logger.info(f"Active solvers: {active_solvers}")
        
        # Build level 2 (CAS9 recording) components - multiply by instances and simulations
        for instance in range(num_gt_instances):
            for sim in range(cas9_simulations_per_gt):
                for tier in active_tiers:
                    structure['level2'].append(f'instance{instance}_sim{sim}_tier{tier}_recording')
        
        # Build level 3 (reconstruction) components - multiply by instances, simulations, tiers, solvers, reconstructions
        reconstructions_per_solver = self.config.get('solvers', {}).get('reconstructions_per_solver', 1)
        for instance in range(num_gt_instances):
            for sim in range(cas9_simulations_per_gt):
                for tier in active_tiers:
                    for solver in active_solvers:
                        for recon_id in range(reconstructions_per_solver):
                            structure['level3'].append(f'instance{instance}_sim{sim}_recon{recon_id}_tier{tier}_{solver}')
        
        logger.info(f"Expected structure: Level1={len(structure['level1'])}, "
                   f"Level2={len(structure['level2'])}, Level3={len(structure['level3'])}")
        
        return structure
    
    def load_level_status(self, level: str) -> Dict[str, Any]:
        """Load status for a specific level."""
        status_file = self.status_dir / f"{level}_status.json"
        
        if not status_file.exists():
            return {}
        
        try:
            with open(status_file, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse {status_file}: {e}")
            return {}
    
    def get_lsf_job_status(self, job_name_pattern: str) -> List[Dict[str, str]]:
        """Query LSF for job status using bjobs command."""
        try:
            cmd = ['bjobs', '-J', job_name_pattern, '-o', 'jobid stat job_name submit_time']
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                # No jobs found is normal
                return []
            
            # Parse bjobs output
            lines = result.stdout.strip().split('\n')
            if len(lines) <= 1:  # Header only
                return []
            
            jobs = []
            for line in lines[1:]:  # Skip header
                parts = line.split()
                if len(parts) >= 4:
                    jobs.append({
                        'job_id': parts[0],
                        'status': parts[1],
                        'job_name': parts[2],
                        'submit_time': ' '.join(parts[3:])
                    })
            
            return jobs
            
        except Exception as e:
            logger.error(f"Failed to query LSF jobs: {e}")
            return []
    
    def check_file_outputs(self) -> Dict[str, Any]:
        """Check for expected output files."""
        num_gt_instances = self.config.get('execution', {}).get('num_gt_instances', 1)
        
        file_status = {
            'gt_trees': {},
            'cas9_instances': {},
            'results': {}
        }
        
        # Check for GT trees - one per instance
        for instance in range(num_gt_instances):
            gt_tree_paths = [
                # Correct naming pattern used by master_gt_worker.py
                self.shared_dir / "gt_trees" / f"gt_tree_instance_{instance}.pkl",
                self.shared_dir / f"gt_tree_instance_{instance}.pkl",
                # Legacy patterns for compatibility
                self.shared_dir / f"gt_instance_{instance:03d}.pkl",
                self.shared_dir / "gt_trees" / f"gt_instance_{instance:03d}.pkl",
                self.shared_dir.parent / f"gt_instance_{instance:03d}.pkl"
            ]
            file_status['gt_trees'][f'instance_{instance}'] = any(path.exists() for path in gt_tree_paths)
        
        # Get cas9_simulations_per_gt from config
        cas9_simulations_per_gt = self.config.get('execution', {}).get('cas9_simulations_per_gt', 1)
        
        # Check Cas9 instances - multiple per GT instance and simulation
        # Get active tiers from config instead of hardcoded 1-4
        active_tiers = list(self.config.get('cas9_tiers', {1: {}, 2: {}, 3: {}, 4: {}}).keys())

        for instance in range(num_gt_instances):
            for sim in range(cas9_simulations_per_gt):
                for tier in active_tiers:
                    # Try multiple naming patterns
                    cas9_file_paths = [
                        # New naming with instance and simulation
                        self.shared_dir / "cas9_instances" / f"instance{instance}_sim{sim}_tier{tier}_instance.pkl",
                        self.shared_dir.parent / "cas9_instances" / f"instance{instance}_sim{sim}_tier{tier}_instance.pkl",
                        # Legacy naming without simulation ID (for backward compatibility)
                        self.shared_dir / "cas9_instances" / f"instance{instance}_tier{tier}_instance.pkl",
                        self.shared_dir.parent / "cas9_instances" / f"instance{instance}_tier{tier}_instance.pkl"
                    ]
                    file_status['cas9_instances'][f'instance{instance}_sim{sim}_tier{tier}'] = any(path.exists() for path in cas9_file_paths)
        
        # Check results - Updated for partitioned parquet format
        active_solvers = self.config.get('solvers', {}).get('enabled', ['nj', 'maxcut', 'greedy', 'vanilla'])
        reconstructions_per_solver = self.config.get('solvers', {}).get('reconstructions_per_solver', 1)

        for instance in range(num_gt_instances):
            for sim in range(cas9_simulations_per_gt):
                for tier in active_tiers:
                    for solver in active_solvers:
                        for recon_id in range(reconstructions_per_solver):
                            result_id = f'instance{instance}_sim{sim}_recon{recon_id}_tier{tier}_{solver}'

                            # Check in partitioned parquet structure first (using equals format)
                            partitioned_dir = self.shared_dir / "partitioned_results" / f"cas9_tier={tier}" / f"solver={solver}"
                            result_found = False

                            if partitioned_dir.exists():
                                try:
                                    parquet_files = list(partitioned_dir.glob("*.parquet"))
                                    if parquet_files:
                                        # Read parquet and check if this specific result exists
                                        try:
                                            combined_df = pl.concat([pl.read_parquet(f) for f in parquet_files])
                                            # Check if our specific reconstruction exists in the data
                                            matching_rows = combined_df.filter(
                                                (pl.col("gt_instance_id") == instance) &
                                                (pl.col("cas9_simulation_id") == sim) &
                                                (pl.col("reconstruction_num") == recon_id) &
                                                (pl.col("cas9_tier") == tier) &
                                                (pl.col("solver") == solver)
                                            )
                                            result_found = len(matching_rows) > 0
                                        except Exception as e:
                                            # If there's an error reading parquet, fallback to old method
                                            logger.debug(f"Error reading partitioned results for {result_id}: {e}")
                                            result_found = False
                                except (FileNotFoundError, PermissionError, OSError) as e:
                                    # If directory scanning fails, fallback to old method
                                    logger.debug(f"Error scanning partitioned directory {partitioned_dir}: {e}")
                                    result_found = False

                            # Fallback to old individual file checking if not found in partitioned results
                            if not result_found:
                                result_file_paths = [
                                    # New naming with full indexing (instance, sim, recon)
                                    self.shared_dir / "results" / f"instance{instance}_sim{sim}_recon{recon_id}_tier{tier}_{solver}_metrics.json",
                                    self.shared_dir.parent / "results" / f"instance{instance}_sim{sim}_recon{recon_id}_tier{tier}_{solver}_metrics.json",
                                    # Previous naming without recon ID
                                    self.shared_dir / "results" / f"instance{instance}_sim{sim}_tier{tier}_{solver}_metrics.json",
                                    self.shared_dir.parent / "results" / f"instance{instance}_sim{sim}_tier{tier}_{solver}_metrics.json",
                                    # Legacy naming without instance/simulation
                                    self.shared_dir / "results" / f"tier{tier}_{solver}_metrics.json",
                                    self.shared_dir.parent / "results" / f"tier{tier}_{solver}_metrics.json",
                                    # Parquet versions
                                    self.shared_dir / "results" / f"instance{instance}_sim{sim}_recon{recon_id}_tier{tier}_{solver}_metrics.parquet",
                                    self.shared_dir.parent / "results" / f"instance{instance}_sim{sim}_recon{recon_id}_tier{tier}_{solver}_metrics.parquet",
                                    self.shared_dir / "results" / f"instance{instance}_sim{sim}_tier{tier}_{solver}_metrics.parquet",
                                    self.shared_dir.parent / "results" / f"instance{instance}_sim{sim}_tier{tier}_{solver}_metrics.parquet",
                                    self.shared_dir / "results" / f"tier{tier}_{solver}_metrics.parquet",
                                    self.shared_dir.parent / "results" / f"tier{tier}_{solver}_metrics.parquet"
                                ]
                                result_found = any(path.exists() for path in result_file_paths)

                            file_status['results'][result_id] = result_found
        
        return file_status
    
    def generate_status_summary(self) -> Dict[str, Any]:
        """Generate comprehensive status summary."""
        summary = {
            'timestamp': datetime.now().isoformat(),
            'levels': {},
            'lsf_jobs': {},
            'file_outputs': self.check_file_outputs(),
            'overall_progress': {}
        }
        
        # Check each level
        for level in ['level1', 'level2', 'level3']:
            level_status = self.load_level_status(level)

            # If no status file exists, infer completion from file outputs
            completed_components = []
            failed_components = []

            if level_status:
                # Use status file data if available
                completed_components = [
                    comp for comp in level_status
                    if level_status[comp].get('status') == 'completed'
                ]
                failed_components = [
                    comp for comp in level_status
                    if level_status[comp].get('status') == 'failed'
                ]
            else:
                # Infer completion from file existence
                file_outputs = summary['file_outputs']

                if level == 'level1':
                    # GT generation - check if GT trees exist
                    if file_outputs.get('gt_trees', {}):
                        completed_count = sum(1 for exists in file_outputs['gt_trees'].values() if exists)
                        if completed_count == len(self.expected_structure[level]):
                            completed_components = self.expected_structure[level][:]

                elif level == 'level2':
                    # Cas9 recording - check if cas9 instances exist
                    if file_outputs.get('cas9_instances', {}):
                        completed_count = sum(1 for exists in file_outputs['cas9_instances'].values() if exists)
                        if completed_count == len(self.expected_structure[level]):
                            completed_components = self.expected_structure[level][:]

                elif level == 'level3':
                    # Reconstruction - check if result files exist
                    if file_outputs.get('results', {}):
                        completed_components = [result_id for result_id, exists in file_outputs['results'].items() if exists]

            summary['levels'][level] = {
                'status_data': level_status,
                'expected_components': self.expected_structure[level],
                'completed_components': completed_components,
                'failed_components': failed_components
            }
        
        # Check LSF job status
        for pattern in ['master_gt_analysis', 'cas9_tier*_analysis', 'reconstruct_tier*']:
            jobs = self.get_lsf_job_status(pattern)
            if jobs:
                summary['lsf_jobs'][pattern] = jobs
        
        # Calculate overall progress based on actual files that exist
        file_outputs = summary['file_outputs']
        num_gt_instances = self.config.get('execution', {}).get('num_gt_instances', 1)
        cas9_simulations_per_gt = self.config.get('execution', {}).get('cas9_simulations_per_gt', 1)

        # Count completed items based on files that exist
        completed_count = 0
        total_expected = 0

        # GT trees (1 per instance)
        total_expected += num_gt_instances
        completed_count += sum(1 for exists in file_outputs['gt_trees'].values() if exists)

        # CAS9 instances (active_tiers × num_instances × cas9_simulations_per_gt)
        active_tiers = list(self.config.get('cas9_tiers', {1: {}, 2: {}, 3: {}, 4: {}}).keys())
        total_expected += len(active_tiers) * num_gt_instances * cas9_simulations_per_gt
        completed_count += sum(1 for exists in file_outputs['cas9_instances'].values() if exists)

        # Results (dynamic based on active solvers, tiers, instances, and simulations)
        active_solvers = self.config.get('solvers', {}).get('enabled', ['nj', 'maxcut', 'greedy', 'vanilla'])
        expected_results = len(active_tiers) * len(active_solvers) * num_gt_instances * cas9_simulations_per_gt
        total_expected += expected_results
        completed_count += sum(1 for exists in file_outputs['results'].values() if exists)

        # Count failures detected from logs (temporary workaround to get failure count)
        detected_failures = []
        logs_dir = self.shared_dir / "logs"
        if logs_dir.exists():
            failure_count = 0
            for err_file in logs_dir.glob("*.err"):
                try:
                    with open(err_file, 'r') as f:
                        content = f.read()
                    if any(pattern in content for pattern in ['TERM_MEMLIMIT', 'User defined signal 2', 'job killed after reaching LSF memory usage limit']):
                        failure_count += 1
                except:
                    continue
            detected_failures.append(failure_count)

        total_failed = sum(detected_failures) if detected_failures else 0

        summary['overall_progress'] = {
            'total_expected': total_expected,
            'total_completed': completed_count,
            'total_failed': total_failed,
            'completion_percentage': (completed_count / total_expected * 100) if total_expected > 0 else 0,
            'failure_percentage': (total_failed / total_expected * 100) if total_expected > 0 else 0.0
        }
        
        return summary
    
    def print_status_report(self, summary: Dict[str, Any]) -> None:
        """Print a human-readable status report."""
        print("\n" + "=" * 70)
        print("📊 LSF CASCADING JOB SYSTEM STATUS REPORT")
        print("=" * 70)
        print(f"🕐 Timestamp: {summary['timestamp']}")
        
        # Overall progress
        progress = summary['overall_progress']
        print(f"\n📈 OVERALL PROGRESS")
        print(f"   Completed: {progress['total_completed']}/{progress['total_expected']} "
              f"({progress['completion_percentage']:.1f}%)")
        print(f"   Failed:    {progress['total_failed']}/{progress['total_expected']} "
              f"({progress['failure_percentage']:.1f}%)")
        
        # Level-by-level status
        level_names = {
            'level1': '🌳 Level 1: GT Generation',
            'level2': '🧬 Level 2: Cas9 Recording',
            'level3': '🔧 Level 3: Reconstruction'
        }

        for level, level_data in summary['levels'].items():
            print(f"\n{level_names.get(level, level.upper())}")
            completed = len(level_data['completed_components'])
            expected = len(level_data['expected_components'])
            failed = len(level_data['failed_components'])

            print(f"   ✅ Completed: {completed}/{expected}")
            if failed > 0:
                print(f"   ❌ Failed: {failed}/{expected}")
                print(f"      Failed components: {', '.join(level_data['failed_components'])}")

        # Show detected failures from logs
        failures = self.detect_failures()
        if failures:
            print(f"\n🚨 FAILURE ANALYSIS")

            # Group failures by type
            failure_groups = {}
            for failure in failures:
                failure_type = failure.get('failure_type', 'unknown')
                if failure_type not in failure_groups:
                    failure_groups[failure_type] = []
                failure_groups[failure_type].append(failure)

            for failure_type, group_failures in failure_groups.items():
                print(f"   {failure_type.upper()}: {len(group_failures)} failures")

                # Show breakdown by tier/solver for reconstruction failures
                if failure_type in ['memory_limit', 'signal_termination']:
                    tier_solver_count = {}
                    for failure in group_failures:
                        if 'job_info' in failure and failure['job_info']:
                            info = failure['job_info']
                            if 'tier' in info and 'solver' in info:
                                key = f"tier{info['tier']}_{info['solver']}"
                                tier_solver_count[key] = tier_solver_count.get(key, 0) + 1

                    if tier_solver_count:
                        sorted_items = sorted(tier_solver_count.items(), key=lambda x: x[1], reverse=True)
                        top_failures = sorted_items[:5]  # Show top 5
                        for combo, count in top_failures:
                            print(f"     • {combo}: {count}")
                        if len(sorted_items) > 5:
                            remaining = sum(count for _, count in sorted_items[5:])
                            print(f"     • others: {remaining}")
        
        # File outputs
        files = summary['file_outputs']
        print(f"\n📁 FILE OUTPUTS")
        
        # GT Trees
        gt_completed = sum(1 for exists in files['gt_trees'].values() if exists)
        gt_total = len(files['gt_trees'])
        print(f"   GT Trees: {gt_completed}/{gt_total}")
        
        # Cas9 Instances
        cas9_completed = sum(1 for exists in files['cas9_instances'].values() if exists)
        cas9_total = len(files['cas9_instances'])
        print(f"   Cas9 Instances: {cas9_completed}/{cas9_total}")
        
        # Results
        results_completed = sum(1 for exists in files['results'].values() if exists)
        results_total = len(files['results'])
        print(f"   Result Files: {results_completed}/{results_total}")
        
        # LSF job status
        if summary['lsf_jobs']:
            print(f"\n⚙️  LSF JOB STATUS")
            for pattern, jobs in summary['lsf_jobs'].items():
                if jobs:
                    statuses = {}
                    for job in jobs:
                        status = job['status']
                        statuses[status] = statuses.get(status, 0) + 1
                    
                    status_str = ', '.join(f"{status}: {count}" for status, count in statuses.items())
                    print(f"   {pattern}: {status_str}")
        
        # Add bjobs snapshot with fixed height
        print(f"\n💻 BJOBS SNAPSHOT")
        try:
            result = subprocess.run(['bjobs'], capture_output=True, text=True)
            if result.returncode == 0 and result.stdout.strip():
                lines = result.stdout.strip().split('\n')
                # Always show header + up to 10 job lines for consistent height
                print(lines[0])  # Header
                job_lines = lines[1:11]  # Up to 10 jobs
                for line in job_lines:
                    print(line)
                # Pad with empty lines to maintain height of 11 lines total
                for _ in range(11 - len(lines)):
                    print("")
                if len(lines) > 11:
                    print(f"   ... and {len(lines) - 11} more jobs (showing first 10)")
            else:
                print("JOBID      USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME")
                for _ in range(10):
                    print("")
                print("   No active jobs")
        except Exception as e:
            print("JOBID      USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME")
            for _ in range(10):
                print("")
            print(f"   Error getting bjobs: {e}")
        
        print("=" * 70 + "\n")
    
    
    def save_status_summary(self, summary: Dict[str, Any]) -> None:
        """Save status summary to file."""
        summary_file = self.status_dir / "overall_status.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        logger.info(f"Status summary saved to: {summary_file}")
    
    def monitor_once(self) -> Dict[str, Any]:
        """Run monitoring once and return status."""
        logger.info("Checking job status...")
        summary = self.generate_status_summary()
        self.print_status_report(summary)
        self.save_status_summary(summary)
        return summary
    
    def monitor_continuous(self, interval_seconds: int = 60) -> None:
        """Monitor continuously with specified interval."""
        logger.info(f"Starting continuous monitoring (interval: {interval_seconds}s)")
        logger.info("Press Ctrl+C to stop monitoring")
        
        # Set up signal handler for graceful shutdown
        def signal_handler(signum, frame):
            logger.info("Monitoring stopped by user")
            sys.exit(0)
        
        signal.signal(signal.SIGINT, signal_handler)
        signal.signal(signal.SIGTERM, signal_handler)
        
        try:
            while True:
                self.monitor_once()
                logger.info(f"Sleeping for {interval_seconds} seconds...")
                time.sleep(interval_seconds)
                
        except KeyboardInterrupt:
            logger.info("Monitoring stopped by user")
            sys.exit(0)
    
    def detect_failures(self) -> List[Dict[str, Any]]:
        """Detect and report job failures."""
        failures = []
        summary = self.generate_status_summary()

        # Check status file failures
        for level, level_data in summary['levels'].items():
            for failed_component in level_data['failed_components']:
                component_data = level_data['status_data'][failed_component]
                failures.append({
                    'level': level,
                    'component': failed_component,
                    'error': component_data.get('details'),
                    'timestamp': component_data.get('timestamp'),
                    'failure_type': 'status_file'
                })

        # Check log files for memory limit failures and other LSF failures
        log_failures = self.detect_log_failures()
        failures.extend(log_failures)

        return failures

    def detect_log_failures(self) -> List[Dict[str, Any]]:
        """Detect failures from LSF log files."""
        failures = []

        # Check logs directory for error files
        logs_dir = self.shared_dir / "logs"
        if not logs_dir.exists():
            return failures

        for err_file in logs_dir.glob("*.err"):
            try:
                with open(err_file, 'r') as f:
                    content = f.read()

                # Check for various LSF failure patterns
                failure_patterns = {
                    'memory_limit': ['TERM_MEMLIMIT', 'job killed after reaching LSF memory usage limit'],
                    'time_limit': ['TERM_RUNLIMIT', 'job killed after reaching LSF time usage limit'],
                    'signal_termination': ['User defined signal 2', 'Terminated'],
                    'system_error': ['ERROR', 'FAILED', 'Exception', 'Traceback']
                }

                job_name = err_file.stem  # Remove .err extension

                # Extract job details from filename (e.g., reconstruct_instance0_sim0_recon0_tier1_maxcut_438825)
                job_info = self.parse_job_name(job_name)

                for failure_type, patterns in failure_patterns.items():
                    if any(pattern in content for pattern in patterns):
                        failure_detail = self.extract_failure_detail(content, failure_type)

                        failures.append({
                            'level': job_info.get('level', 'unknown'),
                            'component': job_name,
                            'error': failure_detail,
                            'timestamp': err_file.stat().st_mtime,
                            'failure_type': failure_type,
                            'job_info': job_info
                        })
                        break  # Only report first failure type found

            except Exception as e:
                logger.debug(f"Error reading log file {err_file}: {e}")
                continue

        return failures

    def parse_job_name(self, job_name: str) -> Dict[str, Any]:
        """Parse job name to extract job information."""
        info = {'level': 'unknown'}

        if job_name.startswith('reconstruct_'):
            info['level'] = 'level3'
            # Try to extract: reconstruct_instance0_sim0_recon0_tier1_maxcut_438825
            parts = job_name.split('_')
            if len(parts) >= 6:
                try:
                    info['instance'] = int(parts[1].replace('instance', ''))
                    info['simulation'] = int(parts[2].replace('sim', ''))
                    info['reconstruction'] = int(parts[3].replace('recon', ''))
                    info['tier'] = int(parts[4].replace('tier', ''))
                    info['solver'] = parts[5]
                except (ValueError, IndexError):
                    pass
        elif job_name.startswith('cas9_'):
            info['level'] = 'level2'
            # Try to extract: cas9_instance0_sim0_tier1_437687
            parts = job_name.split('_')
            if len(parts) >= 4:
                try:
                    info['instance'] = int(parts[1].replace('instance', ''))
                    info['simulation'] = int(parts[2].replace('sim', ''))
                    info['tier'] = int(parts[3].replace('tier', ''))
                except (ValueError, IndexError):
                    pass
        elif 'gt' in job_name.lower():
            info['level'] = 'level1'

        return info

    def extract_failure_detail(self, content: str, failure_type: str) -> str:
        """Extract detailed failure information from log content."""
        if failure_type == 'memory_limit':
            if 'TERM_MEMLIMIT' in content:
                return "Job killed after reaching LSF memory usage limit"
        elif failure_type == 'time_limit':
            if 'TERM_RUNLIMIT' in content:
                return "Job killed after reaching LSF time usage limit"
        elif failure_type == 'signal_termination':
            if 'User defined signal 2' in content:
                return "Job terminated with signal 2 (likely system kill)"
        elif failure_type == 'system_error':
            # Try to extract the actual error message
            lines = content.split('\n')
            for i, line in enumerate(lines):
                if any(keyword in line for keyword in ['ERROR', 'Exception', 'Traceback']):
                    # Return this line and potentially the next few lines for context
                    error_lines = lines[i:i+3]
                    return ' | '.join(error_lines)

        return f"Detected {failure_type} failure"


def main():
    parser = argparse.ArgumentParser(description="Job Monitor - Track LSF cascading job system status")
    parser.add_argument('--shared_dir', required=True,
                       help='Shared directory for job coordination')
    parser.add_argument('--continuous', action='store_true',
                       help='Run continuous monitoring')
    parser.add_argument('--interval', type=int, default=60,
                       help='Monitoring interval in seconds (default: 60)')
    parser.add_argument('--failures_only', action='store_true',
                       help='Only report failures')
    
    args = parser.parse_args()
    
    monitor = JobMonitor(args.shared_dir)
    
    if args.failures_only:
        failures = monitor.detect_failures()
        if failures:
            print("❌ DETECTED FAILURES:")

            # Group failures by type
            failure_groups = {}
            for failure in failures:
                failure_type = failure.get('failure_type', 'unknown')
                if failure_type not in failure_groups:
                    failure_groups[failure_type] = []
                failure_groups[failure_type].append(failure)

            for failure_type, group_failures in failure_groups.items():
                print(f"\n   {failure_type.upper()} ({len(group_failures)} failures):")
                for failure in group_failures[:10]:  # Show first 10 per type
                    component = failure['component']
                    error = failure['error']
                    if 'job_info' in failure and failure['job_info']:
                        info = failure['job_info']
                        if 'tier' in info and 'solver' in info:
                            component = f"tier{info['tier']}_{info['solver']}"
                    print(f"     • {component}: {error}")
                if len(group_failures) > 10:
                    print(f"     ... and {len(group_failures) - 10} more")
        else:
            print("✅ No failures detected")
    elif args.continuous:
        monitor.monitor_continuous(args.interval)
    else:
        monitor.monitor_once()


if __name__ == "__main__":
    main()
