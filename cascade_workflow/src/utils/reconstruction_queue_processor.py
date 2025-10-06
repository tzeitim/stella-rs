#!/usr/bin/env python3
"""
Reconstruction Queue Processor

Processes queued reconstruction jobs from the reconstruction_queue.jsonl file.
This allows Cas9 jobs to complete quickly by queueing reconstruction jobs
instead of waiting for submission slots.
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Any, List, Optional

# Add utils directory to path
utils_dir = Path(__file__).parent
sys.path.append(str(utils_dir))

from job_throttling import JobThrottler, ThrottlingConfig

logger = logging.getLogger(__name__)


class ReconstructionQueueProcessor:
    """Processes queued reconstruction jobs with proper throttling."""

    def __init__(self, shared_dir: Path, throttling_config: Optional[ThrottlingConfig] = None):
        self.shared_dir = Path(shared_dir)
        self.queue_file = self.shared_dir / "reconstruction_queue.jsonl"
        self.processed_file = self.shared_dir / "reconstruction_queue_processed.jsonl"
        self.failed_file = self.shared_dir / "reconstruction_queue_failed.jsonl"

        # Load configuration
        self.config = self.load_config()

        # Use provided throttling config or create default
        if throttling_config is None:
            throttling_config = ThrottlingConfig(
                max_concurrent_reconstruction_jobs=50,
                job_submission_delay=2.0,
                batch_size=10
            )

        self.throttler = JobThrottler(throttling_config, self.shared_dir)

        # Update dynamic throttling config file with command-line values
        self._update_throttling_config_file(throttling_config)

        # Ensure directories exist
        self.shared_dir.mkdir(parents=True, exist_ok=True)

    def load_config(self) -> Dict[str, Any]:
        """Load configuration from cascade_config.yaml."""
        import yaml
        config_path = self.shared_dir / "cascade_config.yaml"
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            logger.warning(f"Could not load config {config_path}: {e}")
            return {}

    def _update_throttling_config_file(self, throttling_config: ThrottlingConfig) -> None:
        """Update the dynamic throttling config file with command-line values."""
        import yaml

        throttling_file = self.shared_dir / "throttling_config.yaml"

        try:
            config_data = {
                'throttling': throttling_config.to_dict(),
                'metadata': {
                    'last_updated': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'updated_by': f'reconstruction_queue_processor.py (PID: {os.getpid()})',
                    'description': 'Dynamic throttling configuration - values set from command-line args'
                }
            }

            with open(throttling_file, 'w') as f:
                yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)

            logger.info(f"Updated throttling config file with command-line values: {throttling_file}")
        except Exception as e:
            logger.warning(f"Failed to update throttling config file: {e}")

    def load_queue(self) -> List[Dict[str, Any]]:
        """Load all queued jobs from the queue file."""
        if not self.queue_file.exists():
            logger.info(f"Queue file {self.queue_file} does not exist")
            return []

        jobs = []
        try:
            with open(self.queue_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        job = json.loads(line)
                        job['queue_line'] = line_num
                        jobs.append(job)
                    except json.JSONDecodeError as e:
                        logger.error(f"Invalid JSON on line {line_num}: {e}")

            logger.info(f"Loaded {len(jobs)} jobs from queue")
            return jobs

        except Exception as e:
            logger.error(f"Failed to load queue file {self.queue_file}: {e}")
            return []

    def should_submit_job(self, job: Dict[str, Any]) -> bool:
        """Check if a job should be submitted (not already processed).

        Checks multiple sources to avoid duplicate submissions:
        1. Processed jobs file
        2. Partitioned results (parquet files)
        3. Legacy JSON results files
        4. Cas9 instance file existence
        """

        # Check if job was already submitted (in processed file)
        if hasattr(self, '_processed_jobs'):
            job_key = (job['cas9_instance_path'], job['solver'], job['tier'], job['instance_id'], job['cas9_simulation_id'])
            if job_key in self._processed_jobs:
                logger.debug(f"Job already processed: {job['solver']} for instance {job['instance_id']}_sim{job['cas9_simulation_id']}")
                return False

        # Check if result already exists in partitioned results (primary storage)
        shared_dir = Path(job['shared_dir'])
        partitioned_dir = shared_dir / "partitioned_results" / f"cas9_tier={job['tier']}" / f"solver={job['solver']}"

        if partitioned_dir.exists():
            # Check parquet files for this result
            import polars as pl
            try:
                parquet_files = list(partitioned_dir.glob("*.parquet"))
                if parquet_files:
                    df = pl.scan_parquet(str(partitioned_dir / "*.parquet")).collect()
                    # Check for exact match on instance_id and cas9_simulation_id
                    matches = df.filter(
                        (pl.col('gt_instance_id') == job['instance_id']) &
                        (pl.col('cas9_simulation_id') == job['cas9_simulation_id'])
                    )
                    if len(matches) > 0:
                        logger.debug(f"Result already exists in partitioned storage: {job['solver']} for "
                                   f"instance {job['instance_id']}_sim{job['cas9_simulation_id']}")
                        return False
            except Exception as e:
                logger.debug(f"Could not check partitioned results: {e}")

        # Check if result already exists in legacy JSON format
        cas9_instance_path = Path(job['cas9_instance_path'])
        solver = job['solver']
        instance_name = cas9_instance_path.stem

        results_dir = shared_dir / "results"
        result_pattern = f"{instance_name}_*_{solver}_metrics.json"

        existing_results = list(results_dir.glob(result_pattern))
        if existing_results:
            logger.debug(f"Result already exists (legacy): {instance_name} {solver}: {existing_results[0]}")
            return False

        # Check if cas9 instance file exists
        if not cas9_instance_path.exists():
            logger.warning(f"Cas9 instance file does not exist: {cas9_instance_path}")
            return False

        return True

    def load_processed_jobs(self) -> set:
        """Load set of already processed jobs."""
        processed_jobs = set()

        if self.processed_file.exists():
            try:
                with open(self.processed_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        try:
                            job = json.loads(line)
                            job_key = (job['cas9_instance_path'], job['solver'], job['tier'], job['instance_id'], job['cas9_simulation_id'])
                            processed_jobs.add(job_key)
                        except json.JSONDecodeError:
                            continue
            except Exception as e:
                logger.warning(f"Could not load processed jobs: {e}")

        logger.info(f"Loaded {len(processed_jobs)} previously processed jobs")
        return processed_jobs

    def submit_reconstruction_job(self, job: Dict[str, Any]) -> Optional[str]:
        """Submit a single reconstruction job."""

        def _do_submit():
            """Internal function to perform the actual submission."""

            # Get LSF configuration for reconstruction jobs
            lsf_config = self.config.get('lsf', {})
            reconstruction_queue = lsf_config.get('queues', {}).get('reconstruction', 'short')
            recon_resources = lsf_config.get('resources', {}).get('reconstruction', {})

            cores = recon_resources.get('cores', 50)
            memory_gb = recon_resources.get('memory_gb', 1.5)

            # Build bsub command with configuration values
            cmd = [
                'bsub',
                '-J', f"reconstruct_instance{job['instance_id']}_sim{job['cas9_simulation_id']}_recon0_tier{job['tier']}_{job['solver']}",
                '-oo', f"{job['shared_dir']}/logs/reconstruct_instance{job['instance_id']}_sim{job['cas9_simulation_id']}_recon0_tier{job['tier']}_{job['solver']}_%J.out",
                '-eo', f"{job['shared_dir']}/logs/reconstruct_instance{job['instance_id']}_sim{job['cas9_simulation_id']}_recon0_tier{job['tier']}_{job['solver']}_%J.err",
                '-q', reconstruction_queue,
                '-n', str(cores), '-R', 'span[hosts=1]',
                '-R', f'rusage[mem={memory_gb}GB]',
                'python', f"{job['shared_dir']}/reconstruction_worker.py",
                '--cas9_instance_path', job['cas9_instance_path'],
                '--solver', job['solver'],
                '--tier', str(job['tier']),
                '--output_dir', job['output_dir'],
                '--shared_dir', job['shared_dir'],
                '--gt_instance_id', str(job['instance_id']),
                '--cas9_simulation_id', str(job['cas9_simulation_id']),
                '--reconstruction_id', '0'
            ]

            # Submit job
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode == 0:
                # Extract job ID from bsub output
                output_lines = result.stdout.strip().split('\n')
                for line in output_lines:
                    if 'Job <' in line and '> is submitted' in line:
                        job_id = line.split('<')[1].split('>')[0]
                        return job_id

                logger.warning(f"Could not extract job ID from bsub output: {result.stdout}")
                return "unknown"
            else:
                raise RuntimeError(f"bsub failed: {result.stderr}")

        try:
            # Use throttling for reconstruction job submission
            job_id = self.throttler.submit_with_throttling(_do_submit, 'reconstruction')

            if job_id:
                logger.info(f"Submitted queued job {job_id}: {job['solver']} for {Path(job['cas9_instance_path']).stem}")
                return job_id
            else:
                logger.error(f"Failed to get slot for {job['solver']} job")
                return None

        except Exception as e:
            logger.error(f"Failed to submit {job['solver']} job: {e}")
            return None

    def mark_job_processed(self, job: Dict[str, Any], job_id: Optional[str], status: str) -> None:
        """Mark a job as processed by writing to processed file."""

        processed_job = job.copy()
        processed_job.update({
            'processed_timestamp': time.time(),
            'submitted_job_id': job_id,
            'processing_status': status
        })

        target_file = self.processed_file if status == 'submitted' else self.failed_file

        with open(target_file, 'a') as f:
            f.write(json.dumps(processed_job) + '\n')

    def process_queue(self, max_jobs: Optional[int] = None, dry_run: bool = False) -> None:
        """Process the reconstruction queue."""

        jobs = self.load_queue()
        if not jobs:
            logger.info("No jobs in queue")
            return

        # Load previously processed jobs to avoid duplicates
        self._processed_jobs = self.load_processed_jobs()

        examined_count = 0
        submitted_count = 0
        skipped_count = 0
        failed_count = 0

        for job in jobs:
            examined_count += 1

            # Check if we should submit this job
            if not self.should_submit_job(job):
                logger.info(f"Skipping job {examined_count}: {job['solver']} for {Path(job['cas9_instance_path']).stem} (already exists or invalid)")
                skipped_count += 1
                continue

            # Check max_jobs limit only for jobs we actually need to process
            if max_jobs and submitted_count + failed_count >= max_jobs:
                logger.info(f"Reached maximum job submission limit ({max_jobs})")
                break

            if dry_run:
                logger.info(f"DRY RUN: Would submit {job['solver']} for {Path(job['cas9_instance_path']).stem}")
                continue

            # Submit the job
            job_id = self.submit_reconstruction_job(job)

            if job_id:
                self.mark_job_processed(job, job_id, 'submitted')
                submitted_count += 1
            else:
                self.mark_job_processed(job, None, 'failed')
                failed_count += 1

        logger.info(f"Queue processing complete: {examined_count} examined, {submitted_count} submitted, {skipped_count} skipped, {failed_count} failed")

    def cleanup_queue(self) -> None:
        """Remove processed jobs from the queue file."""

        if not self.processed_file.exists():
            logger.info("No processed jobs to clean up")
            return

        # Load processed job identifiers
        processed_jobs = set()
        with open(self.processed_file, 'r') as f:
            for line in f:
                try:
                    job = json.loads(line.strip())
                    if 'queue_line' in job:
                        processed_jobs.add(job['queue_line'])
                except json.JSONDecodeError:
                    continue

        # Rewrite queue file without processed jobs
        if not self.queue_file.exists():
            return

        temp_file = self.queue_file.with_suffix('.tmp')
        kept_count = 0

        with open(self.queue_file, 'r') as infile, open(temp_file, 'w') as outfile:
            for line_num, line in enumerate(infile, 1):
                if line_num not in processed_jobs:
                    outfile.write(line)
                    kept_count += 1

        # Replace original file
        temp_file.replace(self.queue_file)

        logger.info(f"Cleaned up queue: {len(processed_jobs)} processed jobs removed, {kept_count} jobs remain")


def main():
    """Main entry point for queue processor."""

    parser = argparse.ArgumentParser(description="Process reconstruction job queue")
    parser.add_argument('shared_dir', help='Shared directory containing the queue file')
    parser.add_argument('--max-jobs', type=int, help='Maximum number of jobs to process')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done without submitting')
    parser.add_argument('--cleanup', action='store_true', help='Clean up processed jobs from queue')
    parser.add_argument('--max-concurrent', type=int, default=50, help='Maximum concurrent reconstruction jobs')
    parser.add_argument('--delay', type=float, default=2.0, help='Delay between job submissions')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Create throttling config
    throttling_config = ThrottlingConfig(
        max_concurrent_reconstruction_jobs=args.max_concurrent,
        job_submission_delay=args.delay,
        batch_size=10
    )

    # Create processor
    processor = ReconstructionQueueProcessor(args.shared_dir, throttling_config)

    if args.cleanup:
        processor.cleanup_queue()
    else:
        processor.process_queue(max_jobs=args.max_jobs, dry_run=args.dry_run)


if __name__ == "__main__":
    main()