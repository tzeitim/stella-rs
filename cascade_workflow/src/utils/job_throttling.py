"""
Job Throttling Utilities for LSF Cascade Workflow

Provides mechanisms to control job submission rate and concurrency to avoid
overwhelming the cluster with too many simultaneous jobs.
"""

import subprocess
import time
import logging
import json
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class ThrottlingConfig:
    """Configuration for job throttling"""
    max_concurrent_cas9_jobs: int = 50
    max_concurrent_reconstruction_jobs: int = 200
    job_submission_delay: float = 1.0
    batch_size: int = 20
    check_interval: int = 30

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            'max_concurrent_cas9_jobs': self.max_concurrent_cas9_jobs,
            'max_concurrent_reconstruction_jobs': self.max_concurrent_reconstruction_jobs,
            'job_submission_delay': self.job_submission_delay,
            'batch_size': self.batch_size,
            'check_interval': self.check_interval
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ThrottlingConfig':
        """Create from dictionary"""
        return cls(
            max_concurrent_cas9_jobs=data.get('max_concurrent_cas9_jobs', 50),
            max_concurrent_reconstruction_jobs=data.get('max_concurrent_reconstruction_jobs', 200),
            job_submission_delay=data.get('job_submission_delay', 1.0),
            batch_size=data.get('batch_size', 20),
            check_interval=data.get('check_interval', 30)
        )


class JobThrottler:
    """
    Manages job submission throttling for LSF cascading workflow.

    Prevents cluster overload by:
    - Limiting concurrent job counts
    - Adding delays between submissions
    - Submitting jobs in controlled batches
    - Supporting dynamic configuration updates
    """

    def __init__(self, config: ThrottlingConfig, shared_dir: Optional[Path] = None):
        self.config = config
        self.shared_dir = shared_dir
        self.last_check_time = 0
        self.cached_job_counts = {"cas9": 0, "reconstruction": 0}
        self.last_config_check = 0
        self.config_check_interval = 300  # Check for config updates every 5 minutes

        # Set up dynamic config file path
        if self.shared_dir:
            self.dynamic_config_file = self.shared_dir / "throttling_config.yaml"
            self._initialize_dynamic_config()
        else:
            self.dynamic_config_file = None

        logger.info(f"Initialized job throttler:")
        logger.info(f"  Max CAS9 jobs: {config.max_concurrent_cas9_jobs}")
        logger.info(f"  Max reconstruction jobs: {config.max_concurrent_reconstruction_jobs}")
        logger.info(f"  Submission delay: {config.job_submission_delay}s")
        logger.info(f"  Batch size: {config.batch_size}")
        if self.dynamic_config_file:
            logger.info(f"  Dynamic config file: {self.dynamic_config_file}")

    def _initialize_dynamic_config(self):
        """Initialize the dynamic configuration file if it doesn't exist"""
        if not self.dynamic_config_file.exists():
            try:
                config_data = {
                    'throttling': self.config.to_dict(),
                    'metadata': {
                        'created_at': time.strftime('%Y-%m-%d %H:%M:%S'),
                        'description': 'Dynamic throttling configuration - edit this file to adjust job limits in real-time'
                    }
                }

                with open(self.dynamic_config_file, 'w') as f:
                    yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)

                logger.info(f"Created dynamic throttling config file: {self.dynamic_config_file}")
            except Exception as e:
                logger.warning(f"Failed to create dynamic config file: {e}")

    def _update_config_from_file(self):
        """Update configuration from dynamic file if it has been modified"""
        if not self.dynamic_config_file or not self.dynamic_config_file.exists():
            return

        try:
            current_time = time.time()

            # Check if enough time has passed since last config check
            if current_time - self.last_config_check < self.config_check_interval:
                return

            self.last_config_check = current_time

            # Read the dynamic config file
            with open(self.dynamic_config_file, 'r') as f:
                dynamic_data = yaml.safe_load(f)

            if not dynamic_data or 'throttling' not in dynamic_data:
                return

            # Update configuration
            new_config = ThrottlingConfig.from_dict(dynamic_data['throttling'])

            # Check if any values actually changed
            if (new_config.max_concurrent_cas9_jobs != self.config.max_concurrent_cas9_jobs or
                new_config.max_concurrent_reconstruction_jobs != self.config.max_concurrent_reconstruction_jobs or
                new_config.job_submission_delay != self.config.job_submission_delay or
                new_config.batch_size != self.config.batch_size or
                new_config.check_interval != self.config.check_interval):

                logger.info("Dynamic throttling configuration updated:")
                logger.info(f"  CAS9 jobs: {self.config.max_concurrent_cas9_jobs} → {new_config.max_concurrent_cas9_jobs}")
                logger.info(f"  Reconstruction jobs: {self.config.max_concurrent_reconstruction_jobs} → {new_config.max_concurrent_reconstruction_jobs}")
                logger.info(f"  Submission delay: {self.config.job_submission_delay}s → {new_config.job_submission_delay}s")
                logger.info(f"  Batch size: {self.config.batch_size} → {new_config.batch_size}")
                logger.info(f"  Check interval: {self.config.check_interval}s → {new_config.check_interval}s")

                self.config = new_config

        except Exception as e:
            logger.warning(f"Failed to update config from dynamic file: {e}")

    def get_running_job_counts(self) -> Dict[str, int]:
        """Get current running job counts by type"""
        current_time = time.time()

        # Use cached counts if checked recently
        if current_time - self.last_check_time < self.config.check_interval:
            return self.cached_job_counts.copy()

        try:
            # Query LSF for running jobs
            result = subprocess.run(
                ['bjobs', '-w', '-u', 'all', '-r'],  # Running jobs only
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode != 0:
                logger.warning(f"bjobs query failed: {result.stderr}")
                return self.cached_job_counts.copy()

            # Parse bjobs output
            lines = result.stdout.strip().split('\n')[1:]  # Skip header
            cas9_count = 0
            reconstruction_count = 0

            for line in lines:
                if not line.strip():
                    continue

                fields = line.split()
                if len(fields) >= 7:
                    job_name = fields[6]  # JOB_NAME field

                    if 'cas9_' in job_name.lower():
                        cas9_count += 1
                    elif 'reconstruct' in job_name.lower():
                        reconstruction_count += 1

            self.cached_job_counts = {
                "cas9": cas9_count,
                "reconstruction": reconstruction_count
            }
            self.last_check_time = current_time

            logger.debug(f"Current running jobs - CAS9: {cas9_count}, Reconstruction: {reconstruction_count}")

        except subprocess.TimeoutExpired:
            logger.warning("bjobs query timed out, using cached counts")
        except Exception as e:
            logger.warning(f"Failed to query job counts: {e}, using cached counts")

        return self.cached_job_counts.copy()

    def can_submit_cas9_job(self) -> bool:
        """Check if we can submit a CAS9 recording job"""
        # Update config from dynamic file first
        self._update_config_from_file()

        counts = self.get_running_job_counts()
        can_submit = counts["cas9"] < self.config.max_concurrent_cas9_jobs

        if not can_submit:
            logger.info(f"CAS9 job submission throttled: {counts['cas9']}/{self.config.max_concurrent_cas9_jobs} running")

        return can_submit

    def can_submit_reconstruction_job(self) -> bool:
        """Check if we can submit a reconstruction job"""
        # Update config from dynamic file first
        self._update_config_from_file()

        counts = self.get_running_job_counts()
        can_submit = counts["reconstruction"] < self.config.max_concurrent_reconstruction_jobs

        if not can_submit:
            logger.info(f"Reconstruction job submission throttled: {counts['reconstruction']}/{self.config.max_concurrent_reconstruction_jobs} running")

        return can_submit

    def wait_for_cas9_slot(self, timeout: int = 3600) -> bool:
        """
        Wait for a CAS9 job slot to become available

        Args:
            timeout: Maximum time to wait in seconds

        Returns:
            True if slot became available, False if timed out
        """
        start_time = time.time()

        while time.time() - start_time < timeout:
            if self.can_submit_cas9_job():
                return True

            logger.info(f"Waiting for CAS9 job slot... ({int(time.time() - start_time)}s elapsed)")
            time.sleep(self.config.check_interval)

        logger.warning(f"Timed out waiting for CAS9 job slot after {timeout}s")
        return False

    def wait_for_reconstruction_slot(self, timeout: int = 1800) -> bool:
        """
        Wait for a reconstruction job slot to become available

        Args:
            timeout: Maximum time to wait in seconds

        Returns:
            True if slot became available, False if timed out
        """
        start_time = time.time()

        while time.time() - start_time < timeout:
            if self.can_submit_reconstruction_job():
                return True

            logger.info(f"Waiting for reconstruction job slot... ({int(time.time() - start_time)}s elapsed)")
            time.sleep(self.config.check_interval)

        logger.warning(f"Timed out waiting for reconstruction job slot after {timeout}s")
        return False

    def submit_with_throttling(self, submit_func, job_type: str, *args, **kwargs):
        """
        Submit a job with throttling applied

        Args:
            submit_func: Function to call for job submission
            job_type: Either 'cas9' or 'reconstruction'
            *args, **kwargs: Arguments to pass to submit_func

        Returns:
            Result from submit_func or None if throttled
        """
        # Check if we can submit
        if job_type == 'cas9':
            if not self.wait_for_cas9_slot():
                logger.error("Failed to get CAS9 job slot, skipping submission")
                return None
        elif job_type == 'reconstruction':
            if not self.wait_for_reconstruction_slot():
                logger.error("Failed to get reconstruction job slot, skipping submission")
                return None
        else:
            raise ValueError(f"Unknown job type: {job_type}")

        # Add submission delay
        if self.config.job_submission_delay > 0:
            time.sleep(self.config.job_submission_delay)

        # Submit the job
        try:
            result = submit_func(*args, **kwargs)
            logger.debug(f"Successfully submitted {job_type} job")
            return result
        except Exception as e:
            logger.error(f"Failed to submit {job_type} job: {e}")
            raise


def create_throttler_from_config(config: Dict[str, Any], shared_dir: Optional[Path] = None) -> JobThrottler:
    """Create JobThrottler from cascade configuration"""
    throttling_config = config.get('lsf', {}).get('throttling', {})

    return JobThrottler(
        ThrottlingConfig(
            max_concurrent_cas9_jobs=throttling_config.get('max_concurrent_cas9_jobs', 50),
            max_concurrent_reconstruction_jobs=throttling_config.get('max_concurrent_reconstruction_jobs', 200),
            job_submission_delay=throttling_config.get('job_submission_delay', 1.0),
            batch_size=throttling_config.get('batch_size', 20),
            check_interval=throttling_config.get('check_interval', 30)
        ),
        shared_dir=shared_dir
    )


# Example usage
if __name__ == "__main__":
    # Test the throttling mechanism
    test_config = ThrottlingConfig(
        max_concurrent_cas9_jobs=5,
        max_concurrent_reconstruction_jobs=10,
        job_submission_delay=0.5,
        batch_size=3,
        check_interval=5
    )

    throttler = JobThrottler(test_config)

    # Test job count checking
    counts = throttler.get_running_job_counts()
    print(f"Current job counts: {counts}")

    # Test submission checks
    print(f"Can submit CAS9 job: {throttler.can_submit_cas9_job()}")
    print(f"Can submit reconstruction job: {throttler.can_submit_reconstruction_job()}")