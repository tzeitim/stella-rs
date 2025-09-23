# Job Throttling Guide

This guide explains how to configure job throttling in the cascade workflow to avoid overwhelming your compute cluster with too many simultaneous jobs.

## Overview

The job throttling system controls:
- **Concurrent job limits**: Maximum number of jobs running simultaneously
- **Submission delays**: Pauses between job submissions
- **Batch sizes**: Number of jobs submitted together
- **Monitoring intervals**: How often to check job counts

## Configuration

Add the `throttling` section to your `cascade_config.yaml` under the `lsf` section:

```yaml
lsf:
  throttling:
    max_concurrent_cas9_jobs: 20          # Max CAS9 recording jobs
    max_concurrent_reconstruction_jobs: 100  # Max reconstruction jobs
    job_submission_delay: 2.0             # Seconds between submissions
    batch_size: 10                        # Jobs per batch
    check_interval: 30                    # Seconds between job count checks
```

## Parameter Tuning

### `max_concurrent_cas9_jobs`
Controls CAS9 recording job concurrency.

**Recommendations:**
- **Small cluster (< 100 cores)**: 5-10 jobs
- **Medium cluster (100-500 cores)**: 10-30 jobs
- **Large cluster (500+ cores)**: 30-50 jobs

CAS9 jobs are typically CPU and memory intensive, so be conservative.

### `max_concurrent_reconstruction_jobs`
Controls reconstruction job concurrency.

**Recommendations:**
- **Small cluster**: 20-50 jobs
- **Medium cluster**: 50-150 jobs
- **Large cluster**: 150-300 jobs

Reconstruction jobs are lighter than CAS9 jobs, so you can allow more.

### `job_submission_delay`
Delay in seconds between individual job submissions.

**Recommendations:**
- **Aggressive throttling**: 5.0-10.0 seconds
- **Moderate throttling**: 1.0-3.0 seconds
- **Light throttling**: 0.5-1.0 seconds

Longer delays reduce cluster scheduler load but slow down job submission.

### `batch_size`
Number of jobs to submit before checking limits again.

**Recommendations:**
- **Conservative**: 5-10 jobs per batch
- **Moderate**: 10-20 jobs per batch
- **Aggressive**: 20+ jobs per batch

Smaller batches provide better control but more overhead.

### `check_interval`
How often to query LSF for current job counts.

**Recommendations:**
- **Frequent checking**: 15-30 seconds
- **Moderate checking**: 30-60 seconds
- **Infrequent checking**: 60+ seconds

More frequent checking provides better control but increases LSF query load.

## Example Configurations

### Conservative (Cluster-Friendly)
```yaml
lsf:
  throttling:
    max_concurrent_cas9_jobs: 5
    max_concurrent_reconstruction_jobs: 25
    job_submission_delay: 10.0
    batch_size: 3
    check_interval: 60
```

**Use when:** Shared cluster, peak hours, or first-time massive runs.

### Moderate (Balanced)
```yaml
lsf:
  throttling:
    max_concurrent_cas9_jobs: 15
    max_concurrent_reconstruction_jobs: 75
    job_submission_delay: 3.0
    batch_size: 8
    check_interval: 45
```

**Use when:** Dedicated time slots or medium-sized experiments.

### Aggressive (Fast Submission)
```yaml
lsf:
  throttling:
    max_concurrent_cas9_jobs: 40
    max_concurrent_reconstruction_jobs: 200
    job_submission_delay: 1.0
    batch_size: 15
    check_interval: 30
```

**Use when:** Dedicated cluster access or off-peak hours.

## Monitoring Throttling

The workflow logs throttling actions:

```
INFO - CAS9 job submission throttled: 20/20 running
INFO - Waiting for CAS9 job slot... (30s elapsed)
INFO - Successfully submitted CAS9 job 12345 for instance 5, sim 2, tier 3
```

### Key Log Messages

- `"job submission throttled"`: Hit concurrent limit, waiting
- `"Waiting for [job type] slot"`: Actively waiting for availability
- `"Successfully submitted"`: Job submitted after throttling check
- `"Timed out waiting"`: Couldn't get slot within timeout

## Calculating Expected Jobs

For planning purposes, calculate total expected jobs:

```
CAS9 Jobs = num_gt_instances × cas9_simulations_per_gt × num_active_tiers

Reconstruction Jobs = CAS9 Jobs × num_enabled_solvers × reconstructions_per_solver
```

**Example:**
- 20 GT instances
- 10 CAS9 simulations per GT
- 4 active tiers
- 5 enabled solvers
- 1 reconstruction per solver

```
CAS9 Jobs = 20 × 10 × 4 = 800 jobs
Reconstruction Jobs = 800 × 5 × 1 = 4,000 jobs
Total = 4,800 jobs
```

## Troubleshooting

### Jobs Not Submitting
**Symptoms:** Workflow stalls, no new jobs appearing
**Causes:**
- Limits too restrictive
- LSF scheduler issues
- Connectivity problems

**Solutions:**
- Increase concurrent limits
- Check `bjobs` command works
- Verify network connectivity

### Cluster Overload
**Symptoms:** Complaints from cluster admins, slow job starts
**Causes:**
- Limits too high
- Insufficient delays
- Large batch sizes

**Solutions:**
- Reduce concurrent limits
- Increase submission delays
- Use smaller batch sizes

### Slow Progress
**Symptoms:** Jobs submitting very slowly
**Causes:**
- Too conservative settings
- Long delays
- Small batch sizes

**Solutions:**
- Increase limits gradually
- Reduce delays
- Increase batch sizes

## Best Practices

1. **Start Conservative**: Begin with low limits and increase gradually
2. **Monitor Initially**: Watch first few hours of large runs closely
3. **Peak Hours**: Use more conservative settings during busy times
4. **Communicate**: Inform cluster admins about large planned runs
5. **Test Small**: Validate throttling with small test runs first
6. **Document**: Keep notes on what settings work for your cluster

## Advanced Features

### Timeout Handling
The throttler will wait up to specified timeouts for job slots:
- CAS9 jobs: 1 hour timeout (3600s)
- Reconstruction jobs: 30 minutes timeout (1800s)

### Job Type Detection
The throttler identifies job types by parsing `bjobs` output for keywords:
- CAS9 jobs: Job names containing `"cas9_"`
- Reconstruction jobs: Job names containing `"reconstruct"`

### Caching
Job count queries are cached for the `check_interval` duration to reduce LSF load.