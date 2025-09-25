# Local Test Infrastructure for Fast Metrics Testing

This directory contains a complete mini-pipeline for testing metrics computation fixes without the overhead of LSF job submission.

## Overview

**Problem:** The full cascade pipeline takes hours to run, making it slow to test metrics fixes.
**Solution:** Local mini-pipeline that runs the same production code but locally, with smart caching.

## Files

- **`test_config_mini.yaml`** - Small-scale configuration (N=50 cells, 2 GT instances, 2 tiers, 3 solvers = 24 tests)
- **`local_test_pipeline.py`** - Main pipeline runner with smart caching
- **`test_metrics_validation.py`** - Comprehensive validation suite for metrics fixes
- **`cache/`** - Cached GT trees and CAS9 instances (auto-generated)
- **`results/`** - Test outputs and validation reports

## Usage

### Basic Testing
```bash
# Run the mini-pipeline (1-2 minutes vs hours)
cd cascade_workflow/src/testing
python local_test_pipeline.py --config test_config_mini.yaml
```

### Comprehensive Validation
```bash
# Test all metrics fixes with validation suite
python test_metrics_validation.py

# Generate detailed report
python test_metrics_validation.py --config test_config_mini.yaml
```

### Custom Testing
```bash
# Use custom config
python local_test_pipeline.py --config your_config.yaml

# Custom cache directory
python local_test_pipeline.py --config test_config_mini.yaml --cache-dir /tmp/test_cache
```

## What It Tests

### 1. Seed Determinism ✅
- Same config → identical results
- Verifies our seed fixes work

### 2. Likelihood Stability ✅
- No more NaN values
- No more -150K artificial values
- Reasonable likelihood ranges

### 3. Parameter Validation ✅
- Missing parameters → proper errors
- Required parameters are checked

### 4. Config Parameter Usage ✅
- `analysis.reconstruction_seed` is used
- `analysis.triplets_trials` affects results
- Branch length parameters are respected

### 5. No Hardcoded Values ✅
- No hardcoded `seed=42`
- No hardcoded `trials=1000`
- Results show appropriate variation

### 6. Lambda Computation Methods ✅
- Both `lam_simulation` and `lam_gt` are calculated
- Values are in reasonable ranges (0 < lam < 10)
- Methods show appropriate correlation but use different data sources
- Source tracking (`phs_lam_source`) works correctly

## Architecture

### Smart Caching System
- **GT Trees:** Cached by `ground_truth` + `random_seed` config hash
- **CAS9 Instances:** Cached by `cas9_tiers` + simulations config hash
- **Reconstruction:** Always fresh (this is what we're testing)

### Production Code Reuse
- Uses actual `GTWorker`, `Cas9RecordingWorker`, `reconstruct_and_calculate_metrics`
- No mocks or simplified versions
- Realistic testing environment

### Fast Execution
- **Full Pipeline:** ~50 cells vs 1000+ cells (20x faster)
- **Smart Caching:** Reuses expensive GT/CAS9 generation
- **Local Execution:** No LSF job queues
- **Expected Runtime:** 1-2 minutes vs hours

## Test Scenarios

### Typical Test Run
- 2 GT instances × 2 simulations × 2 tiers × 3 solvers = **24 reconstructions**
- Tests all major code paths
- Validates all metrics fixes
- Generates comprehensive reports

### Before/After Comparison
- Run with old hardcoded values
- Run with new configurable parameters
- Compare results to verify improvements

## Output

### CSV Results
```
test_results_20250925_131045.csv
- Full metrics for each reconstruction
- Likelihood values, RF distances, etc.
- Metadata for analysis
```

### Validation Report
```markdown
# Metrics Validation Report
Results: 5/5 tests passed
✅ ALL TESTS PASSED

- seed_determinism: ✅ PASS
- likelihood_stability: ✅ PASS
- parameter_validation: ✅ PASS
- config_parameter_usage: ✅ PASS
- no_hardcoded_values: ✅ PASS
```

## Integration with Main Pipeline

This test infrastructure uses the **exact same code** as production:
- Same worker functions
- Same config format
- Same algorithms
- Just bypasses LSF for local execution

**→ Results are directly applicable to production pipeline!**