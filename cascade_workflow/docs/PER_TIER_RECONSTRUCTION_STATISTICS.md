# Per-Tier Reconstruction Statistics Enhancement

## Overview

This enhancement addresses a critical flaw in the original cascade workflow system by enabling **per-tier reconstruction statistics**. Previously, the system only supported a global `reconstructions_per_solver` setting, preventing statistical analysis of the same GT-tier combination.

## Problem Solved

**Before**: All tiers had the same number of reconstructions per solver, regardless of tier quality or statistical needs.

**After**: Each tier can specify its own number of reconstructions per solver, enabling:
- Robust statistics for high-quality tiers (e.g., 5-10 reconstructions)
- Basic statistics for medium-quality tiers (e.g., 2-3 reconstructions)
- Single reconstructions for low-quality tiers (e.g., 1 reconstruction)

## Statistical Benefits

### Confidence Intervals
With multiple reconstructions per tier, you can now compute:
- Mean reconstruction quality metrics
- Standard deviations
- Confidence intervals (e.g., 95% CI for RF distance)

### Variance Analysis
- Assess solver stability across repeated runs
- Identify which solvers are most consistent
- Quantify reconstruction variability by tier quality

### Statistical Significance Testing
- Compare solvers with proper statistical power
- Test hypotheses about reconstruction quality differences
- Validate performance claims with p-values

## Configuration Syntax

### Per-Tier Settings (New)
```yaml
cas9_tiers:
  1:
    name: "Tier 1 - Ultra High Fidelity"
    # ... existing tier config ...
    reconstructions_per_solver: 5  # Tier-specific setting
    
  2:
    name: "Tier 2 - High Fidelity"  
    # ... existing tier config ...
    reconstructions_per_solver: 3  # Different per tier
    
  3:
    name: "Tier 3 - Medium Fidelity"
    # ... existing tier config ...
    reconstructions_per_solver: 2
    
  4:
    name: "Tier 4 - Low Fidelity"
    # ... existing tier config ...
    reconstructions_per_solver: 1

solvers:
  reconstructions_per_solver: 1  # Global fallback
```

### Backward Compatibility
Old configurations without per-tier settings still work:
```yaml
cas9_tiers:
  1:
    name: "Tier 1"
    # ... tier config (no reconstructions_per_solver) ...

solvers:
  reconstructions_per_solver: 2  # Used for all tiers
```

## Implementation Details

### Code Changes
1. **`Cas9TierConfig` dataclass**: Added optional `reconstructions_per_solver` field
2. **`get_reconstructions_for_tier()`**: New method to get tier-specific or fallback count
3. **`total_expected_jobs()`**: Updated to calculate per-tier job counts
4. **Configuration validation**: Maintains backward compatibility

### Job Calculation Formula
**Before**: `num_gt × num_tiers × num_solvers × global_reconstructions`

**After**: `num_gt × Σ(num_solvers × tier_reconstructions)` for each tier

### Example Impact
For a configuration with:
- 10 GT instances
- 4 tiers with reconstructions [5, 3, 2, 1]  
- 5 enabled solvers

**Jobs**: 10 GT + 40 Cas9 + 10×5×(5+3+2+1) = 600 total jobs

## Statistical Power Examples

### High-Quality Tier (Tier 1): 5 reconstructions × 5 solvers = 25 per GT
- **Excellent** statistical power for confidence intervals
- Can detect small differences between solvers
- Robust mean and variance estimates

### Medium-Quality Tier (Tier 2): 3 reconstructions × 5 solvers = 15 per GT  
- **Good** statistical power for comparisons
- Basic confidence intervals possible
- Reasonable variance estimates

### Low-Quality Tier (Tier 4): 1 reconstruction × 5 solvers = 5 per GT
- **Minimal** statistical power
- Point estimates only, no confidence intervals
- Suitable for baseline comparisons

## Configuration Files

### Production Example
- **`cascade_config_medium.yaml`**: Updated with per-tier settings
- Tier 1: 5 reconstructions (high precision)
- Tier 2: 3 reconstructions (good stats)  
- Tier 3: 2 reconstructions (basic stats)
- Tier 4: 1 reconstruction (baseline)

### Statistical Test Example  
- **`cascade_config_statistical_test.yaml`**: Demonstrates statistical focus
- Tier 1: 10 reconstructions (maximum precision)
- Tier 2: 5 reconstructions (good precision)
- Tier 3: 3 reconstructions (basic precision)

## Best Practices

### Reconstruction Count Guidelines
- **High-fidelity tiers**: 5-10 reconstructions for robust statistics
- **Medium-fidelity tiers**: 3-5 reconstructions for good comparisons
- **Low-fidelity tiers**: 1-2 reconstructions for basic assessment

### Statistical Analysis
- Always report confidence intervals when available
- Use appropriate statistical tests for comparisons
- Account for multiple testing corrections
- Validate assumptions (normality, independence)

### Computational Considerations
- Higher reconstruction counts increase job count dramatically
- Balance statistical power with computational resources
- Consider tier quality when allocating reconstructions

## Testing

The enhancement has been thoroughly tested:
- ✅ Per-tier reconstruction counts work correctly
- ✅ Backward compatibility maintained
- ✅ Job counting calculation accurate
- ✅ Configuration validation robust

## Impact

This enhancement transforms the cascade workflow from a simple reconstruction system into a **statistically rigorous benchmarking platform**, enabling:
- Quantitative solver comparisons with confidence intervals
- Tier-specific statistical analysis
- Publication-quality reconstruction quality assessments
- Robust performance validation across different recording qualities