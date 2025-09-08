# Cas9 Simulation Statistics Enhancement

## Overview

This enhancement implements the correct statistical approach for the cascade workflow by enabling **multiple Cas9 simulation generations per ground truth**. This addresses the need for proper statistical analysis by capturing variance in Cas9 recording processes rather than just reconstruction algorithms.

## The Correct Statistical Approach

### What We Generate
**1 Ground Truth** → **N Cas9 Simulations** → **Each reconstructed by each solver**

This structure captures the natural variability in:
- Cas9 recording stochasticity
- Site-specific mutation patterns  
- Recording quality across different simulation runs
- Integration efficiency variations

### Statistical Significance
The key insight is that **recording variability is often larger than reconstruction variability**. By generating multiple Cas9 simulations from the same ground truth, we can:

1. **Assess solver robustness** across different recording realizations
2. **Quantify recording quality impact** on reconstruction success
3. **Compare solver consistency** when faced with recording noise
4. **Generate confidence intervals** for reconstruction metrics across recording conditions

## Configuration

### Global Setting
```yaml
execution:
  cas9_simulations_per_gt: 5  # Generate 5 Cas9 simulations per GT
```

### Job Structure
For each GT instance:
- **1 GT generation job**  
- **`cas9_simulations_per_gt × num_tiers`** Cas9 simulation jobs
- **`cas9_simulations_per_gt × num_tiers × num_solvers`** reconstruction jobs

### Example Impact
With `cas9_simulations_per_gt: 5`, 4 tiers, 3 solvers:
- **Per GT**: 5 × 4 × 3 = **60 reconstructions**
- **Statistical power**: 15 reconstructions per solver per tier
- **Recording variance**: 5 different Cas9 realizations per tier

## Statistical Benefits

### Recording Quality Assessment
```
GT → Cas9_Sim_1 → [Solver_1, Solver_2, Solver_3]
   → Cas9_Sim_2 → [Solver_1, Solver_2, Solver_3]  
   → Cas9_Sim_3 → [Solver_1, Solver_2, Solver_3]
   → Cas9_Sim_4 → [Solver_1, Solver_2, Solver_3]
   → Cas9_Sim_5 → [Solver_1, Solver_2, Solver_3]
```

This allows you to answer:
- **"How consistent is Solver X across different Cas9 recordings?"**
- **"Does recording quality variance dominate reconstruction variance?"**
- **"Which solver is most robust to recording noise?"**

### Confidence Intervals
With 5 Cas9 simulations per tier:
- **Mean reconstruction quality** across recording realizations
- **Standard error** accounting for recording variability
- **95% confidence intervals** for solver performance
- **Effect size estimation** for solver comparisons

### Variance Decomposition
```
Total Variance = Recording Variance + Reconstruction Variance + Interaction
```

You can now separate:
- Variance due to Cas9 recording stochasticity
- Variance due to reconstruction algorithm differences  
- Interaction effects between recording quality and solver performance

## Configuration Examples

### Medium Scale Production
```yaml
execution:
  num_gt_instances: 10
  cas9_simulations_per_gt: 5  # 5 recordings per GT

# Results in:
# - 10 GTs × 5 simulations × 4 tiers × 3 solvers = 600 reconstructions
# - Strong statistical power for recording variance analysis
```

### High Statistical Power
```yaml
execution:
  num_gt_instances: 2  
  cas9_simulations_per_gt: 10  # 10 recordings per GT for maximum statistics

# Results in:  
# - 2 GTs × 10 simulations × 3 tiers × 3 solvers = 180 reconstructions
# - Excellent statistical power with robust confidence intervals
```

### Quick Testing
```yaml
execution:
  num_gt_instances: 1
  cas9_simulations_per_gt: 3  # Minimal statistics for testing

# Results in:
# - 1 GT × 3 simulations × 3 tiers × 3 solvers = 27 reconstructions  
# - Basic variance estimates for development
```

## Analysis Capabilities

### Recording Quality Metrics
- **Recording success rate** across simulation runs
- **Site utilization efficiency** variance
- **Mutation pattern consistency** between runs
- **Information content** stability across recordings

### Solver Robustness Analysis
- **Performance consistency** across recording conditions
- **Recording quality sensitivity** for each solver
- **Failure mode analysis** under different recording scenarios
- **Optimal recording condition** identification per solver

### Statistical Comparisons
- **Pairwise solver comparisons** with proper statistical testing
- **Recording quality impact** quantification  
- **Interaction effects** between solver choice and recording conditions
- **Power analysis** for future experimental design

## Implementation Details

### Code Changes
1. **`cas9_simulations_per_gt` property**: New configuration parameter
2. **Job calculation update**: Accounts for multiple Cas9 simulations per GT
3. **Backward compatibility**: Default value of 1 maintains old behavior
4. **Configuration validation**: Ensures reasonable simulation counts

### File Naming Convention
```
GT_instance_0 → cas9_sim_0_tier_1, cas9_sim_0_tier_2, cas9_sim_0_tier_3, cas9_sim_0_tier_4
             → cas9_sim_1_tier_1, cas9_sim_1_tier_2, cas9_sim_1_tier_3, cas9_sim_1_tier_4
             → cas9_sim_2_tier_1, cas9_sim_2_tier_2, cas9_sim_2_tier_3, cas9_sim_2_tier_4
             → ... (up to cas9_simulations_per_gt)
```

## Comparison with Previous Approach

### ❌ Previous (Incorrect): Multiple Reconstructions per Cas9
```
GT → Cas9_Simulation → [Recon_1, Recon_2, Recon_3, Recon_4, Recon_5]
```
- Only captures reconstruction algorithm variance
- Ignores recording process variability
- Overestimates reconstruction precision
- Misses recording-solver interaction effects

### ✅ Current (Correct): Multiple Cas9 Simulations per GT
```
GT → [Cas9_Sim_1, Cas9_Sim_2, Cas9_Sim_3, Cas9_Sim_4, Cas9_Sim_5] → Reconstructions
```
- Captures both recording AND reconstruction variance
- Reflects real experimental conditions
- Enables recording quality impact assessment
- Provides realistic confidence intervals

## Best Practices

### Simulation Count Guidelines
- **Development/Testing**: 1-3 simulations per GT
- **Basic Statistics**: 3-5 simulations per GT  
- **Robust Analysis**: 5-10 simulations per GT
- **Publication Quality**: 10+ simulations per GT

### Statistical Analysis Recommendations
1. **Always report** recording variance alongside reconstruction variance
2. **Use mixed-effects models** to account for GT and recording random effects
3. **Test for recording-solver interactions** before making solver comparisons
4. **Validate normality assumptions** for confidence interval calculations
5. **Account for multiple testing** when comparing many solver pairs

This approach transforms the cascade workflow into a realistic simulation of the complete phylogenetic reconstruction pipeline, from recording to reconstruction, enabling meaningful statistical inference about solver performance under realistic conditions.