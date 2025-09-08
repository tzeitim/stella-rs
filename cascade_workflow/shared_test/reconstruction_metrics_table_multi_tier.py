#!/usr/bin/env python3
"""
Multi-Tier Cas9 Reconstruction Metrics Analysis

This script implements 4 different tiers of Cas9 simulation regimes with varying
levels of recording fidelity, then applies all reconstruction solvers to analyze
how recording quality affects reconstruction performance.

Tiers:
- Tier 1: Extremely high fidelity (near-perfect recording)
- Tier 2: High fidelity (many sites, varied mutation rates)  
- Tier 3: Medium fidelity (fewer sites, moderate parameters)
- Tier 4: Low fidelity (few sites, uniform mutation rates)
"""

from copy import deepcopy
import numpy as np
import pandas as pd
import cassiopeia as cass
from cassiopeia.simulator import BirthDeathFitnessSimulator, UniformLeafSubsampler, Cas9LineageTracingDataSimulator
import stellars
import convexml
from dataclasses import dataclass
from typing import Dict, Any


@dataclass
class Cas9SimulationTier:
    """Configuration for a Cas9 simulation tier."""
    
    name: str
    description: str
    k: int  # number of cassettes (integrations)
    cassette_size: int  # number of sites per cassette
    m: int  # number of unique mutations per site
    mutation_rates: list  # mutation rates for different sites/cassettes
    state_priors_exponent: float
    
    @property
    def total_sites(self) -> int:
        """Total recording sites across all cassettes."""
        return self.k * self.cassette_size
    
    def get_effective_mutation_rate(self) -> float:
        """Get the effective mutation rate across all sites."""
        return np.mean(self.mutation_rates)


# Define the 4 Cas9 simulation tiers
CAS9_TIERS = {
    1: Cas9SimulationTier(
        name="Tier 1 - Ultra High Fidelity",
        description="Near-perfect recording: 100 integrations × 50 sites with time-scaled mutation rates",
        k=100,  # 100 integrations
        cassette_size=50,  # 50 sites per integration
        m=50,  # 50 unique mutations per site
        # Time-scaled mutation rates: early sites mutate frequently, later sites rarely
        mutation_rates=[2.0 * np.exp(-0.1 * i) for i in range(100)],  # Exponential decay
        state_priors_exponent=1e-6  # Very uniform priors
    ),
    
    2: Cas9SimulationTier(
        name="Tier 2 - High Fidelity", 
        description="Good recording: 50 integrations × 20 sites with varied mutation rates",
        k=50,   # 50 integrations
        cassette_size=20,  # 20 sites per integration  
        m=30,   # 30 unique mutations per site
        # Varied mutation rates across integrations
        mutation_rates=[1.5 * (1.0 + 0.5 * np.sin(0.2 * i)) for i in range(50)],
        state_priors_exponent=1e-5
    ),
    
    3: Cas9SimulationTier(
        name="Tier 3 - Medium Fidelity",
        description="Moderate recording: 20 integrations × 10 sites with some variation",
        k=20,   # 20 integrations
        cassette_size=10,  # 10 sites per integration
        m=20,   # 20 unique mutations per site
        # Some variation in mutation rates
        mutation_rates=[1.0 + 0.3 * (i % 3 - 1) for i in range(20)],  # 0.7, 1.0, 1.3 pattern
        state_priors_exponent=1e-4
    ),
    
    4: Cas9SimulationTier(
        name="Tier 4 - Low Fidelity",
        description="Poor recording: 5 integrations × 3 sites with uniform mutation rates",
        k=5,    # 5 integrations
        cassette_size=3,   # 3 sites per integration
        m=10,   # 10 unique mutations per site
        # Uniform mutation rates (worst case)
        mutation_rates=[0.8] * 5,  # All the same
        state_priors_exponent=1e-3
    )
}


def generate_state_priors(m, exp):
    """Generates state priors for a single cassette."""
    site_priors = np.array([np.random.exponential(exp) for _ in range(m)])
    site_priors /= np.sum(site_priors)
    return {i: site_priors[i] for i in range(m)}


def simulate_gt_with_cas9_tier(tier: Cas9SimulationTier, base_tree=None):
    """
    Simulate a ground truth tree with the specified Cas9 recording tier.
    
    Args:
        tier: Cas9SimulationTier configuration
        base_tree: Optional base tree topology to reuse (for consistency across tiers)
        
    Returns:
        Tuple of (gt_tree, base_tree) where base_tree can be reused for other tiers
    """
    
    # Tree topology configuration (consistent across all tiers)
    tree_config = {
        'N': 1e3,  # number of cells in original tree
        'n': 1e2,  # number of cells in subsampled tree
        'fitness': {
            'birth_waiting_distribution': lambda scale: np.random.exponential(1/scale),
            'initial_birth_scale': 2,
            'death_waiting_distribution': lambda: np.inf,
            'mutation_distribution': lambda: 1 if np.random.uniform() < 0.5 else 0,
            'fitness_distribution': lambda: np.random.normal(0.5, 0.25),
            'fitness_base': 1.1
        }
    }
    
    # Generate or reuse base tree topology
    if base_tree is None:
        topology_simulator = BirthDeathFitnessSimulator(**tree_config['fitness'], num_extant=int(tree_config['N']))
        original_topology = topology_simulator.simulate_tree()
        leaf_subsampler = UniformLeafSubsampler(number_of_leaves=int(tree_config['n']))
        base_tree = leaf_subsampler.subsample_leaves(original_topology)
        base_tree.scale_to_unit_length()
    
    # Create a copy for this tier's simulation
    gt_tree = deepcopy(base_tree)
    
    # Create mutation rates array for each cassette
    cassette_mutation_rates = {}
    for i, base_rate in enumerate(tier.mutation_rates):
        # Convert mutation probability to rate
        mutation_rate = -np.log(1.0 - min(0.99, base_rate * 0.5))  # Cap at 99% probability
        cassette_mutation_rates[i] = mutation_rate
    
    # Simulate Cas9 lineage tracing for each cassette separately
    all_priors = {}  # Store priors for all sites
    
    for cassette_idx in range(tier.k):
        # Generate priors for this cassette
        cassette_priors = generate_state_priors(tier.m, tier.state_priors_exponent)
        
        cassette_simulator = Cas9LineageTracingDataSimulator(
            number_of_cassettes=tier.cassette_size,  # Sites in this cassette
            size_of_cassette=1,  # Each site is size 1
            number_of_states=tier.m,
            mutation_rate=cassette_mutation_rates[cassette_idx],
            state_priors=cassette_priors
        )
        
        # Store priors for this cassette
        for site_in_cassette in range(tier.cassette_size):
            site_idx = cassette_idx * tier.cassette_size + site_in_cassette
            all_priors[site_idx] = cassette_priors
            
        if cassette_idx == 0:
            # First cassette - initialize the character matrix
            cassette_simulator.overlay_data(gt_tree)
            # Store the initial character matrix
            base_char_matrix = gt_tree.character_matrix.copy()
        else:
            # Additional cassettes - extend the character matrix
            cassette_tree = deepcopy(base_tree)
            cassette_simulator.overlay_data(cassette_tree)
            
            # Concatenate character matrices
            new_chars = cassette_tree.character_matrix
            # Rename columns to avoid conflicts
            new_chars.columns = [f"cas{cassette_idx}_site{i}" for i in range(len(new_chars.columns))]
            
            # Merge with existing character matrix
            gt_tree.character_matrix = pd.concat([gt_tree.character_matrix, new_chars], axis=1)
    
    # Set tree parameters
    gt_tree.priors = all_priors
    gt_tree.parameters["stochastic_missing_rate"] = 0
    gt_tree.parameters["heritable_missing_rate"] = 0
    
    # Add tier metadata
    gt_tree.cas9_tier = tier
    gt_tree.cas9_tier_name = tier.name
    gt_tree.cas9_total_sites = tier.total_sites
    
    return gt_tree, base_tree


def get_class_of_solver(solver_name):
    """Returns the solver class based on the solver name."""
    solvers_map = {
        "nj": lambda: cass.solver.NeighborJoiningSolver(add_root=True),
        "maxcut": lambda: cass.solver.MaxCutSolver(),
        "maxcut_greedy": lambda: cass.solver.MaxCutGreedySolver(),
        "greedy": lambda: cass.solver.VanillaGreedySolver(),
        "smj": lambda: cass.solver.SharedMutationJoiningSolver(),
        "spectral": lambda: cass.solver.SpectralSolver(),
        "spectral_greedy": lambda: cass.solver.SpectralGreedySolver()
    }
    
    if solver_name not in solvers_map:
        raise ValueError(f"Unknown solver: {solver_name}")
    
    return solvers_map[solver_name]()


def reconstruct_tree_with_convexml(gt_tree, solver_name):
    """Reconstructs a tree using the specified solver with ConvexML optimization."""
    reconstructed_tree = deepcopy(gt_tree)        
    get_class_of_solver(solver_name).solve(reconstructed_tree)
    
    # Post-process with ConvexML
    reconstructed_tree.reconstruct_ancestral_characters()
    
    tree_newick = reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True)
    leaf_sequences = {leaf: reconstructed_tree.get_character_states(leaf) for leaf in reconstructed_tree.leaves}
    
    min_branch_len = 1e-8
    try:
        branch_lengths = []
        for parent in reconstructed_tree.nodes:
            for child in reconstructed_tree.children(parent):
                bl = reconstructed_tree.get_branch_length(parent, child)
                if bl is not None and bl > 0:
                    branch_lengths.append(bl)
        if branch_lengths:
            min_branch_len = max(min(branch_lengths) * 0.001, 1e-10)
    except:
        pass
    
    convexml_result = convexml.convexml(
        tree_newick=tree_newick,
        leaf_sequences=leaf_sequences,
        minimum_branch_length=min_branch_len,
        resolve_multifurcations_before_branch_length_estimation=True
    )
    
    optimized_tree = cass.data.CassiopeiaTree(
        character_matrix=reconstructed_tree.character_matrix,
        tree=convexml_result["tree_newick"],
        missing_state_indicator=getattr(reconstructed_tree, 'missing_state_indicator', -1)
    )
    
    # Copy over priors and parameters
    if hasattr(reconstructed_tree, 'priors'):
        optimized_tree.priors = reconstructed_tree.priors
    if hasattr(reconstructed_tree, 'parameters'):
        try:
            optimized_tree._parameters = reconstructed_tree._parameters.copy()
        except:
            pass
    
    optimized_tree.reconstruct_ancestral_characters()
    optimized_tree.scale_to_unit_length()
    
    return optimized_tree


def get_lam_and_q(tree, known_priors=True):
    """Gets the mutation rate and collision probability for the tree."""
    cm = tree.character_matrix.to_numpy()
    proportion_mutated = np.sum(cm > 0) / np.sum(cm >= 0)
    lam = -np.log(1.0 - proportion_mutated)

    if known_priors and hasattr(tree, 'priors') and tree.priors:
        # Average collision probability across all sites
        q_values = []
        for site_priors in tree.priors.values():
            priors = np.array(list(site_priors.values()))
            q_values.append(np.sum(priors ** 2))
        q = np.mean(q_values)
    else:
        q = 1.0 / tree.character_matrix.stack().nunique()

    return lam, q


def calculate_metrics_improved(gt_tree, reconstructed_tree):
    """Calculate all metrics using the improved stellars interface."""
    try:
        gt_newick = gt_tree.get_newick(record_branch_lengths=True, record_node_names=True)
        recon_newick = reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True)
        
        # Use improved stellars interface
        rf_distance = stellars.rf_distance_simple(gt_newick, recon_newick, normalize=False)
        triplets_distance = stellars.triplets_distance_simple(gt_newick, recon_newick, method='optimized')
        
        # Other metrics
        parsimony_score = cass.tools.tree_metrics.calculate_parsimony(reconstructed_tree)
        
        with np.errstate(divide='ignore'):
            likelihood_score = cass.tools.tree_metrics.calculate_likelihood_continuous(reconstructed_tree)
        
        # cPHS with explicit parameters (safe approach)
        try:
            lam, q = get_lam_and_q(reconstructed_tree)
            internal_character_states = {}
            for internal_node in reconstructed_tree.internal_nodes:
                node_name = str(internal_node)
                internal_states = reconstructed_tree.get_character_states(internal_node)
                try:
                    if hasattr(internal_states, 'tolist'):
                        internal_character_states[node_name] = internal_states.tolist()
                    else:
                        internal_character_states[node_name] = list(internal_states)
                except Exception:
                    internal_character_states = {}
                    break
            
            cphs_result = stellars.phs(
                tree_newick=recon_newick,
                character_matrix=reconstructed_tree.character_matrix.values.astype(int).tolist(),
                internal_states=internal_character_states,
                mutation_rate=lam,
                collision_probability=q,
                leaf_names=list(reconstructed_tree.character_matrix.index),
                use_internal_states=True
            )
            cphs_score = cphs_result.phs_score
        except Exception as e:
            print(f"Warning: cPHS calculation failed: {e}")
            cphs_score = np.nan
        
        return rf_distance, triplets_distance, parsimony_score, cphs_score, likelihood_score
        
    except Exception as e:
        print(f"Warning: Metrics calculation failed: {e}")
        return np.nan, np.nan, np.nan, np.nan, np.nan


def main():
    """Main function to run multi-tier Cas9 reconstruction analysis with multiple instances."""
    print("Multi-Tier Cas9 Reconstruction Metrics Analysis (4 instances per tier)")
    print("="*80)
    
    # Print tier descriptions
    for tier_num, tier in CAS9_TIERS.items():
        print(f"\n{tier.name}:")
        print(f"  {tier.description}")
        print(f"  Sites: {tier.total_sites} ({tier.k}×{tier.cassette_size})")
        print(f"  Avg mutation rate: {tier.get_effective_mutation_rate():.3f}")
    
    print(f"\n{'='*80}")
    print("Generating 4 instances per tier for statistical analysis...")
    
    # Available solvers
    available_solvers = [
        'nj', 'maxcut', 'maxcut_greedy', 'greedy', 'smj', 
        'spectral', 'spectral_greedy'
    ]
    
    results = []
    
    # Generate 4 instances for each tier
    instances_per_tier = 4
    
    for tier_num, tier_config in CAS9_TIERS.items():
        print(f"\n{'='*60}")
        print(f"Processing {tier_config.name}")
        print(f"{'='*60}")
        
        for instance_num in range(1, instances_per_tier + 1):
            print(f"\nGenerating Instance {instance_num}/4 for {tier_config.name}...")
            
            try:
                # Generate a new Ground Truth tree and apply Cas9 recording layer
                # Note: Each instance gets a completely new GT tree + recording layer
                gt_tree, _ = simulate_gt_with_cas9_tier(tier_config, base_tree=None)
                
                print(f"  Created tree with {gt_tree.character_matrix.shape[1]} recording sites")
                
                # Apply all solvers to this instance
                for solver in available_solvers:
                    print(f"    - Applying solver '{solver}'...")
                    
                    try:
                        # Reconstruct tree
                        reconstructed_tree = reconstruct_tree_with_convexml(gt_tree, solver)
                        
                        # Calculate metrics
                        rf_distance, triplets_distance, parsimony_score, cphs_score, likelihood_score = calculate_metrics_improved(gt_tree, reconstructed_tree)
                        
                        # Store results
                        results.append({
                            'Cas9_tier': tier_num,
                            'Cas9_tier_name': tier_config.name,
                            'instance': instance_num,
                            'recording_sites': tier_config.total_sites,
                            'integrations': tier_config.k,
                            'sites_per_integration': tier_config.cassette_size,
                            'avg_mutation_rate': tier_config.get_effective_mutation_rate(),
                            'solver': solver,
                            'triplets_distance': triplets_distance,
                            'RF_distance': rf_distance,
                            'parsimony_score': parsimony_score,
                            'cPHS': cphs_score,
                            'likelihood_score': likelihood_score
                        })
                        
                    except Exception as e:
                        print(f"      ERROR with {solver}: {e}")
                        # Store failed result
                        results.append({
                            'Cas9_tier': tier_num,
                            'Cas9_tier_name': tier_config.name,
                            'instance': instance_num,
                            'recording_sites': tier_config.total_sites,
                            'integrations': tier_config.k,
                            'sites_per_integration': tier_config.cassette_size,
                            'avg_mutation_rate': tier_config.get_effective_mutation_rate(),
                            'solver': solver,
                            'triplets_distance': np.nan,
                            'RF_distance': np.nan,
                            'parsimony_score': np.nan,
                            'cPHS': np.nan,
                            'likelihood_score': np.nan
                        })
                        
            except Exception as e:
                print(f"  ERROR generating instance {instance_num}: {e}")
                # Create failed entries for all solvers for this instance
                for solver in available_solvers:
                    results.append({
                        'Cas9_tier': tier_num,
                        'Cas9_tier_name': tier_config.name,
                        'instance': instance_num,
                        'recording_sites': tier_config.total_sites,
                        'integrations': tier_config.k,
                        'sites_per_integration': tier_config.cassette_size,
                        'avg_mutation_rate': tier_config.get_effective_mutation_rate(),
                        'solver': solver,
                        'triplets_distance': np.nan,
                        'RF_distance': np.nan,
                        'parsimony_score': np.nan,
                        'cPHS': np.nan,
                        'likelihood_score': np.nan
                    })
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Display results
    print("\n" + "="*100)
    print("MULTI-TIER CAS9 RECONSTRUCTION METRICS TABLE")
    print("="*100)
    print(results_df.to_string(index=False, float_format='%.6f'))
    
    # Statistical analysis by tier and solver
    print(f"\n" + "="*100)
    print("STATISTICAL ANALYSIS BY CAS9 TIER AND SOLVER")
    print("="*100)
    
    numeric_cols = ['triplets_distance', 'RF_distance', 'parsimony_score', 'cPHS', 'likelihood_score']
    
    for tier_num in sorted(results_df['Cas9_tier'].unique()):
        tier_data = results_df[results_df['Cas9_tier'] == tier_num]
        tier_name = tier_data.iloc[0]['Cas9_tier_name']
        sites = tier_data.iloc[0]['recording_sites']
        
        print(f"\n{tier_name} ({sites} sites) - Statistics across 4 instances:")
        print("-" * 80)
        
        # Overall tier statistics (across all solvers and instances)
        print("Overall Tier Performance:")
        for col in numeric_cols:
            if col in tier_data.columns:
                mean_val = tier_data[col].mean()
                std_val = tier_data[col].std()
                min_val = tier_data[col].min()
                max_val = tier_data[col].max()
                print(f"  {col:<20}: mean={mean_val:8.4f} ± {std_val:6.4f} [range: {min_val:8.4f} - {max_val:8.4f}]")
        
        print("\nBy Solver (mean ± std across 4 instances):")
        for solver in sorted(tier_data['solver'].unique()):
            solver_data = tier_data[tier_data['solver'] == solver]
            print(f"  {solver:>15}:", end="")
            for col in ['triplets_distance', 'RF_distance', 'cPHS']:
                if col in solver_data.columns:
                    mean_val = solver_data[col].mean()
                    std_val = solver_data[col].std()
                    print(f" {col.split('_')[0][:8]}={mean_val:6.4f}±{std_val:5.4f}", end="")
            print()
    
    # Best performing solvers by tier
    print(f"\n" + "="*100)
    print("BEST PERFORMING SOLVERS BY TIER (lowest triplets distance)")
    print("="*100)
    
    for tier_num in sorted(results_df['Cas9_tier'].unique()):
        tier_data = results_df[results_df['Cas9_tier'] == tier_num]
        tier_name = tier_data.iloc[0]['Cas9_tier_name']
        
        # Calculate mean performance by solver
        solver_performance = tier_data.groupby('solver')['triplets_distance'].agg(['mean', 'std']).round(6)
        solver_performance = solver_performance.sort_values('mean')
        
        print(f"\n{tier_name}:")
        print("  Rank  Solver              Mean Triplets Distance    Std")
        print("  " + "-" * 55)
        for rank, (solver, row) in enumerate(solver_performance.iterrows(), 1):
            print(f"  {rank:2d}.   {solver:<15}   {row['mean']:12.6f}        {row['std']:8.6f}")
    
    # Variability analysis
    print(f"\n" + "="*100)
    print("VARIABILITY ANALYSIS (Standard Deviations)")
    print("="*100)
    print("Lower standard deviation = more consistent performance across instances")
    
    tier_variability = []
    for tier_num in sorted(results_df['Cas9_tier'].unique()):
        tier_data = results_df[results_df['Cas9_tier'] == tier_num]
        tier_name = tier_data.iloc[0]['Cas9_tier_name']
        
        avg_std = tier_data.groupby('solver')['triplets_distance'].std().mean()
        tier_variability.append((tier_name, avg_std))
        
        print(f"\n{tier_name}: Average std across solvers = {avg_std:.6f}")
        solver_std = tier_data.groupby('solver')['triplets_distance'].std().sort_values()
        for solver, std_val in solver_std.items():
            consistency = "Very consistent" if std_val < 0.01 else "Consistent" if std_val < 0.05 else "Variable"
            print(f"  {solver:<15}: {std_val:8.6f} ({consistency})")
    
    # Save results
    output_file = "multi_tier_cas9_reconstruction_metrics.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\n" + "="*100)
    print(f"Results saved to: {output_file}")
    print(f"Total reconstructions: {len(results_df)} (4 tiers × 4 instances × 7 solvers)")
    
    return results_df


if __name__ == "__main__":
    # Initialize stellars logging
    stellars.init_logging("info")
    
    # Run the multi-tier analysis
    results_df = main()