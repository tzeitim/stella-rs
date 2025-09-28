# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python [conda env:cas]
#     language: python
#     name: conda-env-cas-py
# ---

# %% [markdown]
# June 2025
# Author: Pini Zilber
# Example Notebook for Calculating the cPHS Test Statistic
# Based on joint work with Sebastian Prillo, Nir Yosef, and Boaz Nadler

# %%
from copy import deepcopy
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import binom

import cassiopeia as cass
from cassiopeia.simulator import BirthDeathFitnessSimulator, UniformLeafSubsampler, Cas9LineageTracingDataSimulator

import stellars
import convexml

# Initialize stellars logging at info level to see what's happening
#stellars.init_logging("info")
#stellars.init_logging("debug")
from convexml import to_newick
import pickle
import glob
from pathlib import Path
import argparse
# %%
### Command line argument parsing ###

def parse_arguments():
    parser = argparse.ArgumentParser(description='Phylogenetic reconstruction simulation and analysis')
    parser.add_argument('--mode', choices=['simulate', 'load_pipeline', 'load_results'], default='simulate',
                       help='Mode: simulate (run full simulation), load_pipeline (load objects and re-compute), or load_results (load existing pipeline results)')
    parser.add_argument('--pipeline_dir', default='output/sim',
                       help='Directory containing pipeline results (for load_pipeline/load_results modes)')
    parser.add_argument('--repetitions', type=int, default=3,
                       help='Number of repetitions for simulation mode')
    parser.add_argument('--solvers', nargs='+',
                       default=['nj', 'maxcut', 'maxcut_greedy', 'greedy', 'smj', 'spectral', 'spectral_greedy'],
                       help='List of solvers to use')
    parser.add_argument('--output_prefix', default='simulation_phs_control_results',
                       help='Prefix for output files')
    return parser.parse_args()

# Parse command line arguments
args = parse_arguments()
MODE = args.mode
PIPELINE_DIR = args.pipeline_dir
repetitions = args.repetitions

### configurations for the simulation ###

tree_config = {
    # tree topology configuration
    'N': 1e3, # number of cells in the original tree (before subsampling)
    'n': 1e2, # number of cells in the subsampled tree
    'fitness': { # see cass.simulator.BirthDeathFitnessSimulator
        'birth_waiting_distribution': lambda scale: np.random.exponential(1/scale),
        'initial_birth_scale': 2,
        'death_waiting_distribution': lambda: np.inf, # cells don't die
        'mutation_distribution': lambda: 1 if np.random.uniform() < 0.5 else 0,
        'fitness_distribution': lambda: np.random.normal(0.5, 0.25),
        'fitness_base': 1.1
    },

    # mutations model configuration
    'k': 50, # sequence length
    'm': 50, # number of unique mutations
    'rho': 0.5, # character mutation probability at a leaf
    'state_priors_exponents': 1e-5, # exponent of the random state priors generator
}

# list of solvers to be used in the simulation
# the solvers are defined in cassiopeia/solver/solvers.py
solvers = args.solvers


# %%
### functions for ground-truth tree simulation and tree reconstructions ###

def set_tree_root_to_zeros(tree):
    """
    Sets all character states at the root of the tree to 0.
    """
    tree.set_character_states(tree.root, [0] * tree.n_character)


def get_class_of_solver(solver_name):
    """
    Returns the solver class based on the solver name.
    """
    if solver_name == "nj":
        return cass.solver.NeighborJoiningSolver(add_root=True)
    elif solver_name == "maxcut":
        return cass.solver.MaxCutSolver()
    elif solver_name == "maxcut_greedy":
        return cass.solver.MaxCutGreedySolver()
    elif solver_name == "greedy":
        return cass.solver.VanillaGreedySolver()
    elif solver_name == "smj":
        return cass.solver.SharedMutationJoiningSolver()
    elif solver_name == "spectral":
        return cass.solver.SpectralSolver()
    elif solver_name == "spectral_greedy":
        return cass.solver.SpectralGreedySolver()
    else:
        raise ValueError(f"Unknown solver: {solver_name}")


def generate_state_priors(k, m, exp):
    """
    Generates state priors (q_i) for character states
    """

    # generate random state priors
    state_priors = np.array([
        np.random.exponential(exp)
        for _ in range(m)
    ])
    state_priors /= np.sum(state_priors)
    
    # return a dictionary
    return {i: state_priors[i] for i in range(m)}


def simulate_gt(tree_config):
    """
    Simulates the ground truth tree
    """

    # simulate the original tree topology
    topology_simulator = BirthDeathFitnessSimulator(**tree_config['fitness'], num_extant=int(tree_config['N']))
    original_topology = topology_simulator.simulate_tree()

    # subsample the tree topology to get the GT tree
    leaf_subsampler = UniformLeafSubsampler(number_of_leaves=int(tree_config['n']))
    gt_tree = leaf_subsampler.subsample_leaves(original_topology)
    gt_tree.scale_to_unit_length()

    # generate lineage tracing data sequences for the GT tree
    lam = -np.log(1.0 - tree_config["rho"])
    k, m, exp = tree_config['k'], tree_config['m'], tree_config['state_priors_exponents']
    priors = generate_state_priors(k, m, exp)
    lt_simulator = Cas9LineageTracingDataSimulator(
        number_of_cassettes=k,
        size_of_cassette=1,
        number_of_states=m,
        mutation_rate=lam,
        state_priors=priors
    )
    lt_simulator.overlay_data(gt_tree)
    #set_tree_root_to_zeros(gt_tree)

    # update some parameters of the tree
    # set the same priros for all character locations
    gt_tree.priors = {i: priors for i in range(k)}
    gt_tree.parameters["stochastic_missing_rate"] = 0 # change this if introducing missing data
    gt_tree.parameters["heritable_missing_rate"] = 0 # change this if introducing missing data

    return gt_tree


def reconstruct(gt_tree, solvers):
    """
    Reconstructs the tree from its leaves using different solvers.
    Now uses ConvexML for branch length estimation instead of Cassiopeia's built-in methods.
    """

    reconstructions = []
    for solver in solvers:
        # solve the tree reconstruction problem
        reconstructed_tree = deepcopy(gt_tree)        
        get_class_of_solver(solver).solve(reconstructed_tree)
        
        # post-process the reconstructed tree
        print("next solver...")
        # Note: removed resolve_multifurcations() as ConvexML handles this
        reconstructed_tree.reconstruct_ancestral_characters()
        #set_tree_root_to_zeros(reconstructed_tree)
        
        # Use ConvexML for branch length estimation (replaces deprecated Cassiopeia methods)
        # Get initial tree newick and leaf sequences for ConvexML
        tree_newick = reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True)
        leaf_sequences = {leaf: reconstructed_tree.get_character_states(leaf) for leaf in reconstructed_tree.leaves}
        
        # Estimate minimum branch length from existing tree if possible
        min_branch_len = 1e-8  # Start with smaller default
        try:
            branch_lengths = []
            for parent in reconstructed_tree.nodes:
                for child in reconstructed_tree.children(parent):
                    bl = reconstructed_tree.get_branch_length(parent, child)
                    if bl is not None and bl > 0:
                        branch_lengths.append(bl)
            
            if branch_lengths:
                # Use small fraction to avoid ConvexML "too large" errors
                min_branch_len = min(branch_lengths) * 0.001  # Use 0.1% of minimum
                # Ensure it's not too small either
                min_branch_len = max(min_branch_len, 1e-10)
        except:
            pass  # Use default
        
        # Call ConvexML to optimize branch lengths (replaces IIDExponentialMLE)
        convexml_result = convexml.convexml(
            tree_newick=tree_newick,
            leaf_sequences=leaf_sequences,
            minimum_branch_length=min_branch_len,
            resolve_multifurcations_before_branch_length_estimation=True
        )
        
        # Create new tree object with optimized branch lengths from ConvexML
        optimized_tree = cass.data.CassiopeiaTree(
            character_matrix=reconstructed_tree.character_matrix,
            tree=convexml_result["tree_newick"],
            missing_state_indicator=getattr(reconstructed_tree, 'missing_state_indicator', -1)
        )
        
        # Copy over priors and parameters from the original tree
        if hasattr(reconstructed_tree, 'priors'):
            optimized_tree.priors = reconstructed_tree.priors
        
        # Copy parameters by accessing the underlying dictionary directly
        if hasattr(reconstructed_tree, 'parameters'):
            try:
                # Try to access the internal parameters dict if possible
                optimized_tree._parameters = reconstructed_tree._parameters.copy()
            except:
                # If that doesn't work, just skip parameters copying
                # The tree should still work for basic operations
                pass
        
        # Reconstruct ancestral characters for the new tree with optimized branch lengths
        optimized_tree.reconstruct_ancestral_characters()
        
        # Scale to unit length after ConvexML optimization
        optimized_tree.scale_to_unit_length()
        
        # Store the optimized tree
        reconstructions.append(optimized_tree)
    
    return reconstructions



# %%
### functions for cPHS calculation ###

def calc_phs_and_lca_heights(tree):
    """
    Calculates the pairwise homoplasy scores (PHS) and the heights of the least common ancestors (LCA) for all pairs of leaves in the tree
    """
    cm = tree.character_matrix
    n, k = cm.shape
    res = pd.DataFrame(np.zeros(shape=(n,n)) - 1, index=cm.index, columns=cm.index)
    lca_heights = pd.DataFrame(np.zeros(shape=(n,n)) - 1, index=cm.index, columns=cm.index)
    
    for lca in tree.internal_nodes: # including the root
        lca_states = tree.get_character_states(lca)
        lca_height = tree.get_time(lca)
        children = tree.children(lca)
        for i in range(len(children)):
            ch_1 = children[i]
            for j in range(i + 1, len(children)):
                ch_2 = children[j]
                for leaf_1 in tree.leaves_in_subtree(ch_1):
                    leaf_1_states = tree.get_character_states(leaf_1)
                    for leaf_2 in tree.leaves_in_subtree(ch_2):
                        leaf_2_states = tree.get_character_states(leaf_2)
                        phs = 0
                        for character_id in range(k):
                            if (
                                lca_states[character_id] == 0
                                and leaf_1_states[character_id] > 0
                                and leaf_1_states[character_id]
                                == leaf_2_states[character_id]
                            ):
                                phs += 1
                        res.loc[leaf_1, leaf_2] = phs
                        res.loc[leaf_2, leaf_1] = phs
                        lca_heights.loc[leaf_1, leaf_2] = lca_height
                        lca_heights.loc[leaf_2, leaf_1] = lca_height
    for leaf in tree.leaves:
        res.loc[leaf, leaf] = 0
        lca_heights.loc[leaf, leaf] = 0

    # Make sure all pairs of leaves got iterated over
    if not ((res >= 0).all().all()):
        raise Exception(
            "All pairs of leaves should have been iterated over, "
            "but some were missed!"
        )

    return res, lca_heights


def calc_cPHS(tree, lam, q):
    """
    Calculates the cPHS test statistic for the tree using both Python implementation and stellars
    """
    phs, lca_heights = calc_phs_and_lca_heights(tree)

    ## get parameters values
    cm = tree.character_matrix
    n, k = cm.shape
    m = cm.stack().nunique()
    # estimate mutation rate
    cm_array = tree.character_matrix.to_numpy()
    proportion_mutated = np.sum(cm_array > 0) / np.sum(cm_array >= 0)
    lam_estimate = -np.log(1.0 - proportion_mutated)

    ## calculate p-values using Python implementation
    alpha = np.exp(-lam * lca_heights)
    beta =  (1 - np.exp(-lam * (1 - lca_heights)))
    prob = alpha * beta**2 * q
    prob[lca_heights.to_numpy() == 1] = 1
        # fixing a pathology due to internal nodes with height 1 (with mutationless edges to their children leaves)
    cdf = binom.cdf(phs - 1, k, prob)
    pvalues_matrix = 1 - cdf
    off_diagonal_entries = ~np.eye(pvalues_matrix.shape[0],dtype=bool)
    pvalues = np.array(pvalues_matrix[off_diagonal_entries])
    pvalues[pvalues == 0] = np.finfo(float).eps # zeros cannot be real zeros
    pvalues[pvalues < 1] # ignore p-values of pairs with no homoplasies

    ## adjust p-values for multiple testing
    N = len(pvalues)
    pvalues_sorted = np.sort(pvalues)
    adjusted_pvalues = pvalues_sorted * N / np.arange(1, N+1)
    python_cPHS = min(adjusted_pvalues)

    ## calculate cPHS using stellars optimized implementation
    try:
        # Convert tree to newick format using CassiopeiaTree's method
        tree_newick = tree.get_newick(record_branch_lengths=True, record_node_names=True)
        #print(tree_newick)
        # Convert character matrix to list format expected by stellars
        character_matrix = tree.character_matrix.values.astype(int).tolist()
        
        # CRITICAL FIX: Get leaf names in the same order as character matrix rows
        # This fixes the character state mapping issue that was causing ultra-low p-values
        leaf_names = list(tree.character_matrix.index)  # DataFrame index gives us the leaf names in correct order
        print(f"DEBUG: Using {len(leaf_names)} leaf names for proper character matrix mapping")
        
        # Extract internal character states from Cassiopeia for stellars
        # Create a simple mapping from complex node names to simple IDs
        internal_character_states = {}
        node_name_map = {}
        
        for internal_node in tree.internal_nodes:
            # Use the original node name (no simplification needed now!)
            node_name = str(internal_node)
            
            internal_states = tree.get_character_states(internal_node)
            try:
                if hasattr(internal_states, 'tolist'):
                    internal_character_states[node_name] = internal_states.tolist()
                else:
                    internal_character_states[node_name] = list(internal_states)
            except Exception:
                # Fallback to letting stellars infer states
                internal_character_states = {}
                break
        
        # Get the tree WITH node names preserved in Newick
        tree_newick = tree.get_newick(record_branch_lengths=True, record_node_names=True)
        
        print(f"DEBUG: Extracted {len(internal_character_states)} internal node states")
        if len(internal_character_states) > 0:
            print(f"DEBUG: Example internal nodes: {list(internal_character_states.keys())[:3]}")
            print(f"DEBUG: Tree newick (first 200 chars): {tree_newick[:200]}...")
        
        stellars_result = stellars.phs_optimized(
            tree_newick=tree_newick,
            character_matrix=character_matrix,
            internal_character_states=internal_character_states,  # Use Cassiopeia's internal states
            mutation_rate=lam,  # Use the mutation rate from function parameters
            collision_probability=q,   # Use collision probability from function parameters
            missing_state=-1,
            unedited_state=0,
            use_provided_internal_states=True,  # Use the exact same ancestral states as Cassiopeia
            leaf_names=leaf_names  # CRITICAL FIX: Pass leaf names for correct character matrix mapping
        )
        stellars_cPHS = stellars_result['phs_score']
    except Exception as e:
        print(f"Warning: stellars PHS calculation failed: {e}")
        stellars_cPHS = np.nan

    ## return both cPHS scores
    return (python_cPHS, stellars_cPHS)
    


# %%
### functions for metrics calculation ###

def get_lam_and_q(tree, known_priors=True):
    """
    Gets the mutation rate (lam) and collision probability (q) for the tree
    """
    # estimate the mutation rate (lam) from the tree leaves
    cm = tree.character_matrix.to_numpy()
    proportion_mutated = np.sum(cm > 0) / np.sum(cm >= 0)
    lam = -np.log(1.0 - proportion_mutated)

    # get the collision probability (q) for cPHS.
    # if prior are not known, set q=1/m
    if known_priors:
        priors = np.array(list(tree.priors[0].values())) # get priors of first character, they're all the same
        q = np.sum(priors ** 2)
    else:
        q = 1/tree.character_matrix.stack().nunique()

    return lam, q

def load_pipeline_objects(pipeline_dir):
    """
    Load GT trees and CAS9 instances from pipeline output directory
    """
    pipeline_path = Path(pipeline_dir)

    # Load GT trees
    gt_tree_files = sorted(list(pipeline_path.glob("gt_trees/gt_tree_instance_*.pkl")))
    gt_trees = []
    for gt_file in gt_tree_files:
        with open(gt_file, 'rb') as f:
            gt_tree = pickle.load(f)
        gt_trees.append(gt_tree)

    # Load CAS9 instances (these are the trees with recording data)
    cas9_files = sorted(list(pipeline_path.glob("cas9_instances/instance*_sim*_tier*_instance.pkl")))
    cas9_instances = []
    for cas9_file in cas9_files:
        with open(cas9_file, 'rb') as f:
            cas9_instance = pickle.load(f)
        cas9_instances.append(cas9_instance)

    print(f"Loaded {len(gt_trees)} GT trees and {len(cas9_instances)} CAS9 instances from {pipeline_dir}")

    return gt_trees, cas9_instances

def reconstruct_from_cas9_instance(cas9_tree, solvers):
    """
    Reconstruct trees from a CAS9 instance (tree with recording data) using different solvers
    """
    reconstructions = []
    for solver in solvers:
        # solve the tree reconstruction problem
        reconstructed_tree = deepcopy(cas9_tree)
        get_class_of_solver(solver).solve(reconstructed_tree)

        # post-process the reconstructed tree
        print(f"next solver: {solver}...")
        # Note: removed resolve_multifurcations() as ConvexML handles this
        reconstructed_tree.reconstruct_ancestral_characters()
        #set_tree_root_to_zeros(reconstructed_tree)

        # Use ConvexML for branch length estimation (replaces deprecated Cassiopeia methods)
        # Get initial tree newick and leaf sequences for ConvexML
        tree_newick = reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True)
        leaf_sequences = {leaf: reconstructed_tree.get_character_states(leaf) for leaf in reconstructed_tree.leaves}

        # Estimate minimum branch length from existing tree if possible
        min_branch_len = 1e-8  # Start with smaller default
        try:
            branch_lengths = []
            for parent in reconstructed_tree.nodes:
                for child in reconstructed_tree.children(parent):
                    bl = reconstructed_tree.get_branch_length(parent, child)
                    if bl is not None and bl > 0:
                        branch_lengths.append(bl)

            if branch_lengths:
                # Use small fraction to avoid ConvexML "too large" errors
                min_branch_len = min(branch_lengths) * 0.001  # Use 0.1% of minimum
                # Ensure it's not too small either
                min_branch_len = max(min_branch_len, 1e-10)
        except:
            pass  # Use default

        # Call ConvexML to optimize branch lengths (replaces IIDExponentialMLE)
        convexml_result = convexml.convexml(
            tree_newick=tree_newick,
            leaf_sequences=leaf_sequences,
            minimum_branch_length=min_branch_len,
            resolve_multifurcations_before_branch_length_estimation=True
        )

        # Create new tree object with optimized branch lengths from ConvexML
        optimized_tree = cass.data.CassiopeiaTree(
            character_matrix=reconstructed_tree.character_matrix,
            tree=convexml_result["tree_newick"],
            missing_state_indicator=getattr(reconstructed_tree, 'missing_state_indicator', -1)
        )

        # Copy over priors and parameters from the original tree
        if hasattr(reconstructed_tree, 'priors'):
            optimized_tree.priors = reconstructed_tree.priors

        # Copy parameters by accessing the underlying dictionary directly
        if hasattr(reconstructed_tree, 'parameters'):
            try:
                # Try to access the internal parameters dict if possible
                optimized_tree._parameters = reconstructed_tree._parameters.copy()
            except:
                # If that doesn't work, just skip parameters copying
                # The tree should still work for basic operations
                pass

        # Reconstruct ancestral characters for the new tree with optimized branch lengths
        optimized_tree.reconstruct_ancestral_characters()

        # Scale to unit length after ConvexML optimization
        optimized_tree.scale_to_unit_length()

        # Store the optimized tree
        reconstructions.append(optimized_tree)

    return reconstructions

def load_pipeline_results(pipeline_dir):
    """
    Load existing reconstruction results from pipeline parquet files
    """
    pipeline_path = Path(pipeline_dir)

    # Find all parquet files in partitioned_results
    parquet_files = list(pipeline_path.glob("partitioned_results/**/*.parquet"))

    if not parquet_files:
        raise FileNotFoundError(f"No parquet files found in {pipeline_path}/partitioned_results/")

    print(f"Found {len(parquet_files)} parquet files in {pipeline_dir}")

    # Load and combine all parquet files
    all_results = []
    for parquet_file in parquet_files:
        df = pd.read_parquet(parquet_file)
        all_results.append(df)

    combined_df = pd.concat(all_results, ignore_index=True)

    print(f"Loaded {len(combined_df)} total reconstructions")
    print(f"Available solvers: {sorted(combined_df['solver'].unique())}")
    print(f"GT instances: {sorted(combined_df['gt_instance_id'].unique())}")

    return combined_df

def convert_pipeline_results_to_metrics_format(pipeline_df):
    """
    Convert pipeline results DataFrame to the format used by simulation_phs_control.py
    """

    # Group by gt_instance_id to process each "repetition"
    metrics = []

    for gt_instance in sorted(pipeline_df['gt_instance_id'].unique()):
        instance_data = pipeline_df[pipeline_df['gt_instance_id'] == gt_instance]

        # Create dummy gt_metrics (we'll use the first row's GT values)
        first_row = instance_data.iloc[0]
        gt_metrics = {
            'parsimony': first_row.get('parsimony_score', 0),  # We don't have GT parsimony
            'likelihood': first_row.get('log_likelihood_gt', 0),
            'cPHS': first_row.get('cPHS_gt', 0),
            'cPHS-stellars': first_row.get('cPHS_gt', 0)
        }

        # Create reconstructions_metrics for each solver
        reconstructions_metrics = []

        for _, row in instance_data.iterrows():
            recon_metrics = {
                'parsimony': row.get('parsimony_score', row.get('parsimony_total_mutations', 0)),
                'likelihood': row.get('log_likelihood', 0),
                'cPHS': row.get('cPHS', 0),
                'cPHS-stellars': row.get('cPHS', 0),  # Pipeline uses single cPHS value
                'rf': row.get('RF_distance', 0),
                'triplets': 1 - row.get('triplets_distance', 0)  # Convert distance to correct fraction
            }
            reconstructions_metrics.append(recon_metrics)

        metrics.append({
            'gt_metrics': gt_metrics,
            'reconstructions_metrics': reconstructions_metrics,
            'instance_data': instance_data  # Keep original data for reference
        })

    return metrics

def standardize_columns(df, mode):
    """
    Standardize column names and order for comparison between modes
    """

    # Define the standard column mapping
    if mode == 'load_results':
        # Map pipeline columns to standard names
        column_mapping = {
            'RF_distance': 'rf_distance',
            'cPHS': 'cPHS_statistic',
            'triplets_distance': 'triplets_distance_raw',
            'parsimony_score': 'parsimony_statistic',
            'parsimony_total_mutations': 'parsimony_statistic',  # fallback
            'log_likelihood': 'log_likelihood_statistic',
            'gt_instance_id': 'gt_instance',
            'recording_sites': 'config_k',
            'states_per_site': 'config_m',
            'gt_tree_size': 'config_N',
            'sampled_tree_size': 'config_n',
            'proportion_mutated_simulation': 'config_rho'
        }

        # Apply mapping to existing columns
        df = df.rename(columns={k: v for k, v in column_mapping.items() if k in df.columns})

        # Add missing standard columns that simulation mode creates
        if 'repetition' not in df.columns:
            # Map gt_instance_id to repetition numbers
            df['repetition'] = df['gt_instance'] + 1

        # Standardize metric names to match simulation output
        if 'cPHS_statistic' in df.columns:
            df['cPHS statistic (python)'] = df['cPHS_statistic']
            df['cPHS statistic (stellars)'] = df['cPHS_statistic']

        if 'rf_distance' in df.columns:
            df['rf distance'] = df['rf_distance']

        if 'triplets_distance_raw' in df.columns:
            df['triplets distance'] = df['triplets_distance_raw']

        if 'parsimony_statistic' in df.columns:
            df['parsimony statistic'] = df['parsimony_statistic']

        if 'log_likelihood_statistic' in df.columns:
            df['log-likelihood statistic'] = -df['log_likelihood_statistic']  # Convert to positive as simulation does

        # Add distance columns (set to 0 since we don't have GT comparison in pipeline mode)
        if 'parsimony distance' not in df.columns:
            df['parsimony distance'] = 0
        if 'likelihood distance' not in df.columns:
            df['likelihood distance'] = 0

    # Define standard column order following pipeline structure
    standard_columns = [
        # Core pipeline order (exact original names when available)
        'reconstruction_id',
        'gt_instance_id',
        'cas9_simulation_id',
        'reconstruction_num',
        'cas9_tier',
        'cas9_tier_name',
        'recording_sites',
        'states_per_site',
        'solver',
        'computation_time_seconds',
        'run_name',
        'experiment_id',
        'gt_tree_size',
        'sampled_tree_size',
        'lam_simulation',
        'q_simulation',
        'proportion_mutated_simulation',
        'q_simulation_source',
        'lam_gt',
        'q_gt',
        'proportion_mutated_gt',
        'phs_lam_source',
        'phs_q_source',
        'triplets_distance',
        'RF_distance',
        'parsimony_total_mutations',
        'parsimony_computation_time_ms',
        'parsimony_method_used',
        'parsimony_internal_states_inferred',
        'cPHS',
        'cPHS_gt',
        'likelihood_computation_time_ms',
        'reconstructed_tree_path',
        'worker_cas9_instance_path',
        'worker_solver',
        'worker_tier',
        'worker_tier_name',
        'parsimony_score',
        'log_likelihood',
        'log_likelihood_simulation',
        'log_likelihood_gt',

        # Additional compatibility fields for simulation mode
        'repetition',
        'gt_instance',
        'parsimony statistic',
        'log-likelihood statistic',
        'cPHS statistic (python)',
        'cPHS statistic (stellars)',
        'rf distance',
        'triplets distance',
        'parsimony distance',
        'likelihood distance',
        'config_N',
        'config_n',
        'config_k',
        'config_m',
        'config_rho',
        'config_state_priors_exp'
    ]

    # Select only columns that exist and in standard order
    available_columns = [col for col in standard_columns if col in df.columns]
    df_ordered = df[available_columns].copy()

    return df_ordered

def add_relative_parsimony(df, mode):
    """
    Add relative parsimony scores (solver parsimony / minimum solver parsimony for each tree instance).
    """
    if df.empty:
        print("Warning: Cannot compute relative parsimony - empty DataFrame")
        df['relative_parsimony'] = np.nan
        return df

    # Determine parsimony column name based on mode
    if mode == 'load_results':
        parsimony_col = 'parsimony_score'  # Pipeline uses this column name
    else:
        parsimony_col = 'parsimony statistic'  # Simulation uses this column name

    if parsimony_col not in df.columns:
        print(f"Warning: Cannot compute relative parsimony - missing {parsimony_col} column")
        df['relative_parsimony'] = np.nan
        return df

    print(f"Computing relative parsimony using {parsimony_col} column")

    if mode == 'load_results':
        # For pipeline results, group by tree instance using original identifiers
        # Extract tree instance info from reconstruction_id
        df['tree_instance'] = df['reconstruction_id'].str.extract(r'(instance\d+_sim\d+_.*?_tier\d+)')

        # Handle missing tree_instance by using gt_instance_id + cas9_tier as fallback
        missing_tree_instance = df['tree_instance'].isna()
        if missing_tree_instance.any():
            print(f"Warning: {missing_tree_instance.sum()} rows missing tree_instance, using fallback grouping")
            df.loc[missing_tree_instance, 'tree_instance'] = (
                df.loc[missing_tree_instance, 'gt_instance'].astype(str) + '_' +
                df.loc[missing_tree_instance, 'cas9_tier'].astype(str)
            )

        grouping_col = 'tree_instance'
    else:
        # For simulation modes, group by repetition (each repetition = unique tree instance)
        grouping_col = 'repetition'

    # Compute minimum parsimony for each tree instance
    min_parsimony = df.groupby(grouping_col)[parsimony_col].transform('min')

    # Compute relative parsimony
    df['relative_parsimony'] = df[parsimony_col] / min_parsimony

    # Log statistics
    print("Relative parsimony statistics by solver:")
    relative_stats = df.groupby('solver')['relative_parsimony'].agg(['mean', 'std', 'min', 'max', 'count'])
    print(relative_stats.round(4))

    # Clean up temporary column
    if 'tree_instance' in df.columns and mode == 'load_results':
        df = df.drop('tree_instance', axis=1)

    return df

def calculate_metrics(gt_tree, reconstructions):
    """
    Calculates the metrics for the reconstructed trees
    """

    # get the mutation rate (lam) and collision probability (q) for the cPHS test
    lam, q = get_lam_and_q(gt_tree)

    # calculate metrics for the ground truth tree
    gt_parsimony = cass.tools.tree_metrics.calculate_parsimony(gt_tree)
    with np.errstate(divide='ignore'):
        gt_likelihood = cass.tools.tree_metrics.calculate_likelihood_continuous(gt_tree)
    print("gt tree")
    gt_cPHS_python, gt_cPHS_stellars = calc_cPHS(gt_tree, lam, q)


    gt_metrics = {'parsimony': gt_parsimony, 'likelihood': gt_likelihood, 'cPHS': gt_cPHS_python, 'cPHS-stellars': gt_cPHS_stellars}

    # calculate metrics for the reconstructed trees
    reconstructions_metrics = []
    for tree in reconstructions:
        # parsimony
        parsimony = cass.tools.tree_metrics.calculate_parsimony(tree)
        # likelihood
        with np.errstate(divide='ignore'):
            likelihood = cass.tools.tree_metrics.calculate_likelihood_continuous(tree)
        # cPHS
        print("tree")
        cPHS_python, cPHS_stellars = calc_cPHS(tree, lam, q)
        # RF
        rf, rf_max = cass.critique.compare.robinson_foulds(gt_tree, tree)
        rf = rf / rf_max
        #rf = rf / (2*(gt_tree.n_cell-3))
        # triplets
        triplets = cass.critique.compare.triplets_correct(gt_tree, tree)
        triplets = np.mean(list(triplets[0].values()))
        # store the metrics
        reconstructions_metrics.append({
            'parsimony': parsimony,
            'likelihood': likelihood,
            'cPHS': cPHS_python,
            'cPHS-stellars': cPHS_stellars,
            'rf': rf,
            'triplets': triplets
        })
    
    return gt_metrics, reconstructions_metrics


# %%
### main simulation loop ###

metrics = []

if MODE == 'simulate':
    print("=== SIMULATION MODE ===")
    for i in tqdm(range(repetitions)):
        # simulate the ground truth tree and reconstruct it using different solvers
        gt_tree = simulate_gt(tree_config)
        print("next gt tree ... ")
        reconstructions = reconstruct(gt_tree, solvers)

        # calculate the metrics for the ground truth tree and the reconstructed trees
        gt_metrics, reconstructions_metrics = calculate_metrics(gt_tree, reconstructions)

        # append the new metrics to the existing ones
        metrics.append({
            'gt_metrics': gt_metrics,
            'reconstructions_metrics': reconstructions_metrics
        })

elif MODE == 'load_pipeline':
    print("=== PIPELINE LOADING MODE ===")

    # Load pipeline objects
    gt_trees, cas9_instances = load_pipeline_objects(PIPELINE_DIR)

    print(f"Processing {len(cas9_instances)} CAS9 instances with {len(solvers)} solvers...")

    for i, (gt_tree, cas9_tree) in enumerate(zip(gt_trees, cas9_instances)):
        print(f"\nProcessing instance {i+1}/{len(cas9_instances)}")

        # Reconstruct from the CAS9 instance (which has the recording data)
        reconstructions = reconstruct_from_cas9_instance(cas9_tree, solvers)

        # Calculate metrics using the original GT tree for comparison
        gt_metrics, reconstructions_metrics = calculate_metrics(gt_tree, reconstructions)

        # append the new metrics to the existing ones
        metrics.append({
            'gt_metrics': gt_metrics,
            'reconstructions_metrics': reconstructions_metrics
        })

    # Update repetitions to match actual number of instances processed
    repetitions = len(cas9_instances)

elif MODE == 'load_results':
    print("=== PIPELINE RESULTS LOADING MODE ===")

    # Load existing pipeline results from parquet files
    pipeline_df = load_pipeline_results(PIPELINE_DIR)

    # Filter by requested solvers if specified
    available_solvers = set(pipeline_df['solver'].unique())
    requested_solvers = set(args.solvers)

    if not requested_solvers.issubset(available_solvers):
        missing = requested_solvers - available_solvers
        print(f"Warning: Requested solvers not found in pipeline results: {missing}")
        print(f"Available solvers: {available_solvers}")

    # Filter to only requested solvers
    filtered_df = pipeline_df[pipeline_df['solver'].isin(args.solvers)]

    print(f"Using {len(filtered_df)} reconstructions from {len(args.solvers)} solvers")

    # Directly create standardized results from pipeline data - following pipeline column order
    reconstructions_scores = []
    for _, row in filtered_df.iterrows():
        reconstructions_scores.append({
            # 1. reconstruction_id (follows pipeline order)
            'reconstruction_id': row['reconstruction_id'],
            # 2-4. Instance/repetition info
            'gt_instance_id': row['gt_instance_id'],
            'cas9_simulation_id': row.get('cas9_simulation_id', 0),
            'reconstruction_num': row.get('reconstruction_num', 0),
            # 5-6. CAS9 tier info
            'cas9_tier': row['cas9_tier'],
            'cas9_tier_name': row.get('cas9_tier_name', 'SimulationPHS'),
            # 7-8. Recording configuration
            'recording_sites': row.get('recording_sites', 'N/A'),
            'states_per_site': row.get('states_per_site', 'N/A'),
            # 9. Solver
            'solver': row['solver'],
            # 10-11. Timing and run info
            'computation_time_seconds': row.get('computation_time_seconds', 'N/A'),
            'run_name': row.get('run_name', 'N/A'),
            # 12. Experiment ID
            'experiment_id': row.get('experiment_id', 'N/A'),
            # 13-14. Tree sizes
            'gt_tree_size': row.get('gt_tree_size', 'N/A'),
            'sampled_tree_size': row.get('sampled_tree_size', 'N/A'),
            # 15-16. Simulation parameters
            'lam_simulation': row.get('lam_simulation', 'N/A'),
            'q_simulation': row.get('q_simulation', 'N/A'),
            # 17-18. Mutation rates
            'proportion_mutated_simulation': row.get('proportion_mutated_simulation', 'N/A'),
            'q_simulation_source': row.get('q_simulation_source', 'N/A'),
            # 19-21. GT parameters
            'lam_gt': row.get('lam_gt', 'N/A'),
            'q_gt': row.get('q_gt', 'N/A'),
            'proportion_mutated_gt': row.get('proportion_mutated_gt', 'N/A'),
            # 22-23. PHS source info
            'phs_lam_source': row.get('phs_lam_source', 'N/A'),
            'phs_q_source': row.get('phs_q_source', 'N/A'),
            # 24-25. Distance metrics
            'triplets_distance': row.get('triplets_distance', 0),
            'RF_distance': row.get('RF_distance', 0),
            # 26-29. Parsimony info
            'parsimony_total_mutations': row.get('parsimony_total_mutations', 0),
            'parsimony_computation_time_ms': row.get('parsimony_computation_time_ms', 'N/A'),
            'parsimony_method_used': row.get('parsimony_method_used', 'N/A'),
            'parsimony_internal_states_inferred': row.get('parsimony_internal_states_inferred', True),
            # 30-31. cPHS
            'cPHS': row.get('cPHS', 0),
            'cPHS_gt': row.get('cPHS_gt', 0),
            # 32. Likelihood timing
            'likelihood_computation_time_ms': row.get('likelihood_computation_time_ms', 'N/A'),
            # 33-37. Worker info
            'reconstructed_tree_path': row.get('reconstructed_tree_path', 'N/A'),
            'worker_cas9_instance_path': row.get('worker_cas9_instance_path', 'N/A'),
            'worker_solver': row.get('worker_solver', row['solver']),
            'worker_tier': row.get('worker_tier', row['cas9_tier']),
            'worker_tier_name': row.get('worker_tier_name', 'SimulationPHS'),
            # 38-41. Final metrics
            'parsimony_score': row.get('parsimony_score', row.get('parsimony_total_mutations', 0)),
            'log_likelihood': row.get('log_likelihood', 0),
            'log_likelihood_simulation': row.get('log_likelihood_simulation', 0),
            'log_likelihood_gt': row.get('log_likelihood_gt', 0),

            # Additional fields for compatibility with simulation mode
            'repetition': row['gt_instance_id'] + 1,
            'gt_instance': f"gt_{row['gt_instance_id'] + 1}",
            # Standardized names (aliases)
            'parsimony statistic': row.get('parsimony_score', row.get('parsimony_total_mutations', 0)),
            'log-likelihood statistic': -row.get('log_likelihood', 0),  # Convert to positive
            'cPHS statistic (python)': row.get('cPHS', 0),
            'cPHS statistic (stellars)': row.get('cPHS', 0),
            'rf distance': row.get('RF_distance', 0),
            'triplets distance': row.get('triplets_distance', 0),
            'parsimony distance': 0,  # Not available from pipeline
            'likelihood distance': 0,  # Not available from pipeline
            'config_N': row.get('gt_tree_size', 'N/A'),
            'config_n': row.get('sampled_tree_size', 'N/A'),
            'config_k': row.get('recording_sites', 'N/A'),
            'config_m': row.get('states_per_site', 'N/A'),
            'config_rho': row.get('proportion_mutated_simulation', 'N/A'),
            'config_state_priors_exp': 'from_pipeline',
        })

    # Update repetitions and solvers to match actual data
    repetitions = len(set(row['gt_instance_id'] for _, row in filtered_df.iterrows()))
    solvers = sorted(filtered_df['solver'].unique())

    # Skip the normal metrics processing loop for load_results mode
    metrics = []  # Empty since we've already processed results above

else:
    raise ValueError(f"Invalid MODE: {MODE}. Must be 'simulate', 'load_pipeline', or 'load_results'")


# %%
### output results ###

# prepare the results for output with solver metadata
# For load_results mode, reconstructions_scores is already populated above
if MODE != 'load_results':
    reconstructions_scores = []

for rep_idx, m in enumerate(metrics):
    gt, rs = m['gt_metrics'], m['reconstructions_metrics']
        # rs is a list of reconstructrions (one for each solver), all correspond to the same gt
    for solver_idx, r in enumerate(rs):
        # Handle different solver naming approaches
        if MODE == 'load_results':
            # For load_results mode, get solver name from pipeline data
            if 'instance_data' in m and len(m['instance_data']) > solver_idx:
                solver_name = m['instance_data'].iloc[solver_idx]['solver']
            else:
                solver_name = f"solver_{solver_idx}"
        else:
            # For simulate and load_pipeline modes
            solver_name = solvers[solver_idx]
        # Handle reconstruction_id based on mode
        if MODE == 'load_results' and 'instance_data' in m and len(m['instance_data']) > solver_idx:
            # Use original reconstruction_id from pipeline
            reconstruction_id = m['instance_data'].iloc[solver_idx].get('reconstruction_id', f"{solver_name}_rep{rep_idx + 1}")
        else:
            # Generate reconstruction_id for simulate and load_pipeline modes
            reconstruction_id = f"{solver_name}_rep{rep_idx + 1}"

        reconstructions_scores.append(
            {
                'repetition': rep_idx + 1,
                'solver': solver_name,
                'gt_instance': f"gt_{rep_idx + 1}",
                'cas9_tier': 1,  # Single tier in this simulation
                'reconstruction_id': reconstruction_id,
                'parsimony statistic': r['parsimony'],
                'log-likelihood statistic': -r['likelihood'],
                'cPHS statistic (python)': r['cPHS'],
                'cPHS statistic (stellars)': r['cPHS-stellars'],
                'rf distance': r['rf'],
                'triplets distance': 1 - r['triplets'],
                'parsimony distance': np.max([0, r['parsimony'] / gt['parsimony'] - 1]),
                'likelihood distance': np.max([0, r['likelihood'] / gt['likelihood'] - 1]),
            }
        )

        # Add configuration metadata based on mode
        if MODE == 'simulate':
            # Use the simulation config
            reconstructions_scores[-1].update({
                'config_N': int(tree_config['N']),
                'config_n': int(tree_config['n']),
                'config_k': tree_config['k'],
                'config_m': tree_config['m'],
                'config_rho': tree_config['rho'],
                'config_state_priors_exp': tree_config['state_priors_exponents'],
            })
        elif MODE == 'load_pipeline':
            # Infer config from loaded objects (using first GT tree for this rep)
            current_gt = gt_trees[rep_idx] if 'gt_trees' in locals() else None
            if current_gt is not None:
                cm = current_gt.character_matrix
                reconstructions_scores[-1].update({
                    'config_N': 'N/A',  # Original tree size not available
                    'config_n': len(cm),
                    'config_k': len(cm.columns),
                    'config_m': cm.stack().nunique(),
                    'config_rho': 'inferred',
                    'config_state_priors_exp': 'inferred',
                })
        elif MODE == 'load_results':
            # Use actual config from pipeline results
            if 'instance_data' in m and len(m['instance_data']) > 0:
                pipeline_row = m['instance_data'].iloc[solver_idx]
                reconstructions_scores[-1].update({
                    'config_N': pipeline_row.get('gt_tree_size', 'N/A'),
                    'config_n': pipeline_row.get('sampled_tree_size', 'N/A'),
                    'config_k': pipeline_row.get('recording_sites', 'N/A'),
                    'config_m': pipeline_row.get('states_per_site', 'N/A'),
                    'config_rho': pipeline_row.get('proportion_mutated_simulation', 'N/A'),
                    'config_state_priors_exp': 'from_pipeline',
                    # Add pipeline-specific metadata
                    'lam_simulation': pipeline_row.get('lam_simulation', 'N/A'),
                    'q_simulation': pipeline_row.get('q_simulation', 'N/A'),
                    'lam_gt': pipeline_row.get('lam_gt', 'N/A'),
                    'q_gt': pipeline_row.get('q_gt', 'N/A'),
                    'computation_time_seconds': pipeline_row.get('computation_time_seconds', 'N/A'),
                    'run_name': pipeline_row.get('run_name', 'N/A'),
                })

# display results
results_df = pd.DataFrame(reconstructions_scores)

# Standardize columns for comparison (skip for load_results mode as it's already standardized)
if MODE != 'load_results':
    results_df = standardize_columns(results_df, MODE)

# Add relative parsimony scores (solver parsimony / minimum solver parsimony for each tree instance)
results_df = add_relative_parsimony(results_df, MODE)

print("\n" + "="*80)
print(f"SIMULATION_PHS_CONTROL.PY RESULTS - {MODE.upper()} MODE")
print("="*80)

# Show configuration summary
print(f"\nConfiguration Summary:")
print(f"- Mode: {MODE}")
if MODE == 'load_pipeline':
    print(f"- Pipeline directory: {PIPELINE_DIR}")
print(f"- Repetitions/Instances: {repetitions}")
print(f"- Solvers: {solvers}")

# For simulation mode, show tree_config. For pipeline mode, infer from loaded data
if MODE == 'simulate':
    print(f"- Tree size (N): {int(tree_config['N'])}")
    print(f"- Subsampled size (n): {int(tree_config['n'])}")
    print(f"- Cassettes (k): {tree_config['k']}")
    print(f"- States (m): {tree_config['m']}")
    print(f"- Mutation probability (rho): {tree_config['rho']}")
    print(f"- State priors exponent: {tree_config['state_priors_exponents']}")
elif MODE == 'load_pipeline' and len(results_df) > 0:
    # Infer parameters from the first result
    first_result = results_df.iloc[0]
    if 'config_N' in first_result:
        print(f"- Tree size (N): {first_result.get('config_N', 'N/A')}")
        print(f"- Subsampled size (n): {first_result.get('config_n', 'N/A')}")
        print(f"- Cassettes (k): {first_result.get('config_k', 'N/A')}")
        print(f"- States (m): {first_result.get('config_m', 'N/A')}")
        print(f"- Mutation probability (rho): {first_result.get('config_rho', 'N/A')}")
        print(f"- State priors exponent: {first_result.get('config_state_priors_exp', 'N/A')}")
    else:
        print("- Configuration details: Available from loaded pipeline objects")

# Show detailed results grouped by solver
print(f"\nDetailed Results by Solver:")
print("-" * 80)
for solver in solvers:
    solver_data = results_df[results_df['solver'] == solver]
    print(f"\n{solver.upper()} ({len(solver_data)} reconstructions):")
    print(f"  RF distance:        {solver_data['rf distance'].mean():.4f} ± {solver_data['rf distance'].std():.4f}")
    print(f"  Parsimony:          {solver_data['parsimony statistic'].mean():.2f} ± {solver_data['parsimony statistic'].std():.2f}")
    if 'relative_parsimony' in solver_data.columns:
        print(f"  Relative parsimony: {solver_data['relative_parsimony'].mean():.4f} ± {solver_data['relative_parsimony'].std():.4f}")
    print(f"  Log-likelihood:     {solver_data['log-likelihood statistic'].mean():.2f} ± {solver_data['log-likelihood statistic'].std():.2f}")
    print(f"  cPHS (python):      {solver_data['cPHS statistic (python)'].mean():.6f} ± {solver_data['cPHS statistic (python)'].std():.6f}")
    print(f"  cPHS (stellars):    {solver_data['cPHS statistic (stellars)'].mean():.6f} ± {solver_data['cPHS statistic (stellars)'].std():.6f}")
    print(f"  Triplets distance:  {solver_data['triplets distance'].mean():.4f} ± {solver_data['triplets distance'].std():.4f}")

print(f"\nOverall Summary Statistics:")
print(f"- Total reconstructions: {len(results_df)}")
print(f"- Mean RF distance: {results_df['rf distance'].mean():.4f}")
print(f"- Mean parsimony statistic: {results_df['parsimony statistic'].mean():.2f}")
if 'relative_parsimony' in results_df.columns:
    print(f"- Mean relative parsimony: {results_df['relative_parsimony'].mean():.4f}")
print(f"- Mean log-likelihood statistic: {results_df['log-likelihood statistic'].mean():.2f}")
print(f"- Mean cPHS (python): {results_df['cPHS statistic (python)'].mean():.6f}")
print(f"- Mean cPHS (stellars): {results_df['cPHS statistic (stellars)'].mean():.6f}")
print(f"- Mean triplets distance: {results_df['triplets distance'].mean():.4f}")

# Show comparison table format similar to cascade workflow
print(f"\nComparison Table (similar to output/sim format):")
print("-" * 120)
comparison_cols = ['repetition', 'solver', 'parsimony statistic', 'log-likelihood statistic',
                  'cPHS statistic (stellars)', 'rf distance', 'triplets distance']
if 'relative_parsimony' in results_df.columns:
    comparison_cols.insert(3, 'relative_parsimony')  # Add after parsimony statistic
print(results_df[comparison_cols].to_string(index=False, float_format='%.6f'))

# Save results to parquet file
output_filename = f"{args.output_prefix}_{MODE}_{repetitions}reps.parquet"
try:
    results_df.to_parquet(output_filename, index=False)
    print(f"\n✓ Results saved to: {output_filename}")
    print(f"  File contains {len(results_df)} rows with {len(results_df.columns)} columns")
    print(f"  Columns: {list(results_df.columns)}")
except Exception as e:
    print(f"\n⚠ Warning: Failed to save parquet file: {e}")
    # Fallback to CSV if parquet fails
    csv_filename = f"simulation_phs_control_results_{repetitions}reps.csv"
    results_df.to_csv(csv_filename, index=False)
    print(f"  Saved as CSV instead: {csv_filename}")

# Also create a styled version for Jupyter if running interactively
if 'get_ipython' in globals():
    styled_results = results_df.style.set_table_styles({
        'cPHS statistic (stellars)': [{'selector': 'td', 'props': [('border-right', '3px solid black')]},
                {'selector': 'th', 'props': [('border-right', '3px solid black')]}]
    })
    display(styled_results)

# %%
