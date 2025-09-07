"""
Stellars Usage Guide: Choosing the Right Function

This module provides guidance on when to use convenience functions vs. explicit parameter functions,
and how to avoid common pitfalls with parameter estimation.
"""

# ============================================================================
# FUNCTION SELECTION GUIDE
# ============================================================================

"""
1. DISTANCE CALCULATIONS - Generally Safe for Convenience
==========================================================

✅ SAFE TO USE SIMPLE FUNCTIONS:
    - stellars.rf_distance_simple()
    - stellars.triplets_distance_simple() 
    - stellars.triplets_correctness()

These functions have minimal parameters and are unlikely to fail due to auto-estimation.

Example:
    rf_dist = stellars.rf_distance_simple(tree1, tree2)
    triplets_dist = stellars.triplets_distance_simple(tree1, tree2)


2. PHS CALCULATIONS - CAUTION REQUIRED
=====================================

⚠️  DANGEROUS: stellars.phs_simple() without explicit parameters
❌  Can silently produce incorrect results due to auto-parameter estimation

✅  SAFE: stellars.phs() with explicit parameters
✅  SAFE: stellars.phs_from_cassiopeia() with explicit parameters

RECOMMENDED APPROACH:
    # Calculate parameters manually for reproducibility
    lam, q = calculate_mutation_parameters(tree, char_matrix)
    
    # Use explicit parameters
    result = stellars.phs(
        tree_newick=tree_newick,
        character_matrix=char_matrix,
        mutation_rate=lam,
        collision_probability=q,
        internal_states=internal_states,  # If available
        leaf_names=leaf_names
    )

AVOID:
    # This may give inconsistent results
    p_value = stellars.phs_simple(tree_newick, char_matrix)  # Auto-estimation!


3. PARSIMONY & LIKELIHOOD - Generally Safe
==========================================

✅ SAFE TO USE SIMPLE FUNCTIONS:
    - stellars.parsimony_simple()
    - stellars.likelihood_simple()

These typically don't require critical parameter estimation.

Example:
    score = stellars.parsimony_simple(tree_newick, char_matrix)
    likelihood = stellars.likelihood_simple(tree_newick, char_matrix)
"""


def calculate_mutation_parameters(tree, known_priors=True):
    """
    Manually calculate mutation parameters for PHS calculation.
    
    This is the same approach used in the working reconstruction script.
    Use this when you need explicit control over parameters.
    
    Args:
        tree: Cassiopeia tree object
        known_priors: Whether the tree has known state priors
        
    Returns:
        tuple: (mutation_rate, collision_probability)
        
    Example:
        >>> lam, q = calculate_mutation_parameters(reconstructed_tree)
        >>> result = stellars.phs(tree_newick, char_matrix, 
        ...                       mutation_rate=lam, collision_probability=q)
    """
    import numpy as np
    
    # Estimate mutation rate from tree leaves
    cm = tree.character_matrix.to_numpy()
    proportion_mutated = np.sum(cm > 0) / np.sum(cm >= 0)
    lam = -np.log(1.0 - proportion_mutated)

    # Get collision probability for cPHS
    if known_priors and hasattr(tree, 'priors') and tree.priors:
        priors = np.array(list(tree.priors[0].values()))
        q = np.sum(priors ** 2)
    else:
        q = 1.0 / tree.character_matrix.stack().nunique()

    return lam, q


def extract_internal_states(tree):
    """
    Extract internal character states from a tree.
    
    Use this when you need explicit control over internal states for PHS calculation.
    
    Args:
        tree: Cassiopeia tree object
        
    Returns:
        dict: Mapping of internal node names to character states
        
    Example:
        >>> internal_states = extract_internal_states(reconstructed_tree)
        >>> result = stellars.phs(tree_newick, char_matrix,
        ...                       internal_states=internal_states)
    """
    internal_character_states = {}
    
    for internal_node in tree.internal_nodes:
        node_name = str(internal_node)
        internal_states = tree.get_character_states(internal_node)
        try:
            if hasattr(internal_states, 'tolist'):
                internal_character_states[node_name] = internal_states.tolist()
            else:
                internal_character_states[node_name] = list(internal_states)
        except Exception:
            # If any extraction fails, return empty dict (use Fitch parsimony instead)
            return {}
    
    return internal_character_states


# ============================================================================
# USAGE EXAMPLES
# ============================================================================

def example_safe_phs_calculation():
    """
    Example of safe PHS calculation with explicit parameters.
    This approach avoids the auto-parameter estimation trap.
    """
    # Assuming you have: tree, character_matrix
    tree = None  # Your Cassiopeia tree
    
    # Step 1: Calculate parameters explicitly
    lam, q = calculate_mutation_parameters(tree)
    
    # Step 2: Extract internal states if needed
    internal_states = extract_internal_states(tree)
    
    # Step 3: Get tree structure
    tree_newick = tree.get_newick(record_branch_lengths=True, record_node_names=True)
    character_matrix = tree.character_matrix.values.astype(int).tolist()
    leaf_names = list(tree.character_matrix.index)
    
    # Step 4: Call stellars with explicit parameters
    import stellars
    result = stellars.phs(
        tree_newick=tree_newick,
        character_matrix=character_matrix,
        internal_states=internal_states,
        mutation_rate=lam,
        collision_probability=q,
        leaf_names=leaf_names,
        use_internal_states=True
    )
    
    print(f"cPHS p-value: {result.phs_score:.2e}")
    print(f"Parameters used: λ={lam:.4f}, q={q:.6f}")
    print(f"Internal states: {len(internal_states)} nodes")
    
    return result


def example_convenience_with_warnings():
    """
    Example of using convenience functions with proper warnings enabled.
    """
    import stellars
    
    # This will warn about parameter estimation
    result = stellars.phs_from_cassiopeia(
        tree,
        warn_on_estimation=True,     # Enable warnings
        validate_parameters=True     # Enable validation
    )
    
    # Check what parameters were estimated
    params = result.get('parameters_used', {})
    if params.get('mutation_rate_estimated'):
        print(f"⚠️  Mutation rate was estimated: {params['mutation_rate']:.4f}")
    if params.get('collision_probability_estimated'):
        print(f"⚠️  Collision probability was estimated: {params['collision_probability']:.6f}")
    
    print(f"cPHS p-value: {result['phs_score']:.2e}")


# ============================================================================
# MIGRATION GUIDE
# ============================================================================

"""
MIGRATING FROM PROBLEMATIC CODE:

❌ OLD (Dangerous):
    p_value = stellars.phs_simple(tree_newick, char_matrix)

✅ NEW (Safe):
    # Option 1: Explicit parameters (Recommended)
    lam, q = calculate_mutation_parameters(tree)
    result = stellars.phs(tree_newick, char_matrix, mutation_rate=lam, collision_probability=q)
    p_value = result.phs_score
    
    # Option 2: Convenience with warnings
    result = stellars.phs_from_cassiopeia(tree, warn_on_estimation=True)
    p_value = result['phs_score']

KEY PRINCIPLE: 
When in doubt, provide explicit parameters rather than relying on auto-estimation.
This ensures reproducible, reliable results.
"""