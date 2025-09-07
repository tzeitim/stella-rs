
from typing import Optional, Dict, Any, Union
import warnings

# Import low-level Rust functions (for backward compatibility)
from stellars.stellars import (
        triplets_correct, 
        triplets_correct_fast,
        triplets_correct_ultra,
        triplets_correct_parallel,
        triplets_correct_optimized,
        phs_optimized,
        robinson_foulds_distance,
        robinson_foulds_distance_normalized,
        parsimony_score,
        parsimony_p_value,
        likelihood_score,
        likelihood_distance,
        init_logging,
        )

# Import new structured results and high-level API
from .results import (
    RFDistanceResult, 
    TripletsResult, 
    PHSResult, 
    ParsimonyResult, 
    LikelihoodResult
)

from .distances import (
    rf_distance,
    rf_distance_simple,
    triplets_distance,
    triplets_distance_simple,
    triplets_correctness,
    robinson_foulds_distance as robinson_foulds,  # Backward compatibility alias
    robinson_foulds_simple
)

from .metrics import (
    phs,
    phs_simple,
    phs_significance_test,
    parsimony,
    parsimony_simple,
    likelihood,
    likelihood_simple,
    likelihood_distance
)

# Import usage guidance functions
from .usage_guide import (
    calculate_mutation_parameters,
    extract_internal_states,
    example_safe_phs_calculation
)

# Legacy aliases for backward compatibility
rf_distance_normalized = robinson_foulds_distance_normalized
triplets_correct_distance = triplets_distance


# Convenience functions for Cassiopeia trees
def phs_from_cassiopeia(
    tree,
    mutation_rate: Optional[float] = None,
    collision_probability: Optional[float] = None,
    missing_state: int = -1,
    unedited_state: int = 0,
    use_provided_internal_states: bool = True,
    max_threads: Optional[int] = None,
    reconstruct_ancestral: bool = False,
    compute_branch_lengths: bool = True,
    minimum_branch_length: Optional[float] = None,
    validate_parameters: bool = True,
    warn_on_estimation: bool = True
) -> Dict[str, Any]:
    """
    Calculate corrected Pairwise Homoplasy Score (cPHS) directly from a Cassiopeia tree.
    
    This convenience function automatically extracts all necessary data from a Cassiopeia
    tree object and computes the cPHS score using Stellars' optimized implementation.
    
    Args:
        tree: Cassiopeia tree object with character matrix and internal states
        mutation_rate: Mutation rate parameter (λ). If None, will be estimated automatically
        collision_probability: Collision probability parameter (q). If None, will be estimated automatically  
        missing_state: Value representing missing/unknown character states (default: -1)
        unedited_state: Value representing unedited/ancestral character states (default: 0)
        use_provided_internal_states: Whether to use tree's internal states or reconstruct with Fitch parsimony (default: True)
        max_threads: Maximum number of threads for parallel computation (default: None = auto)
        reconstruct_ancestral: Whether to call tree.reconstruct_ancestral_characters() before analysis (default: True)
        compute_branch_lengths: Whether to compute branch lengths using ConvexML after ancestral reconstruction (default: True)
        minimum_branch_length: Minimum branch length for ConvexML optimization. If None, attempts to estimate from tree (default: None)
        
    Returns:
        Dict containing:
            - phs_score: The corrected pairwise homoplasy p-value
            - total_pairs: Number of leaf pairs analyzed  
            - computation_time_ms: Time taken for computation
            - method_used: Method identifier
            - parallel_chunks_used: Number of parallel chunks used
            
    Example:
        >>> import cassiopeia as cass
        >>> import stellars
        >>> 
        >>> # Load or create your Cassiopeia tree
        >>> tree = cass.data.load_tree_from_pickle("my_tree.pkl")
        >>> 
        >>> # Compute cPHS score with automatic parameter estimation
        >>> result = stellars.phs_from_cassiopeia(tree)
        >>> print(f"cPHS p-value: {result['phs_score']:.2e}")
        >>> 
        >>> # Or specify parameters explicitly
        >>> result = stellars.phs_from_cassiopeia(
        ...     tree, 
        ...     mutation_rate=0.3, 
        ...     collision_probability=0.1
        ... )
    """
    
    # Validate input
    if not hasattr(tree, 'character_matrix'):
        raise ValueError("Tree object must have a 'character_matrix' attribute")
    if not hasattr(tree, 'get_newick'):
        raise ValueError("Tree object must have a 'get_newick' method")
        
    # Handle ancestral character reconstruction
    if reconstruct_ancestral:
        if hasattr(tree, 'reconstruct_ancestral_characters'):
            # Check if ancestral characters are already reconstructed
            try:
                # Try to get internal nodes - if they have character states, they're already reconstructed
                internal_nodes = getattr(tree, 'internal_nodes', [])
                if internal_nodes:
                    # Check if first internal node has character states
                    try:
                        first_internal = internal_nodes[0]
                        existing_states = tree.get_character_states(first_internal)
                        if existing_states is not None and len(existing_states) > 0:
                            warnings.warn(
                                "Tree already appears to have reconstructed ancestral characters. "
                                "reconstruct_ancestral=True will overwrite existing reconstructions."
                            )
                    except:
                        pass  # If we can't check, just proceed with reconstruction
                
                # Perform reconstruction
                tree.reconstruct_ancestral_characters()
            except Exception as e:
                warnings.warn(f"Failed to reconstruct ancestral characters: {e}")
        else:
            warnings.warn("Tree does not have 'reconstruct_ancestral_characters' method, skipping ancestral reconstruction")
    
    # Compute branch lengths using ConvexML if requested
    if compute_branch_lengths:
        try:
            import convexml
            
            # Get current tree newick and leaf sequences
            current_newick = tree.get_newick(record_branch_lengths=True, record_node_names=True)
            leaf_sequences = {leaf: tree.get_character_states(leaf) for leaf in tree.leaves}
            
            # Estimate minimum branch length if not provided
            min_branch_len = minimum_branch_length
            if min_branch_len is None:
                try:
                    # Try to estimate from existing branch lengths
                    branch_lengths = []
                    for parent in tree.nodes:
                        for child in tree.children(parent):
                            bl = tree.get_branch_length(parent, child)
                            if bl is not None and bl > 0:
                                branch_lengths.append(bl)
                    
                    if branch_lengths:
                        min_branch_len = min(branch_lengths) * 0.1  # Use 10% of minimum as safety margin
                    else:
                        min_branch_len = 1e-6  # Default fallback
                except:
                    min_branch_len = 1e-6  # Default fallback
            
            # Call ConvexML to optimize branch lengths
            convexml_result = convexml.convexml(
                tree_newick=current_newick,
                leaf_sequences=leaf_sequences,
                minimum_branch_length=min_branch_len
            )
            
            # Update the tree with optimized branch lengths
            optimized_newick = convexml_result["tree_newick"]
            
            # Create new tree object with optimized branch lengths
            import cassiopeia as cas
            optimized_tree = cas.data.CassiopeiaTree(
                character_matrix=tree.character_matrix,
                tree=optimized_newick,
                missing_state_indicator=getattr(tree, 'missing_state_indicator', missing_state)
            )
            
            # Replace the original tree reference for subsequent operations
            tree = optimized_tree
            
        except ImportError:
            warnings.warn("ConvexML not available. Install with 'pip install convexml' to enable branch length optimization. Skipping compute_branch_lengths.")
        except Exception as e:
            warnings.warn(f"Failed to compute branch lengths with ConvexML: {e}. Continuing with original branch lengths.")
        
    # Extract tree structure with branch lengths and node names
    try:
        tree_newick = tree.get_newick(record_branch_lengths=True, record_node_names=True)
    except Exception as e:
        raise ValueError(f"Failed to extract Newick string from tree: {e}")
    
    # Extract character matrix and ensure proper ordering
    try:
        character_matrix = tree.character_matrix.values.astype(int).tolist()
        leaf_names = list(tree.character_matrix.index)  # Critical: preserves DataFrame order
    except Exception as e:
        raise ValueError(f"Failed to extract character matrix from tree: {e}")
        
    # Validate character matrix structure
    if not character_matrix:
        raise ValueError("Character matrix is empty")
    if not leaf_names:
        raise ValueError("No leaf names found in character matrix")
    if len(character_matrix) != len(leaf_names):
        raise ValueError(f"Character matrix rows ({len(character_matrix)}) don't match leaf names ({len(leaf_names)})")
    
    # Extract internal character states if requested
    internal_character_states = {}
    if use_provided_internal_states:
        try:
            # Get all internal nodes from the tree
            internal_nodes = getattr(tree, 'internal_nodes', [])
            if not internal_nodes:
                # Fallback: try to get internal nodes differently
                all_nodes = getattr(tree, 'nodes', [])
                if all_nodes:
                    # Filter to internal nodes (nodes that aren't leaves)
                    leaf_set = set(leaf_names)
                    internal_nodes = [node for node in all_nodes if str(node) not in leaf_set]
                
            for internal_node in internal_nodes:
                node_name = str(internal_node)
                try:
                    internal_states = tree.get_character_states(internal_node)
                    if hasattr(internal_states, 'tolist'):
                        internal_character_states[node_name] = internal_states.tolist()
                    else:
                        internal_character_states[node_name] = list(internal_states)
                except Exception as node_error:
                    warnings.warn(f"Could not extract character states for internal node '{node_name}': {node_error}")
                    
        except Exception as e:
            warnings.warn(f"Failed to extract internal character states: {e}. Will use Fitch parsimony instead.")
            use_provided_internal_states = False
            internal_character_states = {}
    
    # Validate and estimate parameters if needed
    original_mutation_rate = mutation_rate
    original_collision_probability = collision_probability
    
    if mutation_rate is None:
        # Estimate mutation rate from character matrix
        cm = tree.character_matrix.to_numpy()
        proportion_mutated = np.sum(cm > 0) / np.sum(cm >= 0)
        mutation_rate = -np.log(1.0 - proportion_mutated) if proportion_mutated < 1.0 else 1.0
        
        if warn_on_estimation:
            warnings.warn(f"Mutation rate not provided, estimated as {mutation_rate:.4f} from character matrix. "
                         f"For reproducible results, consider providing explicit mutation_rate parameter.")
    
    if collision_probability is None:
        # Estimate collision probability from priors if available
        if hasattr(tree, 'priors') and tree.priors:
            try:
                priors = np.array(list(tree.priors[0].values()))
                collision_probability = np.sum(priors ** 2)
            except Exception:
                collision_probability = 1.0 / tree.character_matrix.stack().nunique()
        else:
            collision_probability = 1.0 / tree.character_matrix.stack().nunique()
            
        if warn_on_estimation:
            warnings.warn(f"Collision probability not provided, estimated as {collision_probability:.6f}. "
                         f"For reproducible results, consider providing explicit collision_probability parameter.")
    
    # Validate estimated parameters
    if validate_parameters:
        if not (0.0 < mutation_rate <= 10.0):  # Reasonable range for mutation rate
            warnings.warn(f"Mutation rate {mutation_rate:.4f} seems unusually {'high' if mutation_rate > 10.0 else 'low'}. "
                         f"This may indicate parameter estimation issues.")
        
        if not (0.0 < collision_probability <= 1.0):
            warnings.warn(f"Collision probability {collision_probability:.6f} is outside valid range [0, 1]. "
                         f"This may indicate parameter estimation issues.")
        
        # Check for potentially problematic scenarios
        if len(internal_character_states) == 0 and use_provided_internal_states:
            warnings.warn("No internal character states found, but use_provided_internal_states=True. "
                         "This may lead to suboptimal results. Consider setting use_provided_internal_states=False "
                         "or ensuring ancestral characters are reconstructed.")

    # Call the optimized PHS function
    try:
        result = phs_optimized(
            tree_newick=tree_newick,
            character_matrix=character_matrix,
            internal_character_states=internal_character_states,
            mutation_rate=mutation_rate,
            collision_probability=collision_probability,
            missing_state=missing_state,
            unedited_state=unedited_state,
            use_provided_internal_states=use_provided_internal_states,
            leaf_names=leaf_names,  # Critical: ensures correct character mapping
            max_threads=max_threads
        )
        
        # Add parameter info to result for transparency
        result['parameters_used'] = {
            'mutation_rate': mutation_rate,
            'collision_probability': collision_probability,
            'mutation_rate_estimated': original_mutation_rate is None,
            'collision_probability_estimated': original_collision_probability is None,
            'internal_states_count': len(internal_character_states)
        }
        
        return result
        
    except Exception as e:
        raise RuntimeError(f"Stellars PHS calculation failed: {e}")


def quick_phs(
    tree,
    **kwargs
) -> float:
    """
    Quick convenience function that returns just the cPHS p-value from a Cassiopeia tree.
    
    This is a simplified wrapper around phs_from_cassiopeia() that returns only the
    p-value for cases where you just need the statistical result.
    
    Args:
        tree: Cassiopeia tree object
        **kwargs: Additional arguments passed to phs_from_cassiopeia()
        
    Returns:
        float: The corrected pairwise homoplasy p-value
        
    Example:
        >>> import cassiopeia as cass
        >>> import stellars
        >>> 
        >>> tree = cass.data.load_tree_from_pickle("my_tree.pkl")
        >>> p_value = stellars.quick_phs(tree, mutation_rate=0.3, collision_probability=0.1)
        >>> 
        >>> if p_value < 0.05:
        ...     print(f"Significant homoplasies detected (p = {p_value:.2e})")
        ... else:
        ...     print(f"No significant homoplasies (p = {p_value:.3f})")
    """
    result = phs_from_cassiopeia(tree, **kwargs)
    return result['phs_score']


def batch_phs_analysis(
    trees: list,
    tree_names: Optional[list] = None,
    **kwargs
) -> Dict[str, Dict[str, Any]]:
    """
    Analyze multiple Cassiopeia trees and return cPHS results for each.
    
    Args:
        trees: List of Cassiopeia tree objects
        tree_names: Optional list of names for the trees. If None, will use indices.
        **kwargs: Additional arguments passed to phs_from_cassiopeia()
        
    Returns:
        Dict mapping tree names to their PHS analysis results
        
    Example:
        >>> trees = [tree1, tree2, tree3]
        >>> names = ["reconstruction_1", "reconstruction_2", "reconstruction_3"]
        >>> results = stellars.batch_phs_analysis(trees, names, mutation_rate=0.3)
        >>> 
        >>> for name, result in results.items():
        ...     print(f"{name}: p-value = {result['phs_score']:.2e}")
    """
    if tree_names is None:
        tree_names = [f"tree_{i}" for i in range(len(trees))]
    
    if len(trees) != len(tree_names):
        raise ValueError("Number of trees must match number of tree names")
    
    results = {}
    for tree, name in zip(trees, tree_names):
        try:
            result = phs_from_cassiopeia(tree, **kwargs)
            results[name] = result
        except Exception as e:
            warnings.warn(f"Failed to analyze tree '{name}': {e}")
            results[name] = {"error": str(e), "phs_score": None}
    
    return results


def phs_significance_test(
    tree,
    alpha: float = 0.05,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform a significance test for homoplasies in a Cassiopeia tree.
    
    Args:
        tree: Cassiopeia tree object
        alpha: Significance threshold (default: 0.05)
        **kwargs: Additional arguments passed to phs_from_cassiopeia()
        
    Returns:
        Dict containing:
            - phs_score: The p-value
            - is_significant: Boolean indicating if result is significant
            - alpha: The significance threshold used
            - interpretation: Text interpretation of the result
            - total_pairs: Number of leaf pairs analyzed
            - computation_time_ms: Computation time
            
    Example:
        >>> result = stellars.phs_significance_test(tree, alpha=0.01)
        >>> print(result['interpretation'])
        >>> if result['is_significant']:
        ...     print("Reject null hypothesis - significant homoplasies detected")
    """
    phs_result = phs_from_cassiopeia(tree, **kwargs)
    p_value = phs_result['phs_score']
    is_significant = p_value < alpha
    
    # Generate interpretation
    if p_value < 0.001:
        interpretation = f"Highly significant evidence of excess homoplasies (p = {p_value:.2e})"
    elif p_value < alpha:
        interpretation = f"Significant evidence of excess homoplasies (p = {p_value:.3f})"
    elif p_value < 0.1:
        interpretation = f"Marginally significant evidence of homoplasies (p = {p_value:.3f})"
    else:
        interpretation = f"No significant evidence of excess homoplasies (p = {p_value:.3f})"
    
    return {
        'phs_score': p_value,
        'is_significant': is_significant,
        'alpha': alpha,
        'interpretation': interpretation,
        'total_pairs': phs_result['total_pairs'],
        'computation_time_ms': phs_result['computation_time_ms'],
        'method_used': phs_result['method_used']
    }
