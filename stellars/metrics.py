"""
Phylogenetic tree metrics with improved interface.

This module provides clean, Pythonic wrappers for parsimony, likelihood, 
and PHS calculations with standardized return types.
"""

from typing import Optional, Dict, List, Any, Union
from .results import PHSResult, ParsimonyResult, LikelihoodResult
from . import stellars


def phs(tree_newick: str,
        character_matrix: List[List[int]],
        internal_states: Optional[Dict[str, List[int]]] = None,
        mutation_rate: Optional[float] = None,
        collision_probability: Optional[float] = None,
        missing_state: int = -1,
        unedited_state: int = 0,
        max_threads: Optional[int] = None,
        use_internal_states: bool = True,
        leaf_names: Optional[List[str]] = None) -> PHSResult:
    """
    Calculate corrected Pairwise Homoplasy Score (cPHS) for a tree.
    
    Args:
        tree_newick: Tree in Newick format
        character_matrix: Character matrix as list of lists
        internal_states: Internal node character states (optional)
        mutation_rate: Mutation rate parameter λ (estimated if None)
        collision_probability: Collision probability parameter q (estimated if None)
        missing_state: Value representing missing data (default: -1)
        unedited_state: Value representing unedited state (default: 0)
        max_threads: Maximum threads for parallel computation
        use_internal_states: Whether to use provided internal states
        leaf_names: Names of leaves for proper character matrix mapping
        
    Returns:
        PHSResult with comprehensive PHS information
        
    Examples:
        >>> result = stellars.phs(tree_newick, char_matrix)
        >>> print(f"cPHS p-value: {result.phs_score:.2e}")
        >>> 
        >>> # With significance assessment
        >>> result = stellars.phs(tree_newick, char_matrix).assess_significance(0.01)
        >>> print(result.interpretation)
    """
    raw_result = stellars.phs_optimized(
        tree_newick=tree_newick,
        character_matrix=character_matrix,
        internal_character_states=internal_states or {},
        mutation_rate=mutation_rate,
        collision_probability=collision_probability,
        missing_state=missing_state,
        unedited_state=unedited_state,
        max_threads=max_threads,
        use_provided_internal_states=use_internal_states,
        leaf_names=leaf_names
    )
    
    return PHSResult(
        phs_score=raw_result['phs_score'],
        total_pairs=raw_result.get('total_pairs', 0),
        computation_time_ms=raw_result.get('computation_time_ms', 0.0),
        method_used=raw_result.get('method_used', 'phs_optimized'),
        parallel_chunks_used=raw_result.get('parallel_chunks_used')
    )


def phs_simple(tree_newick: str, character_matrix: List[List[int]], 
               warn_on_auto_params: bool = True, **kwargs) -> float:
    """
    Simple PHS calculation that returns just the p-value.
    
    CAUTION: This function may auto-estimate mutation_rate and collision_probability
    if not provided, which can lead to inconsistent results. For production use,
    consider providing explicit parameters or using phs() for full control.
    
    Args:
        tree_newick: Tree in Newick format
        character_matrix: Character matrix as list of lists
        warn_on_auto_params: If True, warn when parameters are auto-estimated
        **kwargs: Additional arguments passed to phs()
        
    Returns:
        float: The cPHS p-value
        
    Examples:
        >>> # Basic usage (may auto-estimate parameters)
        >>> p_value = stellars.phs_simple(tree_newick, char_matrix)
        >>> if p_value < 0.05:
        ...     print(f"Significant homoplasies detected (p = {p_value:.2e})")
        
        >>> # Recommended: provide explicit parameters
        >>> p_value = stellars.phs_simple(tree_newick, char_matrix,
        ...                               mutation_rate=0.3, collision_probability=0.1)
    """
    # Check if critical parameters are missing
    if warn_on_auto_params:
        missing_params = []
        if 'mutation_rate' not in kwargs or kwargs['mutation_rate'] is None:
            missing_params.append('mutation_rate')
        if 'collision_probability' not in kwargs or kwargs['collision_probability'] is None:
            missing_params.append('collision_probability')
        
        if missing_params:
            import warnings
            warnings.warn(f"Parameters {missing_params} not provided to phs_simple(). "
                         f"Auto-estimation may lead to inconsistent results. "
                         f"Consider providing explicit values or using phs() for full control.")
    
    try:
        result = phs(tree_newick, character_matrix, **kwargs)
        return result.phs_score
    except Exception as e:
        raise RuntimeError(f"PHS calculation failed. This may be due to auto-parameter estimation issues. "
                          f"Try providing explicit mutation_rate and collision_probability. Original error: {e}")


def phs_significance_test(tree_newick: str, 
                         character_matrix: List[List[int]],
                         alpha: float = 0.05,
                         **kwargs) -> PHSResult:
    """
    Perform significance test for homoplasies with interpretation.
    
    Args:
        tree_newick: Tree in Newick format
        character_matrix: Character matrix as list of lists
        alpha: Significance threshold (default: 0.05)
        **kwargs: Additional arguments passed to phs()
        
    Returns:
        PHSResult with significance assessment and interpretation
        
    Examples:
        >>> result = stellars.phs_significance_test(tree_newick, char_matrix, alpha=0.01)
        >>> print(result.interpretation)
        >>> if result.is_significant:
        ...     print("Significant homoplasies detected!")
    """
    result = phs(tree_newick, character_matrix, **kwargs)
    return result.assess_significance(alpha)


def parsimony(tree_newick: str,
              character_matrix: List[List[int]],
              internal_states: Optional[Dict[str, List[int]]] = None,
              missing_state: int = -1,
              unedited_state: int = 0,
              calculate_p_value: bool = False,
              mutation_rate: Optional[float] = None) -> ParsimonyResult:
    """
    Calculate parsimony score for a tree.
    
    Args:
        tree_newick: Tree in Newick format
        character_matrix: Character matrix as list of lists
        internal_states: Internal node character states (optional)
        missing_state: Value representing missing data
        unedited_state: Value representing unedited state
        calculate_p_value: Whether to also calculate statistical p-value
        mutation_rate: Mutation rate for p-value calculation (if needed)
        
    Returns:
        ParsimonyResult with parsimony information
        
    Examples:
        >>> result = stellars.parsimony(tree_newick, char_matrix)
        >>> print(f"Parsimony score: {result.score}")
        >>> 
        >>> # With p-value calculation
        >>> result = stellars.parsimony(tree_newick, char_matrix, 
        ...                           calculate_p_value=True, mutation_rate=0.3)
        >>> print(f"Parsimony p-value: {result.p_value:.3f}")
    """
    raw_result = stellars.parsimony_score(
        tree_newick=tree_newick,
        character_matrix=character_matrix,
        internal_character_states=internal_states,
        missing_state=missing_state,
        unedited_state=unedited_state
    )
    
    p_value = None
    if calculate_p_value:
        p_result = stellars.parsimony_p_value(
            tree_newick=tree_newick,
            character_matrix=character_matrix,
            internal_character_states=internal_states,
            mutation_rate=mutation_rate,
            missing_state=missing_state,
            unedited_state=unedited_state
        )
        p_value = p_result['p_value']
    
    return ParsimonyResult(
        score=raw_result['parsimony_score'],
        total_mutations=raw_result.get('total_mutations', 0),
        p_value=p_value,
        computation_time_ms=raw_result.get('computation_time_ms', 0.0),
        method_used=raw_result.get('method_used', 'parsimony'),
        internal_states_inferred=raw_result.get('internal_states_inferred', False)
    )


def parsimony_simple(tree_newick: str, character_matrix: List[List[int]], **kwargs) -> int:
    """
    Simple parsimony calculation that returns just the score.
    
    Args:
        tree_newick: Tree in Newick format
        character_matrix: Character matrix as list of lists
        **kwargs: Additional arguments passed to parsimony()
        
    Returns:
        int: The parsimony score
        
    Examples:
        >>> score = stellars.parsimony_simple(tree_newick, char_matrix)
        >>> print(f"Parsimony score: {score}")
    """
    result = parsimony(tree_newick, character_matrix, **kwargs)
    return result.score


def likelihood(tree_newick: str,
               character_matrix: List[List[int]],
               internal_states: Optional[Dict[str, List[int]]] = None,
               mutation_rate: Optional[float] = None,
               collision_probability: Optional[float] = None,
               missing_state: int = -1,
               unedited_state: int = 0) -> LikelihoodResult:
    """
    Calculate likelihood for a tree.
    
    Args:
        tree_newick: Tree in Newick format
        character_matrix: Character matrix as list of lists
        internal_states: Internal node character states (optional)
        mutation_rate: Mutation rate parameter (estimated if None)
        collision_probability: Collision probability parameter (estimated if None)
        missing_state: Value representing missing data
        unedited_state: Value representing unedited state
        
    Returns:
        LikelihoodResult with likelihood information
        
    Examples:
        >>> result = stellars.likelihood(tree_newick, char_matrix)
        >>> print(f"Log-likelihood: {result.log_likelihood:.2f}")
        >>> print(f"Likelihood: {result.likelihood:.2e}")
    """
    raw_result = stellars.likelihood_score(
        tree_newick=tree_newick,
        character_matrix=character_matrix,
        internal_character_states=internal_states,
        mutation_rate=mutation_rate,
        collision_probability=collision_probability,
        missing_state=missing_state,
        unedited_state=unedited_state
    )
    
    return LikelihoodResult(
        log_likelihood=raw_result['log_likelihood'],
        likelihood=raw_result['likelihood'],
        n_characters=raw_result.get('n_characters', 0),
        n_leaves=raw_result.get('n_leaves', 0),
        computation_time_ms=raw_result.get('computation_time_ms', 0.0),
        method_used=raw_result.get('method_used', 'likelihood')
    )


def likelihood_simple(tree_newick: str, character_matrix: List[List[int]], 
                     return_log: bool = True, **kwargs) -> float:
    """
    Simple likelihood calculation that returns just the likelihood value.
    
    Args:
        tree_newick: Tree in Newick format
        character_matrix: Character matrix as list of lists
        return_log: If True, return log-likelihood; if False, return raw likelihood
        **kwargs: Additional arguments passed to likelihood()
        
    Returns:
        float: The likelihood or log-likelihood value
        
    Examples:
        >>> log_likelihood = stellars.likelihood_simple(tree_newick, char_matrix)
        >>> likelihood = stellars.likelihood_simple(tree_newick, char_matrix, return_log=False)
    """
    result = likelihood(tree_newick, character_matrix, **kwargs)
    return result.log_likelihood if return_log else result.likelihood


def likelihood_distance(reconstructed_tree_newick: str,
                       ground_truth_tree_newick: str,
                       character_matrix: List[List[int]],
                       reconstructed_states: Optional[Dict[str, List[int]]] = None,
                       ground_truth_states: Optional[Dict[str, List[int]]] = None,
                       mutation_rate: Optional[float] = None,
                       collision_probability: Optional[float] = None,
                       missing_state: int = -1,
                       unedited_state: int = 0) -> float:
    """
    Calculate likelihood distance between reconstructed and ground truth trees.
    
    Args:
        reconstructed_tree_newick: Reconstructed tree in Newick format
        ground_truth_tree_newick: Ground truth tree in Newick format
        character_matrix: Character matrix as list of lists
        reconstructed_states: Internal states of reconstructed tree
        ground_truth_states: Internal states of ground truth tree
        mutation_rate: Mutation rate parameter
        collision_probability: Collision probability parameter
        missing_state: Value representing missing data
        unedited_state: Value representing unedited state
        
    Returns:
        float: The likelihood distance
        
    Examples:
        >>> distance = stellars.likelihood_distance(recon_tree, gt_tree, char_matrix)
        >>> print(f"Likelihood distance: {distance:.3f}")
    """
    raw_result = stellars.likelihood_distance(
        reconstructed_tree_newick=reconstructed_tree_newick,
        ground_truth_tree_newick=ground_truth_tree_newick,
        character_matrix=character_matrix,
        internal_states_reconstructed=reconstructed_states,
        internal_states_ground_truth=ground_truth_states,
        mutation_rate=mutation_rate,
        collision_probability=collision_probability,
        missing_state=missing_state,
        unedited_state=unedited_state
    )
    
    return raw_result['likelihood_distance']