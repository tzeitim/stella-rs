"""
Distance calculation functions with improved, standardized interface.

This module provides clean, Pythonic wrappers around the core stellars distance
functions with consistent return types and naming conventions.
"""

from typing import Union, Optional
from .results import RFDistanceResult, TripletsResult
from . import stellars


def rf_distance(tree1_newick: str, tree2_newick: str, normalize: bool = True) -> Union[RFDistanceResult, float]:
    """
    Calculate Robinson-Foulds distance between two trees.
    
    Args:
        tree1_newick: First tree in Newick format
        tree2_newick: Second tree in Newick format  
        normalize: If True, return normalized distance. If False, return raw distance.
        
    Returns:
        RFDistanceResult with comprehensive distance information
        
    Examples:
        >>> # Get normalized RF distance (most common use case)
        >>> distance = stellars.rf_distance(tree1, tree2)
        >>> print(f"Normalized RF distance: {distance.normalized}")
        
        >>> # Get raw RF distance  
        >>> result = stellars.rf_distance(tree1, tree2, normalize=False)
        >>> print(f"Raw RF distance: {result.distance}")
    """
    if normalize:
        raw_result = stellars.robinson_foulds_distance_normalized(tree1_newick, tree2_newick)
        return RFDistanceResult(
            distance=raw_result['rf_distance'],
            normalized_distance=raw_result['normalized_rf_distance'],
            max_possible_distance=raw_result.get('max_possible_rf'),
            computation_time_ms=raw_result.get('computation_time_ms', 0.0),
            method_used=raw_result.get('method_used', 'robinson_foulds_normalized'),
            tree1_leaves=raw_result.get('tree1_leaves'),
            tree2_leaves=raw_result.get('tree2_leaves')
        )
    else:
        raw_result = stellars.robinson_foulds_distance(tree1_newick, tree2_newick)
        return RFDistanceResult(
            distance=raw_result['rf_distance'],
            computation_time_ms=raw_result.get('computation_time_ms', 0.0),
            method_used=raw_result.get('method_used', 'robinson_foulds')
        )


def rf_distance_simple(tree1_newick: str, tree2_newick: str, normalize: bool = True) -> float:
    """
    Simple RF distance calculation that returns just the distance value.
    
    Args:
        tree1_newick: First tree in Newick format
        tree2_newick: Second tree in Newick format
        normalize: If True, return normalized distance
        
    Returns:
        float: The RF distance value
        
    Examples:
        >>> distance = stellars.rf_distance_simple(tree1, tree2)
        >>> print(f"RF distance: {distance:.3f}")
    """
    result = rf_distance(tree1_newick, tree2_newick, normalize=normalize)
    return result.normalized if normalize else result.distance


def triplets_distance(tree1_newick: str, tree2_newick: str, 
                     method: str = 'optimized',
                     trials: int = 10000,
                     min_depth: int = 1,
                     seed: Optional[int] = None,
                     max_threads: Optional[int] = None) -> TripletsResult:
    """
    Calculate triplets distance between two trees using the specified method.
    
    Args:
        tree1_newick: First tree in Newick format
        tree2_newick: Second tree in Newick format
        method: Method to use ('basic', 'fast', 'ultra', 'parallel', 'optimized')
        trials: Number of trials for sampling-based methods
        min_depth: Minimum depth for triplets consideration
        seed: Random seed for reproducibility
        max_threads: Maximum threads for parallel methods
        
    Returns:
        TripletsResult with comprehensive triplets information
        
    Examples:
        >>> # Most common usage - get triplets distance
        >>> result = stellars.triplets_distance(tree1, tree2)
        >>> print(f"Triplets distance: {result.distance:.3f}")
        >>> print(f"Triplets correctness: {result.correctness:.3f}")
        
        >>> # Use specific method
        >>> result = stellars.triplets_distance(tree1, tree2, method='parallel', max_threads=8)
    """
    # Select the appropriate function based on method
    method_map = {
        'basic': stellars.triplets_correct,
        'fast': stellars.triplets_correct_fast,
        'ultra': stellars.triplets_correct_ultra,
        'parallel': stellars.triplets_correct_parallel,
        'optimized': stellars.triplets_correct_optimized
    }
    
    if method not in method_map:
        raise ValueError(f"Unknown method '{method}'. Choose from: {list(method_map.keys())}")
    
    func = method_map[method]
    
    # Call the appropriate function with relevant parameters
    if method in ['parallel', 'optimized']:
        raw_result = func(tree1_newick, tree2_newick, trials, min_depth, seed, max_threads)
    else:
        raw_result = func(tree1_newick, tree2_newick, trials, min_depth, seed)
    
    return TripletsResult.from_raw_result(raw_result)


def triplets_distance_simple(tree1_newick: str, tree2_newick: str, **kwargs) -> float:
    """
    Simple triplets distance calculation that returns just the distance value.
    
    Args:
        tree1_newick: First tree in Newick format  
        tree2_newick: Second tree in Newick format
        **kwargs: Additional arguments passed to triplets_distance()
        
    Returns:
        float: The triplets distance value (1 - correctness)
        
    Examples:
        >>> distance = stellars.triplets_distance_simple(tree1, tree2)
        >>> print(f"Triplets distance: {distance:.3f}")
    """
    result = triplets_distance(tree1_newick, tree2_newick, **kwargs)
    return result.distance


def triplets_correctness(tree1_newick: str, tree2_newick: str, **kwargs) -> float:
    """
    Calculate triplets correctness between two trees.
    
    Args:
        tree1_newick: First tree in Newick format
        tree2_newick: Second tree in Newick format  
        **kwargs: Additional arguments passed to triplets_distance()
        
    Returns:
        float: The proportion of triplets that are correct
        
    Examples:
        >>> correctness = stellars.triplets_correctness(tree1, tree2)
        >>> print(f"Triplets correctness: {correctness:.3f}")
    """
    result = triplets_distance(tree1_newick, tree2_newick, **kwargs)
    return result.correctness


# Convenient aliases for backward compatibility and discoverability
robinson_foulds_distance = rf_distance
robinson_foulds_simple = rf_distance_simple
triplets_correct_distance = triplets_distance