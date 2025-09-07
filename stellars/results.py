"""
Structured result classes for stellars computations.

This module provides dataclasses that standardize the return types of stellars
functions, making the API more Pythonic and user-friendly.
"""

from dataclasses import dataclass
from typing import Optional, Dict, Any, Union


@dataclass
class RFDistanceResult:
    """Result of Robinson-Foulds distance calculation."""
    
    distance: float
    """Raw RF distance between trees"""
    
    normalized_distance: Optional[float] = None  
    """RF distance normalized by maximum possible RF distance"""
    
    max_possible_distance: Optional[int] = None
    """Maximum possible RF distance for trees with this number of leaves"""
    
    computation_time_ms: float = 0.0
    """Time taken for computation in milliseconds"""
    
    method_used: str = "robinson_foulds"
    """Method identifier used for computation"""
    
    tree1_leaves: Optional[int] = None
    """Number of leaves in first tree"""
    
    tree2_leaves: Optional[int] = None  
    """Number of leaves in second tree"""

    @property
    def normalized(self) -> float:
        """Convenience property for normalized distance."""
        return self.normalized_distance if self.normalized_distance is not None else self.distance


@dataclass  
class TripletsResult:
    """Result of triplets correctness calculation."""
    
    correctness: float
    """Proportion of triplets that are correct"""
    
    distance: float
    """Triplets distance (1 - correctness)"""
    
    all_correct: Dict[int, float]  
    """All triplets correctness by depth"""
    
    resolvable_correct: Dict[int, float]
    """Resolvable triplets correctness by depth"""
    
    unresolved_correct: Dict[int, float] 
    """Unresolved triplets correctness by depth"""
    
    proportion_unresolvable: Dict[int, float]
    """Proportion of unresolvable triplets by depth"""
    
    total_triplets_checked: Optional[int] = None
    """Total number of triplets evaluated"""
    
    method_used: str = "triplets_correct"
    """Method identifier used for computation"""
    
    computation_time_ms: float = 0.0
    """Time taken for computation in milliseconds"""
    
    parallel_chunks_used: Optional[int] = None
    """Number of parallel chunks used (if applicable)"""

    @classmethod
    def from_raw_result(cls, raw_result: Dict[str, Any]) -> 'TripletsResult':
        """Create TripletsResult from raw stellars output."""
        # Extract the primary correctness value (use depth 0 if available)
        all_correct = raw_result.get('all_triplets_correct', {})
        primary_correctness = all_correct.get(0, 0.0) if all_correct else 0.0
        
        return cls(
            correctness=primary_correctness,
            distance=1.0 - primary_correctness,
            all_correct=all_correct,
            resolvable_correct=raw_result.get('resolvable_triplets_correct', {}),
            unresolved_correct=raw_result.get('unresolved_triplets_correct', {}), 
            proportion_unresolvable=raw_result.get('proportion_unresolvable', {}),
            total_triplets_checked=raw_result.get('total_triplets_checked'),
            method_used=raw_result.get('method_used', 'triplets_correct'),
            parallel_chunks_used=raw_result.get('parallel_chunks_used')
        )


@dataclass
class PHSResult:
    """Result of Pairwise Homoplasy Score calculation."""
    
    phs_score: float
    """The corrected pairwise homoplasy p-value"""
    
    is_significant: Optional[bool] = None
    """Whether the result is statistically significant (if alpha was provided)"""
    
    alpha: Optional[float] = None
    """Significance threshold used (if applicable)"""
    
    total_pairs: int = 0
    """Number of leaf pairs analyzed"""
    
    computation_time_ms: float = 0.0
    """Time taken for computation in milliseconds"""
    
    method_used: str = "phs_optimized"
    """Method identifier used for computation"""
    
    parallel_chunks_used: Optional[int] = None
    """Number of parallel chunks used (if applicable)"""
    
    interpretation: Optional[str] = None
    """Human-readable interpretation of the result"""

    def assess_significance(self, alpha: float = 0.05) -> 'PHSResult':
        """Add significance assessment to the result."""
        is_significant = self.phs_score < alpha
        
        # Generate interpretation
        if self.phs_score < 0.001:
            interpretation = f"Highly significant evidence of excess homoplasies (p = {self.phs_score:.2e})"
        elif self.phs_score < alpha:
            interpretation = f"Significant evidence of excess homoplasies (p = {self.phs_score:.3f})"
        elif self.phs_score < 0.1:
            interpretation = f"Marginally significant evidence of homoplasies (p = {self.phs_score:.3f})"
        else:
            interpretation = f"No significant evidence of excess homoplasies (p = {self.phs_score:.3f})"
        
        # Return a copy with significance info
        return PHSResult(
            phs_score=self.phs_score,
            is_significant=is_significant,
            alpha=alpha,
            total_pairs=self.total_pairs,
            computation_time_ms=self.computation_time_ms,
            method_used=self.method_used,
            parallel_chunks_used=self.parallel_chunks_used,
            interpretation=interpretation
        )


@dataclass
class ParsimonyResult:
    """Result of parsimony score calculation."""
    
    score: int
    """The parsimony score (number of mutations)"""
    
    total_mutations: int = 0
    """Total number of mutations counted"""
    
    p_value: Optional[float] = None
    """Statistical p-value (if calculated)"""
    
    computation_time_ms: float = 0.0
    """Time taken for computation in milliseconds"""
    
    method_used: str = "parsimony"
    """Method identifier used for computation"""
    
    internal_states_inferred: bool = False
    """Whether internal states were inferred or provided"""


@dataclass
class LikelihoodResult:
    """Result of likelihood calculation."""
    
    log_likelihood: float
    """Log-likelihood value"""
    
    likelihood: float
    """Raw likelihood value (exp(log_likelihood))"""
    
    distance: Optional[float] = None
    """Likelihood distance (if comparing to ground truth)"""
    
    n_characters: int = 0
    """Number of characters analyzed"""
    
    n_leaves: int = 0
    """Number of leaves in the tree"""
    
    computation_time_ms: float = 0.0
    """Time taken for computation in milliseconds"""
    
    method_used: str = "likelihood"
    """Method identifier used for computation"""


# Union type for all result types
AnyResult = Union[RFDistanceResult, TripletsResult, PHSResult, ParsimonyResult, LikelihoodResult]