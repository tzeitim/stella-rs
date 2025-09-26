#!/usr/bin/env python3
"""
Core metrics computation functions extracted from reconstruction_worker.py.
This module provides the canonical implementation of all phylogenetic metrics
used in the cascade workflow, ensuring consistency across analysis and re-analysis.
"""

import logging
import numpy as np
import os
import time
from typing import Dict, Any, Optional, List, Tuple
import cassiopeia as cass
import stellars

logger = logging.getLogger(__name__)


def extract_parameters_from_tree(tree: cass.data.CassiopeiaTree,
                                source_name: str = "tree") -> Dict[str, Any]:
    """
    Extract mutation parameters from a Cassiopeia tree.
    Returns both simulation-derived and ground truth parameters if available.

    This implements the parameter extraction logic from reconstruction_worker.py lines 232-256.
    """
    params = {}

    # Extract from character matrix (simulation method)
    cm = tree.character_matrix.to_numpy()
    proportion_mutated = np.sum(cm > 0) / np.sum(cm >= 0)
    lam_from_matrix = -np.log(1.0 - proportion_mutated) if proportion_mutated < 1.0 else np.inf

    # Extract collision probability from priors or estimate uniformly
    if hasattr(tree, 'priors') and tree.priors:
        priors = np.array(list(tree.priors[0].values()))
        q_from_priors = np.sum(priors ** 2)
        q_source = "priors"
    else:
        n_states = cm[cm >= 0].max() if (cm >= 0).any() else 1
        q_from_priors = 1.0 / max(n_states, 1)
        q_source = "uniform_estimate"

    params.update({
        f'lam_from_{source_name}': lam_from_matrix,
        f'q_from_{source_name}': q_from_priors,
        f'proportion_mutated_{source_name}': proportion_mutated,
        f'q_{source_name}_source': q_source
    })

    # Extract stored ground truth parameters if available
    if hasattr(tree, 'parameters'):
        gt_params = tree.parameters
        params.update({
            'lam_gt': gt_params.get('lam_true', gt_params.get('lam_gt')),
            'q_gt': gt_params.get('q_true', gt_params.get('q_gt')),
            'proportion_mutated_gt': gt_params.get('proportion_mutated_gt')
        })

    return params


def calculate_phs_scores(tree: cass.data.CassiopeiaTree,
                        lam_sim: float, q_sim: float,
                        lam_gt: Optional[float] = None, q_gt: Optional[float] = None) -> Dict[str, float]:
    """
    Calculate PHS scores using both simulation and ground truth parameters.

    This implements the corrected PHS calculation from reconstruction_worker.py lines 348-417
    with proper character matrix mapping and internal state handling.
    """
    results = {}

    # Get tree structure and character data
    tree_newick = tree.get_newick(record_branch_lengths=True, record_node_names=True)

    # CRITICAL FIX: Extract leaf names for proper character matrix mapping (line 354)
    leaf_names = list(tree.character_matrix.index)
    logger.debug(f"Using {len(leaf_names)} leaf names for proper character matrix mapping")

    # Correct character matrix format (line 358)
    character_matrix = tree.character_matrix.values.astype(int).tolist()

    # Extract internal character states from Cassiopeia (lines 361-374)
    internal_character_states = {}
    try:
        for internal_node in tree.internal_nodes:
            node_name = str(internal_node)
            internal_states = tree.get_character_states(internal_node)
            if hasattr(internal_states, 'tolist'):
                internal_character_states[node_name] = internal_states.tolist()
            else:
                internal_character_states[node_name] = list(internal_states)
    except Exception as e:
        logger.warning(f"Failed to extract internal character states: {e}")
        internal_character_states = {}

    logger.debug(f"Extracted {len(internal_character_states)} internal node character states")

    # Calculate PHS with simulation parameters
    try:
        phs_result_sim = stellars.phs_optimized(
            tree_newick=tree_newick,
            character_matrix=character_matrix,
            internal_character_states=internal_character_states,
            mutation_rate=lam_sim,
            collision_probability=q_sim,
            missing_state=-1,
            unedited_state=0,
            use_provided_internal_states=True,
            leaf_names=leaf_names  # CRITICAL FIX: Proper character matrix mapping
        )
        results['phs_sim'] = phs_result_sim['phs_score']
    except Exception as e:
        logger.warning(f"PHS calculation with simulation params failed: {e}")
        results['phs_sim'] = np.nan

    # Calculate PHS with ground truth parameters if available
    if lam_gt is not None and q_gt is not None:
        try:
            phs_result_gt = stellars.phs_optimized(
                tree_newick=tree_newick,
                character_matrix=character_matrix,
                internal_character_states=internal_character_states,
                mutation_rate=lam_gt,
                collision_probability=q_gt,
                missing_state=-1,
                unedited_state=0,
                use_provided_internal_states=True,
                leaf_names=leaf_names
            )
            results['phs_gt'] = phs_result_gt['phs_score']
        except Exception as e:
            logger.warning(f"PHS calculation with GT params failed: {e}")
            results['phs_gt'] = np.nan
    else:
        results['phs_gt'] = np.nan

    return results


def calculate_tree_distances(reconstructed_tree: cass.data.CassiopeiaTree,
                           reference_tree: cass.data.CassiopeiaTree,
                           config: Optional[Dict[str, Any]] = None) -> Dict[str, float]:
    """
    Calculate RF distance and triplets correctness between two trees.

    This implements the distance calculations from reconstruction_worker.py lines 296-331.
    """
    results = {}

    # Robinson-Foulds distance (lines 323-331)
    try:
        rf_result = stellars.rf_distance(
            reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True),
            reference_tree.get_newick(record_branch_lengths=True, record_node_names=True)
        )
        results['RF_distance'] = rf_result.normalized_distance
    except Exception as e:
        logger.warning(f"RF distance calculation failed: {e}")
        results['RF_distance'] = np.nan

    # Triplets correctness (lines 296-320)
    try:
        # Get configuration parameters
        max_threads = os.environ.get('LSB_MAX_NUM_PROCESSORS', '16')
        triplets_trials = config.get('analysis', {}).get('triplets_trials', 1000) if config else 1000

        # Generate reproducible seed
        base_seed = config.get('analysis', {}).get('reconstruction_seed', 42) if config else 42
        triplets_seed = (base_seed + hash(str(reconstructed_tree))) % (2**31)

        triplets_result = stellars.triplets_correct_parallel(
            reference_tree.get_newick(record_branch_lengths=True, record_node_names=True),
            reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True),
            number_of_trials=triplets_trials,
            min_triplets_at_depth=1,
            seed=triplets_seed,
            max_threads=int(max_threads)
        )
        results['triplets_distance'] = 1.0 - triplets_result['all_triplets_correct'][0]
    except Exception as e:
        logger.warning(f"Triplets calculation failed: {e}")
        results['triplets_distance'] = np.nan

    return results


def calculate_parsimony_score(tree: cass.data.CassiopeiaTree) -> Dict[str, Any]:
    """
    Calculate parsimony score for a tree.

    This implements the parsimony calculation from reconstruction_worker.py lines 335-346.
    """
    try:
        parsimony_result = stellars.parsimony_score(
            tree_newick=tree.get_newick(record_branch_lengths=True, record_node_names=True),
            character_matrix=tree.character_matrix.to_numpy(),
            missing_state=-1,
            unedited_state=0
        )
        return parsimony_result
    except Exception as e:
        logger.warning(f"Parsimony calculation failed: {e}")
        return {'parsimony_score': np.nan}


def calculate_likelihood_score(tree: cass.data.CassiopeiaTree) -> Dict[str, Any]:
    """
    Calculate likelihood score for a tree using Cassiopeia's implementation.

    This implements the likelihood calculation from reconstruction_worker.py lines 419-463.
    """
    try:
        # Parameter validation - needs at least ONE of the required parameters (same logic as reconstruction_worker.py)
        required_params = ['heritable_missing_rate', 'stochastic_missing_rate', 'stochastic_missing_probability']
        has_required_param = False

        for param in required_params:
            if hasattr(tree, 'parameters') and param in tree.parameters:
                has_required_param = True
                break

        if not has_required_param:
            missing_params = [p for p in required_params if not hasattr(tree, 'parameters') or p not in tree.parameters]
            logger.warning(f"Missing required parameters for likelihood calculation: {missing_params}")
            logger.warning(f"Available parameters: {list(tree.parameters.keys()) if hasattr(tree, 'parameters') else 'None'}")
            raise ValueError(f"Tree missing required parameters: {missing_params}")

        # Calculate likelihood using Cassiopeia's robust implementation
        with np.errstate(divide='ignore'):
            likelihood_value = cass.tools.tree_metrics.calculate_likelihood_continuous(tree)

        # Validate result
        if np.isnan(likelihood_value) or np.isinf(likelihood_value):
            raise ValueError(f"Likelihood calculation returned invalid value: {likelihood_value}")

        return {'log_likelihood': -likelihood_value}

    except Exception as e:
        logger.warning(f"Likelihood calculation failed: {e}")
        return {'log_likelihood': np.nan}


def compute_all_metrics(reconstructed_tree: cass.data.CassiopeiaTree,
                       reference_tree: cass.data.CassiopeiaTree,
                       cas9_tree: Optional[cass.data.CassiopeiaTree] = None,
                       config: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Compute all phylogenetic metrics for a reconstructed tree.

    This is the main entry point that orchestrates all metric calculations,
    equivalent to the metrics computation in reconstruction_worker.py.

    Args:
        reconstructed_tree: The reconstructed phylogenetic tree
        reference_tree: The ground truth reference tree
        cas9_tree: Original CAS9 tree for parameter extraction (optional)
        config: Configuration dictionary (optional)

    Returns:
        Dictionary containing all computed metrics
    """
    results = {}

    # Extract parameters
    if cas9_tree is not None:
        params = extract_parameters_from_tree(cas9_tree, "simulation")
        results.update(params)

        # Use CAS9-derived parameters for computations
        lam_sim = params['lam_from_simulation']
        q_sim = params['q_from_simulation']
    else:
        # Fall back to reconstructed tree parameters
        params = extract_parameters_from_tree(reconstructed_tree, "reconstructed")
        results.update(params)
        lam_sim = params['lam_from_reconstructed']
        q_sim = params['q_from_reconstructed']

    # Get ground truth parameters
    lam_gt = results.get('lam_gt')
    q_gt = results.get('q_gt')

    # Calculate PHS scores
    phs_scores = calculate_phs_scores(reconstructed_tree, lam_sim, q_sim, lam_gt, q_gt)
    results.update(phs_scores)

    # Calculate tree distances
    distances = calculate_tree_distances(reconstructed_tree, reference_tree, config)
    results.update(distances)

    # Calculate parsimony score
    parsimony = calculate_parsimony_score(reconstructed_tree)
    results.update(parsimony)

    # Calculate likelihood score
    likelihood = calculate_likelihood_score(reconstructed_tree)
    results.update(likelihood)

    return results


def create_metrics_rows(base_metrics: Dict[str, Any],
                       reconstruction_id: str,
                       solver: str,
                       tier_info: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Create standardized metrics rows in pipeline format.

    Returns two rows: one with simulation parameters, one with GT parameters (if available).
    This maintains compatibility with the existing pipeline output format.
    """
    rows = []

    # Base row structure
    base_row = {
        'reconstruction_id': reconstruction_id,
        'solver': solver,
        **tier_info,
        **base_metrics
    }

    # Row 1: Simulation parameters
    row_sim = base_row.copy()
    row_sim.update({
        'cPHS': base_metrics.get('phs_sim', np.nan),
        'cPHS_gt': base_metrics.get('phs_gt', np.nan),
        'phs_lam_source': 'simulation',
        'phs_q_source': 'simulation'
    })
    rows.append(row_sim)

    # Row 2: Ground truth parameters (if available)
    if base_metrics.get('phs_gt') is not None and not np.isnan(base_metrics.get('phs_gt', np.nan)):
        row_gt = base_row.copy()
        row_gt.update({
            'cPHS': base_metrics.get('phs_sim', np.nan),  # Keep simulation as primary
            'cPHS_gt': base_metrics.get('phs_gt', np.nan),
            'phs_lam_source': 'ground_truth',
            'phs_q_source': 'ground_truth'
        })
        rows.append(row_gt)

    return rows