#!/usr/bin/env python3
"""
Shared metrics calculation functions used by both reconstruction pipeline and re-analysis.
Extracted from reconstruction_worker.py to avoid duplication.
"""

import logging
import numpy as np
import os
import cassiopeia as cass
import stellars

logger = logging.getLogger(__name__)


def calculate_metrics_for_trees(reconstructed_tree: cass.data.CassiopeiaTree,
                               reference_tree: cass.data.CassiopeiaTree,
                               config: dict,
                               gt_instance_id: int = 0,
                               cas9_simulation_id: int = 0,
                               reconstruction_id: int = 0,
                               solver_name: str = "unknown",
                               lam_sim: float = None,
                               q_sim: float = None,
                               lam_gt: float = None,
                               q_gt: float = None,
                               tier_num: int = None) -> dict:
    """
    Calculate all metrics between reconstructed and reference trees.
    This is extracted from reconstruction_worker.py to be shared between pipeline and re-analysis.
    """
    metrics = {}

    # Copy parameters from reference to reconstructed tree if needed
    if hasattr(reference_tree, 'parameters'):
        if not hasattr(reconstructed_tree, 'parameters'):
            reconstructed_tree.parameters = {}

        # Copy essential parameters for likelihood calculation (from reconstruction_worker.py lines 195-210)
        required_params = [
            'heritable_missing_rate', 'stochastic_missing_rate', 'stochastic_missing_probability',
            'mutation_rate', 'mutation_rates', 'state_priors', 'lam_gt', 'q_gt',
            'proportion_mutated_gt', 'cas9_simulation_parameters'
        ]

        for param_name in required_params:
            if param_name in reference_tree.parameters:
                reconstructed_tree.parameters[param_name] = reference_tree.parameters[param_name]

    # PHS calculations (from reconstruction_worker.py lines 350-415)
    try:
        max_threads = os.environ.get('LSB_MAX_NUM_PROCESSORS', '16')
        logger.info(f"Starting PHS calculations with {max_threads} threads")

        # Prepare data for PHS calculation (same as reconstruction_worker.py)
        tree_newick = reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True)
        leaf_names = list(reconstructed_tree.character_matrix.index)
        character_matrix = reconstructed_tree.character_matrix.values.astype(int).tolist()
        logger.info(f"Prepared character matrix: {len(leaf_names)} leaves, {len(character_matrix[0]) if character_matrix else 0} characters")

        # Extract internal character states
        internal_character_states = {}
        for internal_node in reconstructed_tree.internal_nodes:
            node_name = str(internal_node)
            try:
                internal_states = reconstructed_tree.get_character_states(internal_node)
                if hasattr(internal_states, 'tolist'):
                    internal_character_states[node_name] = internal_states.tolist()
                else:
                    internal_character_states[node_name] = list(internal_states)
            except Exception:
                internal_character_states = {}
                break
        logger.info(f"Extracted {len(internal_character_states)} internal character states")

        # Simulation parameters PHS
        if lam_sim is not None and q_sim is not None:
            logger.info(f"Computing cPHS with simulation parameters: lam={lam_sim:.4f}, q={q_sim:.6f}")
            phs_result_sim = stellars.phs_optimized(
                tree_newick=tree_newick,
                character_matrix=character_matrix,
                internal_character_states=internal_character_states,
                mutation_rate=lam_sim,
                collision_probability=q_sim,
                leaf_names=leaf_names,
                max_threads=int(max_threads)
            )
            metrics['cPHS_simulation'] = phs_result_sim['phs_score']
            logger.info(f"cPHS simulation result: {metrics['cPHS_simulation']:.6e}")
        else:
            logger.info("Skipping cPHS simulation (missing parameters)")
            metrics['cPHS_simulation'] = np.nan

        # Ground truth parameters PHS
        if lam_gt is not None and q_gt is not None:
            logger.info(f"Computing cPHS with ground truth parameters: lam={lam_gt:.4f}, q={q_gt:.6f}")
            phs_gt_result = stellars.phs_optimized(
                tree_newick=tree_newick,
                character_matrix=character_matrix,
                internal_character_states=internal_character_states,
                mutation_rate=lam_gt,
                collision_probability=q_gt,
                leaf_names=leaf_names,
                max_threads=int(max_threads)
            )
            metrics['cPHS_gt'] = phs_gt_result['phs_score']
            logger.info(f"cPHS ground truth result: {metrics['cPHS_gt']:.6e}")
        else:
            logger.info("Skipping cPHS ground truth (missing parameters)")
            metrics['cPHS_gt'] = np.nan

    except Exception as e:
        logger.warning(f"PHS calculation failed: {e}")
        metrics['cPHS_simulation'] = np.nan
        metrics['cPHS_gt'] = np.nan

    # Triplets distance (from reconstruction_worker.py lines 295-320)
    try:
        triplets_trials = config.get('analysis', {}).get('triplets_trials', 1000)
        base_seed = config.get('analysis', {}).get('reconstruction_seed',
                              config.get('execution', {}).get('random_seed', 42))
        triplets_seed = (base_seed + hash((gt_instance_id, cas9_simulation_id, reconstruction_id, solver_name))) % (2**31)
        logger.info(f"Computing triplets distance with {triplets_trials} trials, seed={triplets_seed}")

        triplets_result = stellars.triplets_correct_parallel(
            reference_tree.get_newick(record_branch_lengths=True, record_node_names=True),
            reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True),
            number_of_trials=triplets_trials,
            min_triplets_at_depth=1,
            seed=triplets_seed,
            max_threads=int(max_threads)
        )
        metrics['triplets_distance'] = 1.0 - triplets_result['all_triplets_correct'][0]
        logger.info(f"Triplets distance result: {metrics['triplets_distance']:.4f}")
    except Exception as e:
        logger.warning(f"Triplets calculation failed: {e}")
        metrics['triplets_distance'] = np.nan

    # Robinson-Foulds distance (from reconstruction_worker.py lines 322-331)
    try:
        logger.info("Computing Robinson-Foulds distance")
        rf_result = stellars.rf_distance(
            reference_tree.get_newick(record_branch_lengths=True, record_node_names=True),
            reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True)
        )
        metrics['RF_distance'] = rf_result.normalized_distance
        logger.info(f"RF distance result: {metrics['RF_distance']:.4f}")
    except Exception as e:
        logger.warning(f"RF distance calculation failed: {e}")
        metrics['RF_distance'] = np.nan

    # Parsimony score (from reconstruction_worker.py lines 333-345)
    try:
        logger.info("Computing parsimony score")
        parsimony_result = stellars.parsimony_score(
            tree_newick=reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True),
            character_matrix=reconstructed_tree.character_matrix.to_numpy(),
            missing_state=-1,
            unedited_state=0
        )
        metrics['parsimony_score'] = parsimony_result
        logger.info(f"Parsimony score computed: {parsimony_result.get('parsimony_score', 'N/A')}")
    except Exception as e:
        logger.warning(f"Parsimony calculation failed: {e}")
        metrics['parsimony_score'] = {'parsimony_score': np.nan}

    # Likelihood calculation following simulation_phs.py golden standard
    try:
        logger.info("Computing likelihood score")
        # Parameter validation - needs at least ONE of the required parameters
        required_params = ['heritable_missing_rate', 'stochastic_missing_rate', 'stochastic_missing_probability']
        has_required_param = False

        for param in required_params:
            if hasattr(reconstructed_tree, 'parameters') and param in reconstructed_tree.parameters:
                has_required_param = True
                break

        if not has_required_param:
            missing_params = [p for p in required_params if not hasattr(reconstructed_tree, 'parameters') or p not in reconstructed_tree.parameters]
            logger.warning(f"Missing required parameters for likelihood calculation: {missing_params}")
            logger.warning(f"Available parameters: {list(reconstructed_tree.parameters.keys()) if hasattr(reconstructed_tree, 'parameters') else 'None'}")
            raise ValueError(f"Tree missing required parameters: {missing_params}")

        # Following simulation_phs.py approach: use only k cassette characters for likelihood
        # Extract k from the tier configuration to match simulation_phs.py line 155
        k = None

        # Get k from tier configuration
        if tier_num is not None and config and 'cas9_tiers' in config:
            tier_config = config['cas9_tiers'].get(tier_num)
            if tier_config and 'k' in tier_config:
                k = tier_config['k']
                logger.info(f"Using k={k} cassettes from tier {tier_num} config")

        # Fallback: try to determine k from priors
        if k is None and hasattr(reconstructed_tree, 'priors') and reconstructed_tree.priors:
            # If priors exist, k should be the number of cassette characters that have priors
            # Following simulation_phs.py: gt_tree.priors = {i: priors for i in range(k)}
            k = len(reconstructed_tree.priors)
            logger.info(f"Detected k={k} cassettes from priors")

        # Final fallback: estimate from character matrix size
        if k is None:
            k = min(100, reconstructed_tree.character_matrix.shape[1])
            logger.warning(f"No tier config or priors found, using fallback k={k}")

        # Create likelihood tree with only k cassette characters (simulation_phs.py approach)
        cassette_char_matrix = reconstructed_tree.character_matrix.iloc[:, :k]

        likelihood_tree = cass.data.CassiopeiaTree(
            character_matrix=cassette_char_matrix,
            tree=reconstructed_tree.get_newick(),
            missing_state_indicator=reconstructed_tree.missing_state_indicator if hasattr(reconstructed_tree, 'missing_state_indicator') else -1
        )

        # Copy parameters
        for param, value in reconstructed_tree.parameters.items():
            likelihood_tree.parameters[param] = value

        # Set priors following simulation_phs.py approach: only for k cassettes
        if hasattr(reconstructed_tree, 'priors') and reconstructed_tree.priors:
            likelihood_tree.priors = {i: reconstructed_tree.priors[i] for i in range(k) if i in reconstructed_tree.priors}
            logger.info(f"Set priors for {len(likelihood_tree.priors)} cassettes")

        logger.info(f"Created likelihood tree with {k} cassettes (following simulation_phs.py)")

        with np.errstate(divide='ignore'):
            likelihood_value = cass.tools.tree_metrics.calculate_likelihood_continuous(likelihood_tree)

        # Validate likelihood result
        if np.isnan(likelihood_value) or np.isinf(likelihood_value):
            logger.warning(f"Invalid likelihood value: {likelihood_value}")
            raise ValueError(f"Likelihood calculation returned invalid value: {likelihood_value}")

        likelihood_result = {'log_likelihood': -likelihood_value}
        metrics['likelihood_score_simulation'] = likelihood_result
        metrics['likelihood_score_gt'] = likelihood_result
        metrics['likelihood_score'] = likelihood_result
        logger.info(f"Likelihood result: {-likelihood_value:.2f}")

    except Exception as e:
        logger.warning(f"Likelihood calculation failed: {e}")
        metrics['likelihood_score'] = np.nan
        metrics['likelihood_score_simulation'] = np.nan
        metrics['likelihood_score_gt'] = np.nan

    return metrics