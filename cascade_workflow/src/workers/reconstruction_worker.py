#!/usr/bin/env python3
"""
Reconstruction Worker - Level 3 of LSF Cascading Job System

Performs phylogenetic reconstruction using a specific solver on a Cas9 instance
and calculates all metrics. This is the final level that produces the actual results.
"""

import argparse
import json
import logging
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Dict, Any
import numpy as np
import polars as pl

# Fix scipy.errstate compatibility issue for spectral solver
import scipy
if not hasattr(scipy, 'errstate'):
    scipy.errstate = np.errstate

import cassiopeia as cass
from copy import deepcopy
import stellars
from convexml import convexml
from solver_config import get_solver_class, SOLVERS
import yaml
import sys

# Add utils directory to path - handle both original location and copied location
current_dir = Path(__file__).parent
utils_candidates = [
    current_dir.parent / "cascade_workflow" / "src" / "utils",  # Original structure
    current_dir.parent / "utils",  # Organized structure
    current_dir / "utils",  # If utils is copied alongside
    current_dir,  # If utils files are copied to same directory
]

for utils_path in utils_candidates:
    if utils_path.exists() and (utils_path / "partitioned_results_writer.py").exists():
        sys.path.append(str(utils_path))
        break
else:
    # Fallback: try to import directly (assume it's available in PYTHONPATH)
    pass

try:
    from partitioned_results_writer import (
        PartitionedResultsWriter,
        ConcurrentPartitionedWriter,
        create_simulation_partition_config
    )
except ImportError as e:
    # If import fails, check if utils files are in current directory
    if (current_dir / "partitioned_results_writer.py").exists():
        sys.path.append(str(current_dir))
        from partitioned_results_writer import (
            PartitionedResultsWriter,
            ConcurrentPartitionedWriter,
            create_simulation_partition_config
        )
    else:
        raise ImportError(f"Could not import partitioned_results_writer: {e}. "
                         f"Checked paths: {[str(p) for p in utils_candidates]}")

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# CAS9 Tier Configuration
class Cas9SimulationTier:
    def __init__(self, name: str, recording_sites: int, states_per_site: int, 
                 mutation_rate: float, description: str = ""):
        self.name = name
        self.recording_sites = recording_sites
        self.states_per_site = states_per_site
        self.mutation_rate = mutation_rate
        self.description = description

CAS9_TIERS = {
    1: Cas9SimulationTier("Ultra High Fidelity", 5000, 50, 0.5),
    2: Cas9SimulationTier("High Fidelity", 1000, 30, 0.4),
    3: Cas9SimulationTier("Medium Fidelity", 200, 20, 0.3),
    4: Cas9SimulationTier("Low Fidelity", 15, 10, 0.2)
}

# Solver definitions
# SOLVERS now imported from solver_config.py

def get_class_of_solver(solver_name):
    """Returns the solver class based on the solver name"""
    if solver_name == "nj":
        return cass.solver.NeighborJoiningSolver(add_root=True)
    elif solver_name == "maxcut":
        return cass.solver.MaxCutSolver()
    elif solver_name == "greedy":
        return cass.solver.VanillaGreedySolver()
    elif solver_name == "vanilla":
        return cass.solver.VanillaGreedySolver()
    elif solver_name == "spectral":
        return cass.solver.SpectralSolver()
    elif solver_name == "smj":
        return cass.solver.SharedMutationJoiningSolver()
    elif solver_name == "dmj":
        return cass.solver.DistanceSolver.DistanceSolver()
    elif solver_name == "ilp":
        return cass.solver.ILPSolver()
    else:
        raise ValueError(f"Unknown solver: {solver_name}")

def reconstruct_and_calculate_metrics(cas9_tree, solver_name: str, tier_num: int, 
                                     gt_instance_id: int = 0, cas9_simulation_id: int = 0, 
                                     reconstruction_id: int = 0):
    """Reconstruct tree with specified solver and calculate comprehensive metrics"""
    logger.info(f"Starting reconstruction with {solver_name} solver for tier {tier_num} "
               f"(GT:{gt_instance_id}, Sim:{cas9_simulation_id}, Recon:{reconstruction_id})")
    
    start_time = time.time()
    
    try:
        # Create a copy for reconstruction
        reconstructed_tree = deepcopy(cas9_tree)
        
        # Perform reconstruction
        solver = get_class_of_solver(solver_name)
        solver.solve(reconstructed_tree)
        
        # Post-process the reconstructed tree
        reconstructed_tree.reconstruct_ancestral_characters()
        
        # Use ConvexML for branch length estimation
        tree_newick = reconstructed_tree.get_newick(record_branch_lengths=True, record_node_names=True)
        leaf_sequences = {leaf: reconstructed_tree.get_character_states(leaf) for leaf in reconstructed_tree.leaves}
        
        # Estimate minimum branch length
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
        
        # Apply ConvexML branch length optimization
        try:
            convexml_result = convexml(
                tree_newick=tree_newick,
                leaf_sequences=leaf_sequences,
                minimum_branch_length=min_branch_len,
                resolve_multifurcations_before_branch_length_estimation=True
            )
            
            # Create optimized tree
            optimized_tree = cass.data.CassiopeiaTree(
                character_matrix=reconstructed_tree.character_matrix,
                tree=convexml_result["tree_newick"],
                missing_state_indicator=getattr(reconstructed_tree, 'missing_state_indicator', -1)
            )
        except (json.JSONDecodeError, ValueError, KeyError) as e:
            logger.warning(f"ConvexML branch length optimization failed: {e}. Using basic reconstruction.")
            # Fall back to basic reconstruction without branch length optimization
            optimized_tree = reconstructed_tree
            optimized_tree.scale_to_unit_length()
        
        # Copy priors and parameters
        if hasattr(reconstructed_tree, 'priors'):
            optimized_tree.priors = reconstructed_tree.priors
        
        # Reconstruct ancestral characters and scale
        optimized_tree.reconstruct_ancestral_characters()
        optimized_tree.scale_to_unit_length()
        
        # Calculate metrics
        computation_time = time.time() - start_time
        
        # Get parameters for PHS calculation - calculate BOTH methods for comparison

        # Method 1: From simulation (cas9_tree) - current implementation
        cm_sim = cas9_tree.character_matrix.to_numpy()
        proportion_mutated_sim = np.sum(cm_sim > 0) / np.sum(cm_sim >= 0)
        lam_from_simulation = -np.log(1.0 - proportion_mutated_sim)

        if hasattr(cas9_tree, 'priors') and cas9_tree.priors:
            priors_sim = np.array(list(cas9_tree.priors[0].values()))
            q_from_simulation = np.sum(priors_sim ** 2)
            q_simulation_source = "priors"
        else:
            q_from_simulation = 1.0 / cas9_tree.character_matrix.stack().nunique()
            q_simulation_source = "uniform_estimate"

        # Method 2: From ground truth (propagated through parameters)
        lam_from_gt = None
        q_from_gt = None
        proportion_mutated_gt = None

        if hasattr(cas9_tree, 'parameters'):
            lam_from_gt = cas9_tree.parameters.get("lam_gt")
            q_from_gt = cas9_tree.parameters.get("q_gt")
            proportion_mutated_gt = cas9_tree.parameters.get("proportion_mutated_gt")

        # Use simulation values for actual computation (maintaining current behavior)
        # but track both for comparison
        lam = lam_from_simulation
        q = q_from_simulation
        
        metrics = {
            'reconstruction_id': f"instance{gt_instance_id}_sim{cas9_simulation_id}_recon{reconstruction_id}_tier{tier_num}_{solver_name}",
            'gt_instance_id': gt_instance_id,
            'cas9_simulation_id': cas9_simulation_id,
            'reconstruction_num': reconstruction_id,
            'cas9_tier': tier_num,
            'cas9_tier_name': f'Tier {tier_num}',  # Simplified for now
            'recording_sites': 100,  # Default fallback
            'states_per_site': 10,   # Default fallback
            'solver': solver_name,
            'computation_time_seconds': computation_time,

            # PHS parameters from simulation (used for computation)
            'lam_simulation': lam_from_simulation,
            'q_simulation': q_from_simulation,
            'proportion_mutated_simulation': proportion_mutated_sim,
            'q_simulation_source': q_simulation_source,

            # PHS parameters from ground truth (for comparison)
            'lam_gt': lam_from_gt,
            'q_gt': q_from_gt,
            'proportion_mutated_gt': proportion_mutated_gt,

            # Track which parameters were actually used for PHS calculation
            'phs_lam_source': 'simulation',  # Could be changed to 'gt' in future
            'phs_q_source': 'simulation',    # Could be changed to 'gt' in future
        }
        
        try:
            # Triplets correctness - use all available cores
            max_threads = os.environ.get('LSB_MAX_NUM_PROCESSORS', '20')
            triplets_result = stellars.triplets_correct_parallel(
                cas9_tree.get_newick(record_branch_lengths=True, record_node_names=True),
                optimized_tree.get_newick(record_branch_lengths=True, record_node_names=True),
                number_of_trials=1000,
                min_triplets_at_depth=1,
                seed=42,
                max_threads=int(max_threads)
            )
            metrics['triplets_distance'] = 1.0 - triplets_result['all_triplets_correct'][0]
        except Exception as e:
            logger.warning(f"Triplets calculation failed: {e}")
            metrics['triplets_distance'] = np.nan
        
        try:
            # Robinson-Foulds distance
            rf_result = stellars.rf_distance(
                cas9_tree.get_newick(record_branch_lengths=True, record_node_names=True),
                optimized_tree.get_newick(record_branch_lengths=True, record_node_names=True)
            )
            metrics['RF_distance'] = rf_result.normalized_distance
        except Exception as e:
            logger.warning(f"RF distance calculation failed: {e}")
            metrics['RF_distance'] = np.nan
        
        try:
            # Parsimony score using stellars (Rust-based, faster)
            metrics['parsimony_score'] = stellars.parsimony_score(
                tree_newick=optimized_tree.get_newick(record_branch_lengths=True, record_node_names=True),
                character_matrix=optimized_tree.character_matrix.to_numpy(),
                missing_state=-1,
                unedited_state=0
            )
        except Exception as e:
            logger.warning(f"Parsimony calculation failed: {e}")
            metrics['parsimony_score'] = np.nan
        
        try:
            # cPHS calculation using CORRECTED approach from simulation_phs.py (lines 336-381)
            max_threads = os.environ.get('LSB_MAX_NUM_PROCESSORS', '20')
            tree_newick = optimized_tree.get_newick(record_branch_lengths=True, record_node_names=True)

            # CRITICAL FIX: Extract leaf names for proper character matrix mapping (simulation_phs.py line 340)
            # This fixes the character state mapping issue that was causing ultra-low p-values
            leaf_names = list(optimized_tree.character_matrix.index)
            logger.info(f"Using {len(leaf_names)} leaf names for proper character matrix mapping")

            # Correct character matrix format (simulation_phs.py line 336)
            character_matrix = optimized_tree.character_matrix.values.astype(int).tolist()

            # Extract internal character states from Cassiopeia (simulation_phs.py lines 343-361)
            internal_character_states = {}
            for internal_node in optimized_tree.internal_nodes:
                node_name = str(internal_node)
                try:
                    internal_states = optimized_tree.get_character_states(internal_node)
                    if hasattr(internal_states, 'tolist'):
                        internal_character_states[node_name] = internal_states.tolist()
                    else:
                        internal_character_states[node_name] = list(internal_states)
                except Exception:
                    # If extraction fails for any node, clear all and let stellars infer
                    internal_character_states = {}
                    break

            logger.info(f"Extracted {len(internal_character_states)} internal node character states")

            # Calculate cPHS with simulation parameters using CORRECTED method (simulation_phs.py line 371)
            phs_result_sim = stellars.phs_optimized(
                tree_newick=tree_newick,
                character_matrix=character_matrix,
                internal_character_states=internal_character_states,
                mutation_rate=lam_from_simulation,
                collision_probability=q_from_simulation,
                missing_state=-1,
                unedited_state=0,
                use_provided_internal_states=True,  # Use Cassiopeia's exact internal states
                leaf_names=leaf_names  # CRITICAL FIX: Proper character matrix mapping
            )
            metrics['cPHS_simulation'] = phs_result_sim['phs_score']

            # Calculate cPHS with ground truth parameters using CORRECTED method if available
            if lam_from_gt is not None and q_from_gt is not None:
                phs_result_gt = stellars.phs_optimized(
                    tree_newick=tree_newick,
                    character_matrix=character_matrix,
                    internal_character_states=internal_character_states,
                    mutation_rate=lam_from_gt,
                    collision_probability=q_from_gt,
                    missing_state=-1,
                    unedited_state=0,
                    use_provided_internal_states=True,
                    leaf_names=leaf_names
                )
                metrics['cPHS_gt'] = phs_result_gt['phs_score']
            else:
                metrics['cPHS_gt'] = np.nan

            # Debug logging
            logger.info(f"CORRECTED cPHS - Simulation: {metrics['cPHS_simulation']:.2e}, "
                       f"Ground Truth: {metrics.get('cPHS_gt', 'N/A')}")

            # Legacy compatibility field
            metrics['cPHS'] = metrics['cPHS_simulation']

        except Exception as e:
            logger.warning(f"CORRECTED cPHS calculation failed: {e}")
            metrics['cPHS'] = np.nan
            metrics['cPHS_simulation'] = np.nan
            metrics['cPHS_gt'] = np.nan
        
        try:
            # Likelihood score using stellars - calculate for both parameter sets
            tree_newick = optimized_tree.get_newick(record_branch_lengths=True, record_node_names=True)
            character_matrix = optimized_tree.character_matrix.to_numpy()

            # Calculate likelihood with simulation parameters
            likelihood_result_sim = stellars.likelihood_score(
                tree_newick=tree_newick,
                character_matrix=character_matrix,
                mutation_rate=lam_from_simulation,
                collision_probability=q_from_simulation,
                missing_state=-1,
                unedited_state=0
            )
            metrics['likelihood_score_simulation'] = likelihood_result_sim

            # Calculate likelihood with ground truth parameters if available
            if lam_from_gt is not None and q_from_gt is not None:
                likelihood_result_gt = stellars.likelihood_score(
                    tree_newick=tree_newick,
                    character_matrix=character_matrix,
                    mutation_rate=lam_from_gt,
                    collision_probability=q_from_gt,
                    missing_state=-1,
                    unedited_state=0
                )
                metrics['likelihood_score_gt'] = likelihood_result_gt
            else:
                metrics['likelihood_score_gt'] = np.nan

            # Keep legacy likelihood_score field for backward compatibility (using simulation parameters)
            metrics['likelihood_score'] = likelihood_result_sim

        except Exception as e:
            logger.warning(f"Likelihood calculation failed: {e}")
            metrics['likelihood_score'] = np.nan
            metrics['likelihood_score_simulation'] = np.nan
            metrics['likelihood_score_gt'] = np.nan
        
        logger.info(f"Reconstruction completed successfully: {solver_name} on tier {tier_num}")
        triplets_val = metrics.get('triplets_distance', 'N/A')
        triplets_str = f"{triplets_val:.3f}" if triplets_val != 'N/A' else 'N/A'
        logger.info(f"Key metrics - RF: {metrics.get('RF_distance', 'N/A')}, "
                   f"Triplets: {triplets_str}, "
                   f"cPHS: {metrics.get('cPHS', 'N/A')}")
        
        return metrics
        
    except Exception as e:
        logger.error(f"Reconstruction failed: {e}")
        computation_time = time.time() - start_time

        # Return result with same schema as successful reconstruction to avoid DataFrame mismatch
        # Use fallback values for tier config
        return {
            'reconstruction_id': f"instance{gt_instance_id}_sim{cas9_simulation_id}_recon{reconstruction_id}_tier{tier_num}_{solver_name}",
            'gt_instance_id': gt_instance_id,
            'cas9_simulation_id': cas9_simulation_id,
            'reconstruction_num': reconstruction_id,
            'cas9_tier': tier_num,
            'cas9_tier_name': tier_config.name,
            'recording_sites': tier_config.recording_sites,
            'states_per_site': tier_config.states_per_site,
            'solver': solver_name,
            'computation_time_seconds': computation_time,
            'lam_simulation': None,
            'q_simulation': None,
            'proportion_mutated_simulation': None,
            'q_simulation_source': None,
            'lam_gt': None,
            'q_gt': None,
            'proportion_mutated_gt': None,
            'phs_lam_source': None,
            'phs_q_source': None,
            'triplets_distance': None,
            'RF_distance': None,
            'cPHS': None,
            'cPHS_simulation': None,
            'cPHS_gt': None,
            'parsimony_score': None,
            'total_mutations': None,
            'parsimony_computation_time_ms': None,
            'parsimony_method': None,
            'log_likelihood': None,
            'likelihood': None,
            'likelihood_computation_time_ms': None,
            'likelihood_method': None,
            'likelihood_score': None,
            'likelihood_score_simulation': None,
            'likelihood_score_gt': None,
            'n_characters': None,
            'n_leaves': None,
            'error': str(e),
            'status': 'failed'
        }


class ReconstructionWorker:
    """Worker that performs reconstruction and calculates metrics."""
    
    def __init__(self, cas9_instance_path: str, solver: str, tier: int, 
                 output_dir: str, shared_dir: str, gt_instance_id: int = 0, 
                 cas9_simulation_id: int = 0, reconstruction_id: int = 0):
        self.cas9_instance_path = Path(cas9_instance_path)
        self.solver = solver
        self.tier = tier
        self.output_dir = Path(output_dir)
        self.shared_dir = Path(shared_dir)
        self.gt_instance_id = gt_instance_id
        self.cas9_simulation_id = cas9_simulation_id
        self.reconstruction_id = reconstruction_id
        
        # Ensure directories exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate inputs
        if solver not in SOLVERS:
            raise ValueError(f"Invalid solver: {solver}. Must be one of {SOLVERS}")
        # Load config to validate tier against available tiers
        config_path = Path(shared_dir) / "cascade_config.yaml"
        config = {}
        if config_path.exists():
            try:
                import yaml
                with open(config_path, 'r') as f:
                    config = yaml.safe_load(f)
            except Exception as e:
                logger.warning(f"Could not load config {config_path}: {e}")

        config_tiers = list(config.get('cas9_tiers', {}).keys())
        if config_tiers and tier not in config_tiers:
            raise ValueError(f"Invalid tier: {tier}. Must be one of {config_tiers}")
            
        # Initialize partitioned results writer
        partition_config = create_simulation_partition_config(
            output_dir=self.shared_dir / "partitioned_results",
            buffer_size=100  # Smaller buffer for more frequent flushes in concurrent env
        )
        self.partitioned_writer = ConcurrentPartitionedWriter(partition_config)
    
    def load_cas9_instance(self):
        """Load the Cas9 instance tree."""
        logger.info(f"Loading Cas9 instance from: {self.cas9_instance_path}")
        
        if not self.cas9_instance_path.exists():
            raise FileNotFoundError(f"Cas9 instance file not found: {self.cas9_instance_path}")
        
        with open(self.cas9_instance_path, 'rb') as f:
            cas9_tree = pickle.load(f)
        
        return cas9_tree
    
    def perform_reconstruction_and_analysis(self):
        """Perform reconstruction and calculate all metrics."""
        logger.info(f"Performing reconstruction with {self.solver} for Tier {self.tier}...")
        
        start_time = time.time()
        
        try:
            # Load Cas9 instance
            cas9_tree = self.load_cas9_instance()
            
            # Get tier configuration for metadata from config
            config_path = self.shared_dir / "cascade_config.yaml"
            config = {}
            if config_path.exists():
                try:
                    import yaml
                    with open(config_path, 'r') as f:
                        config = yaml.safe_load(f)
                except Exception as e:
                    logger.warning(f"Could not load config {config_path}: {e}")

            tier_config_dict = config.get('cas9_tiers', {}).get(self.tier, {})
            if not tier_config_dict:
                raise ValueError(f"Tier {self.tier} configuration not found in config")

            # Create a simple tier config object for metadata
            class TierConfig:
                def __init__(self, tier_num, **kwargs):
                    self.tier_num = tier_num
                    self.name = kwargs.get('name', f'Tier {tier_num}')
                    self.k = kwargs.get('k', 10)
                    self.cassette_size = kwargs.get('cassette_size', 3)
                    self.recording_sites = self.k * self.cassette_size
                    self.states_per_site = kwargs.get('m', 10)

            tier_config = TierConfig(self.tier, **tier_config_dict)
            
            # Perform reconstruction and calculate metrics
            result = reconstruct_and_calculate_metrics(
                cas9_tree=cas9_tree,
                solver_name=self.solver,
                tier_num=self.tier,
                gt_instance_id=self.gt_instance_id,
                cas9_simulation_id=self.cas9_simulation_id,
                reconstruction_id=self.reconstruction_id
            )
            
            # Add timing information
            result['computation_time_seconds'] = time.time() - start_time
            result['worker_info'] = {
                'cas9_instance_path': str(self.cas9_instance_path),
                'solver': self.solver,
                'tier': self.tier,
                'tier_name': tier_config.name
            }
            
            # Save to partitioned parquet structure - create separate rows for each parameter source
            try:
                # Create two separate result rows - one for simulation, one for ground truth

                # 1. Simulation parameters result
                result_simulation = result.copy()
                result_simulation['phs_lam_source'] = 'simulation'
                result_simulation['phs_q_source'] = 'simulation'
                result_simulation['cPHS'] = result.get('cPHS_simulation', result.get('cPHS'))  # Use simulation PHS

                flattened_result_sim = self.flatten_result_for_parquet(result_simulation)
                self.partitioned_writer.add_result(flattened_result_sim)

                # 2. Ground truth parameters result (only if GT parameters are available)
                if result.get('cPHS_gt') is not None and not np.isnan(result.get('cPHS_gt', np.nan)):
                    result_gt = result.copy()
                    result_gt['phs_lam_source'] = 'ground_truth'
                    result_gt['phs_q_source'] = 'ground_truth'
                    result_gt['cPHS'] = result.get('cPHS_gt')  # Use ground truth PHS

                    flattened_result_gt = self.flatten_result_for_parquet(result_gt)
                    self.partitioned_writer.add_result(flattened_result_gt)

                    logger.info(f"Results added to partitioned storage (tier={self.tier}, solver={self.solver}) - both simulation and ground truth")
                else:
                    logger.info(f"Result added to partitioned storage (tier={self.tier}, solver={self.solver}) - simulation only (GT parameters unavailable)")

            except Exception as e:
                logger.warning(f"Failed to save to partitioned storage: {e}")
            
            logger.info(f"Results saved to partitioned storage")
            triplets_val = result.get('triplets_distance', 'N/A')
            cphs_val = result.get('cPHS', 'N/A')
            triplets_str = f"{triplets_val:.4f}" if triplets_val != 'N/A' else 'N/A'
            cphs_str = f"{cphs_val:.6f}" if cphs_val != 'N/A' else 'N/A'
            logger.info(f"Key metrics - RF: {result.get('RF_distance', 'N/A')}, "
                       f"Triplets: {triplets_str}, "
                       f"cPHS: {cphs_str}")
            
            # Update status
            # self.update_status("level3", f"tier{self.tier}_{self.solver}", "completed", {
            #     'computation_time': result['computation_time_seconds'],
            #     'rf_distance': result.get('RF_distance'),
            #     'triplets_distance': result.get('triplets_distance'),
            #     'cphs_pvalue': result.get('cPHS')
            # })
            
            return result
            
        except Exception as e:
            logger.error(f"Failed reconstruction with {self.solver} for Tier {self.tier}: {e}")
            # self.update_status("level3", f"tier{self.tier}_{self.solver}", "failed", {
            #     'error': str(e),
            #     'computation_time': time.time() - start_time
            # })
            raise
    
    def flatten_result_for_parquet(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Flatten nested result structure for efficient parquet storage."""
        flattened = {}
        
        # Copy basic fields
        for key, value in result.items():
            if key in ['worker_info', 'parsimony_score', 'likelihood_score', 'likelihood_score_simulation', 'likelihood_score_gt']:
                continue  # Handle these separately
            flattened[key] = value
        
        # Flatten worker_info
        if 'worker_info' in result:
            for key, value in result['worker_info'].items():
                flattened[f'worker_{key}'] = value
        
        # Flatten parsimony_score if it's a dict (stellars format)
        if 'parsimony_score' in result:
            pars = result['parsimony_score']
            if isinstance(pars, dict):
                flattened['parsimony_score'] = pars.get('parsimony_score', pars)
                flattened['total_mutations'] = pars.get('total_mutations', None)
                flattened['parsimony_computation_time_ms'] = pars.get('computation_time_ms', None)
                flattened['parsimony_method'] = pars.get('method_used', None)
            else:
                flattened['parsimony_score'] = pars
        
        # Flatten likelihood_score if it's a dict (stellars format)
        if 'likelihood_score' in result:
            like = result['likelihood_score']
            if isinstance(like, dict):
                flattened['log_likelihood'] = like.get('log_likelihood', like)
                flattened['likelihood'] = like.get('likelihood', None)
                flattened['likelihood_computation_time_ms'] = like.get('computation_time_ms', None)
                flattened['likelihood_method'] = like.get('method_used', None)
                flattened['n_characters'] = like.get('n_characters', None)
                flattened['n_leaves'] = like.get('n_leaves', None)
            else:
                flattened['log_likelihood'] = like

        # Flatten likelihood_score_simulation if it's a dict (stellars format)
        if 'likelihood_score_simulation' in result:
            like_sim = result['likelihood_score_simulation']
            if isinstance(like_sim, dict):
                flattened['log_likelihood_simulation'] = like_sim.get('log_likelihood', like_sim)
                flattened['likelihood_simulation'] = like_sim.get('likelihood', None)
                flattened['likelihood_computation_time_ms_simulation'] = like_sim.get('computation_time_ms', None)
                flattened['likelihood_method_simulation'] = like_sim.get('method_used', None)
                flattened['n_characters_simulation'] = like_sim.get('n_characters', None)
                flattened['n_leaves_simulation'] = like_sim.get('n_leaves', None)
            else:
                flattened['log_likelihood_simulation'] = like_sim

        # Flatten likelihood_score_gt if it's a dict (stellars format)
        if 'likelihood_score_gt' in result:
            like_gt = result['likelihood_score_gt']
            if isinstance(like_gt, dict):
                flattened['log_likelihood_gt'] = like_gt.get('log_likelihood', like_gt)
                flattened['likelihood_gt'] = like_gt.get('likelihood', None)
                flattened['likelihood_computation_time_ms_gt'] = like_gt.get('computation_time_ms', None)
                flattened['likelihood_method_gt'] = like_gt.get('method_used', None)
                flattened['n_characters_gt'] = like_gt.get('n_characters', None)
                flattened['n_leaves_gt'] = like_gt.get('n_leaves', None)
            else:
                flattened['log_likelihood_gt'] = like_gt

        return flattened
    
    def update_status(self, level: str, component: str, status: str, details: Any = None) -> None:
        """Update job status in shared status file."""
        status_file = self.shared_dir / "status" / f"{level}_status.json"
        temp_file = status_file.with_suffix('.json.tmp')
        
        max_retries = 3
        for attempt in range(max_retries):
            try:
                # Load existing status
                if status_file.exists():
                    with open(status_file, 'r') as f:
                        status_data = json.load(f)
                else:
                    status_data = {}
                
                # Update status
                status_data[component] = {
                    'status': status,
                    'timestamp': str(time.time()),
                    'details': details
                }
                
                # Atomic write using temporary file
                with open(temp_file, 'w') as f:
                    json.dump(status_data, f, indent=2)
                
                # Atomic move
                temp_file.rename(status_file)
                return
                
            except (json.JSONDecodeError, FileNotFoundError, PermissionError) as e:
                logger.warning(f"Status update attempt {attempt + 1} failed: {e}")
                if attempt < max_retries - 1:
                    time.sleep(0.1 * (attempt + 1))  # Exponential backoff
                else:
                    logger.error(f"Failed to update status after {max_retries} attempts: {e}")
            finally:
                # Clean up temp file if it exists
                if temp_file.exists():
                    try:
                        temp_file.unlink()
                    except:
                        pass
    
    def aggregate_to_consolidated(self, flattened_result: Dict[str, Any]):
        """Legacy method - now handled by partitioned writer"""
        # This method is kept for compatibility but partitioned writer handles aggregation
        pass
    
    def run(self) -> None:
        """Execute the reconstruction workflow."""
        logger.info(f"=== Reconstruction Worker Starting - Tier {self.tier}, Solver {self.solver} ===")
        
        try:
            # Perform reconstruction and analysis
            result = self.perform_reconstruction_and_analysis()
            
            logger.info(f"=== Reconstruction Worker Completed Successfully - "
                       f"Tier {self.tier}, Solver {self.solver} ===")
            
            return result
            
        except Exception as e:
            logger.error(f"Reconstruction Worker failed for Tier {self.tier}, Solver {self.solver}: {e}")
            sys.exit(1)
        finally:
            # Ensure partitioned writer buffer is flushed
            if hasattr(self, 'partitioned_writer'):
                self.partitioned_writer.flush()


def main():
    parser = argparse.ArgumentParser(description="Reconstruction Worker - Perform reconstruction and calculate metrics")
    parser.add_argument('--cas9_instance_path', required=True, help='Path to Cas9 instance tree')
    parser.add_argument('--solver', required=True, choices=SOLVERS, help='Reconstruction solver to use')
    parser.add_argument('--tier', type=int, required=True, help='Cas9 recording tier (1-4)')
    parser.add_argument('--output_dir', required=True, help='Directory to save results')
    parser.add_argument('--shared_dir', required=True, 
                       help='Shared directory for job coordination')
    parser.add_argument('--gt_instance_id', type=int, default=0, 
                       help='Ground truth instance ID')
    parser.add_argument('--cas9_simulation_id', type=int, default=0,
                       help='Cas9 simulation ID for this GT instance')
    parser.add_argument('--reconstruction_id', type=int, default=0,
                       help='Reconstruction ID for this Cas9 simulation')
    
    args = parser.parse_args()
    
    worker = ReconstructionWorker(
        args.cas9_instance_path,
        args.solver,
        args.tier,
        args.output_dir,
        args.shared_dir,
        args.gt_instance_id,
        args.cas9_simulation_id,
        args.reconstruction_id
    )
    worker.run()


if __name__ == "__main__":
    main()
