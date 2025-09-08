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
import cassiopeia as cass
from copy import deepcopy
import stellars
from convexml import convexml
from solver_config import get_solver_class, SOLVERS
import yaml

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
        
        # Get parameters for PHS calculation
        cm = cas9_tree.character_matrix.to_numpy()
        proportion_mutated = np.sum(cm > 0) / np.sum(cm >= 0)
        lam = -np.log(1.0 - proportion_mutated)
        
        # Get collision probability
        if hasattr(cas9_tree, 'priors') and cas9_tree.priors:
            priors = np.array(list(cas9_tree.priors[0].values()))
            q = np.sum(priors ** 2)
        else:
            q = 1.0 / cas9_tree.character_matrix.stack().nunique()
        
        metrics = {
            'reconstruction_id': f"instance{gt_instance_id}_sim{cas9_simulation_id}_recon{reconstruction_id}_tier{tier_num}_{solver_name}",
            'gt_instance_id': gt_instance_id,
            'cas9_simulation_id': cas9_simulation_id,
            'reconstruction_num': reconstruction_id,
            'cas9_tier': tier_num,
            'cas9_tier_name': CAS9_TIERS[tier_num].name,
            'recording_sites': CAS9_TIERS[tier_num].recording_sites,
            'solver': solver_name,
            'computation_time_seconds': computation_time,
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
            # cPHS calculation using stellars - use all available cores
            max_threads = os.environ.get('LSB_MAX_NUM_PROCESSORS', '20')
            phs_result = stellars.phs_from_cassiopeia(
                optimized_tree,
                mutation_rate=lam,
                collision_probability=q,
                reconstruct_ancestral=False,  # Already done
                max_threads=int(max_threads)
            )
            metrics['cPHS'] = phs_result['phs_score']
        except Exception as e:
            logger.warning(f"cPHS calculation failed: {e}")
            metrics['cPHS'] = np.nan
        
        try:
            # Likelihood score using stellars (more reliable than Cassiopeia's)
            metrics['likelihood_score'] = stellars.likelihood_score(
                tree_newick=optimized_tree.get_newick(record_branch_lengths=True, record_node_names=True),
                character_matrix=optimized_tree.character_matrix.to_numpy(),
                mutation_rate=lam,
                collision_probability=q,
                missing_state=-1,  # Standard missing state indicator
                unedited_state=0   # Standard unedited state
            )
        except Exception as e:
            logger.warning(f"Likelihood calculation failed: {e}")
            metrics['likelihood_score'] = np.nan
        
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
        return {
            'reconstruction_id': f"instance{gt_instance_id}_sim{cas9_simulation_id}_recon{reconstruction_id}_tier{tier_num}_{solver_name}",
            'gt_instance_id': gt_instance_id,
            'cas9_simulation_id': cas9_simulation_id,
            'reconstruction_num': reconstruction_id,
            'cas9_tier': tier_num,
            'solver': solver_name,
            'error': str(e),
            'computation_time_seconds': computation_time,
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
        if tier not in CAS9_TIERS:
            raise ValueError(f"Invalid tier: {tier}. Must be one of {list(CAS9_TIERS.keys())}")
    
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
            
            # Get tier configuration for metadata
            tier_config = CAS9_TIERS[self.tier]
            
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
            
            # Save result as JSON with proper indexing
            result_filename = f"instance{self.gt_instance_id}_sim{self.cas9_simulation_id}_recon{self.reconstruction_id}_tier{self.tier}_{self.solver}_metrics.json"
            result_path = self.output_dir / result_filename
            
            with open(result_path, 'w') as f:
                json.dump(result, f, indent=2, default=str)
                
            # Also save as Parquet for efficient analysis
            try:
                # Flatten the result for tabular format
                flattened_result = self.flatten_result_for_parquet(result)
                df = pl.DataFrame([flattened_result])
                
                parquet_filename = f"instance{self.gt_instance_id}_sim{self.cas9_simulation_id}_recon{self.reconstruction_id}_tier{self.tier}_{self.solver}_metrics.parquet"
                parquet_path = self.output_dir / parquet_filename
                df.write_parquet(parquet_path)
                
                logger.info(f"Parquet saved to: {parquet_path}")
            except Exception as e:
                logger.warning(f"Failed to save parquet: {e}")
            
            logger.info(f"Results saved to: {result_path}")
            triplets_val = result.get('triplets_distance', 'N/A')
            cphs_val = result.get('cPHS', 'N/A')
            triplets_str = f"{triplets_val:.4f}" if triplets_val != 'N/A' else 'N/A'
            cphs_str = f"{cphs_val:.6f}" if cphs_val != 'N/A' else 'N/A'
            logger.info(f"Key metrics - RF: {result.get('RF_distance', 'N/A')}, "
                       f"Triplets: {triplets_str}, "
                       f"cPHS: {cphs_str}")
            
            # Update status
            # self.update_status("level3", f"tier{self.tier}_{self.solver}", "completed", {
            #     'result_path': str(result_path),
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
            if key in ['worker_info', 'parsimony_score', 'likelihood_score']:
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
