use phylo::prelude::*;
use phylo::node::NodeID;
// use rayon::prelude::*;  // Removed: not currently used
use fxhash::FxHashMap;
use std::collections::HashMap;
use statrs::distribution::{Binomial, DiscreteCDF};
use crate::parsimony::fitch_parsimony;
use log::{debug, info, warn};

/// PHS calculation result using optimized techniques
#[derive(Debug, Clone)]
pub struct PHSResult {
    pub phs_score: f64,
    pub total_pairs: usize,
    pub computation_time_ms: f64,
    pub method_used: String,
    pub parallel_chunks_used: usize,
}

/// Simple PHS calculation following Python implementation exactly
struct PHSTreeData {
    // Character states for leaves (leaf_name -> character_states)
    leaf_character_states: FxHashMap<String, Vec<i32>>,
    
    // Inferred character states for internal nodes (node_id -> character_states)
    internal_character_states: FxHashMap<NodeID, Vec<i32>>,
    
    // Parameters
    missing_state: i32,
    unedited_state: i32,
    n_characters: usize,
}

impl PHSTreeData {
    fn new(
        tree: &mut PhyloTree,
        character_matrix: Vec<Vec<i32>>,
        internal_character_states: HashMap<String, Vec<i32>>,
        missing_state: i32,
        unedited_state: i32,
        use_provided_internal_states: bool,
        leaf_names: Option<Vec<String>>,
    ) -> Self {
        tree.precompute_constant_time_lca();
        
        // Build leaf character states map
        // PROPER FIX: Map character states by leaf name instead of assuming index correspondence
        let mut leaf_character_states = FxHashMap::default();
        
        if let Some(names) = &leaf_names {
            // Use provided leaf names to map character_matrix rows to leaf names
            info!("Using provided leaf names for character matrix mapping");
            if names.len() != character_matrix.len() {
                warn!("Leaf names length ({}) does not match character matrix length ({})", 
                    names.len(), character_matrix.len());
            }
            
            for (idx, leaf_name) in names.iter().enumerate() {
                if idx < character_matrix.len() {
                    leaf_character_states.insert(leaf_name.clone(), character_matrix[idx].clone());
                    debug!("CORRECT MAPPING: Leaf '{}' gets character_matrix[{}] = {:?}", 
                        leaf_name, idx, &character_matrix[idx][..5.min(character_matrix[idx].len())]);
                }
            }
            
            info!("Successfully mapped {} leaf character states by name", leaf_character_states.len());
        } else {
            // Fallback to old behavior with warning
            warn!("CRITICAL BUG: No leaf names provided! Using tree traversal order (likely incorrect)");
            warn!("This may cause homoplasy overcounting due to character state mismatches");
            
            for (idx, leaf) in tree.get_leaves().enumerate() {
                if let Some(taxa) = leaf.get_taxa() {
                    if idx < character_matrix.len() {
                        leaf_character_states.insert(taxa.clone(), character_matrix[idx].clone());
                        debug!("FALLBACK MAPPING: Leaf '{}' at tree position {} gets character_matrix[{}]", 
                            taxa, idx, idx);
                    }
                }
            }
        }
        
        let n_characters = if character_matrix.is_empty() { 
            0 
        } else { 
            character_matrix[0].len() 
        };
        
        let internal_char_states = if use_provided_internal_states && !internal_character_states.is_empty() {
            info!("Using provided internal character states from Cassiopeia");
            map_cassiopeia_to_phylo_states(tree, internal_character_states, n_characters, unedited_state)
        } else {
            info!("Reconstructing ancestral states using Fitch parsimony algorithm");
            if !internal_character_states.is_empty() {
                debug!("Python provided {} internal states, but using Fitch reconstruction because use_provided_internal_states=false", 
                    internal_character_states.len());
            }
            fitch_parsimony(tree, &character_matrix, -1, unedited_state)
        };
        
        debug!("Using {} internal nodes", internal_char_states.len());
        for (node_id, states) in &internal_char_states {
            debug!("rigi Node {}: {:?}", node_id, states);
        }
        
        PHSTreeData {
            leaf_character_states,
            internal_character_states: internal_char_states,
            missing_state,
            unedited_state,
            n_characters,
        }
    }
}

/// Map Cassiopeia internal node states to phylo-rs NodeIDs using node names from Newick
/// This uses the actual node names that were preserved during Newick parsing
fn map_cassiopeia_to_phylo_states(
    tree: &PhyloTree,
    cassiopeia_states: HashMap<String, Vec<i32>>,
    n_characters: usize,
    unedited_state: i32,
) -> FxHashMap<NodeID, Vec<i32>> {
    let mut phylo_states = FxHashMap::default();
    
    debug!("Mapping {} Cassiopeia internal states to phylo-rs NodeIDs using node names", cassiopeia_states.len());
    
    // Create a mapping from node names (taxa) to NodeIDs in phylo-rs tree
    let mut name_to_node_id = HashMap::new();
    
    // Get all nodes and create name -> NodeID mapping
    for node in tree.get_nodes() {
        if !node.is_leaf() {
            if let Some(taxa) = node.get_taxa() {
                name_to_node_id.insert(taxa.clone(), node.get_id());
                debug!("Found internal node with name '{}' mapped to NodeID {}", taxa, node.get_id());
            } else {
                debug!("Found internal node with NodeID {} but no name", node.get_id());
            }
            
            // Check if branch lengths are available
            if let Some(parent_id) = node.get_parent() {
                if let Some(branch_length) = tree.get_edge_weight(parent_id, node.get_id()) {
                    debug!("  Branch length from parent {} to node {}: {}", parent_id, node.get_id(), branch_length);
                } else {
                    debug!("  No branch length available from parent {} to node {}", parent_id, node.get_id());
                }
            } else {
                debug!("  Node {} is root (no parent)", node.get_id());
            }
        }
    }
    
    debug!("Created name mapping for {} internal nodes in phylo-rs tree", name_to_node_id.len());
    
    // Map Cassiopeia states to phylo-rs NodeIDs using the name mapping
    let mut successful_mappings = 0;
    let total_cassiopeia_states = cassiopeia_states.len();
    
    for (cassiopeia_name, states) in cassiopeia_states {
        if let Some(&node_id) = name_to_node_id.get(&cassiopeia_name) {
            // Ensure the states vector has the right length
            let mut padded_states = states;
            padded_states.resize(n_characters, unedited_state);
            
            phylo_states.insert(node_id, padded_states);
            successful_mappings += 1;
            
            debug!("Successfully mapped Cassiopeia node '{}' to phylo-rs NodeID {} with {} characters", 
                cassiopeia_name, node_id, n_characters);
        } else {
            warn!("Could not find phylo-rs node with name '{}' (available names: {:?})", 
                cassiopeia_name, name_to_node_id.keys().take(5).collect::<Vec<_>>());
        }
    }
    
    info!("Successfully mapped {}/{} Cassiopeia internal states to phylo-rs using node names", 
        successful_mappings, total_cassiopeia_states);
    
    phylo_states
}

/// Calculate the distance from root matching Python/Cassiopeia tree.get_time()
/// Python uses raw distance from root on trees that are already unit-length scaled
fn calculate_distance_from_root(tree: &PhyloTree, node_id: NodeID) -> f64 {
    // Calculate distance from root to this node (matching Python's tree.get_time())
    let mut current_node_id = node_id;
    let mut total_distance = 0.0;
    let mut path_details = Vec::new();
    
    // Traverse from current node to root, summing branch lengths
    while let Some(parent_id) = tree.get_node_parent_id(current_node_id) {
        if let Some(branch_length) = tree.get_edge_weight(parent_id, current_node_id) {
            total_distance += branch_length as f64;
            path_details.push((current_node_id, parent_id, branch_length));
        } else {
            path_details.push((current_node_id, parent_id, 0.0));
        }
        current_node_id = parent_id;
    }
    
    // Only log for debugging specific cases - can be disabled in production
    if log::log_enabled!(log::Level::Debug) && !path_details.is_empty() {
        debug!("Node {} distance calculation:", node_id);
        for (child, parent, length) in &path_details {
            debug!("  {} -> {} : length = {:.6}", child, parent, length);
        }
        debug!("  Total distance from root: {:.6}", total_distance);
    }
    
    total_distance
}

/// Calculate the maximum distance from root to any leaf in the tree
fn calculate_max_root_to_leaf_distance(tree: &PhyloTree) -> f64 {
    let mut max_distance: f64 = 0.0;
    
    for leaf in tree.get_leaves() {
        let mut current_node_id = leaf.get_id();
        let mut total_distance = 0.0;
        
        // Traverse from leaf to root, summing branch lengths
        while let Some(parent_id) = tree.get_node_parent_id(current_node_id) {
            if let Some(branch_length) = tree.get_edge_weight(parent_id, current_node_id) {
                total_distance += branch_length as f64;
            }
            current_node_id = parent_id;
        }
        
        max_distance = max_distance.max(total_distance);
    }
    
    max_distance
}

/// Calculate homoplasies following the Python implementation exactly
/// Python counts each leaf pair exactly once with their LCA
fn calculate_phs_and_lca_heights_direct(
    tree: &PhyloTree,
    tree_data: &PHSTreeData,
) -> (Vec<usize>, Vec<f64>) {
    let mut homoplasy_counts = Vec::new();
    let mut lca_heights = Vec::new();
    
    // Get all leaves once
    let all_leaves: Vec<_> = tree.get_leaves().collect();
    
    info!("=== STELLARS DETAILED HOMOPLASY CALCULATION ===");
    info!("Tree has {} leaves, processing {} pairs", all_leaves.len(), all_leaves.len() * (all_leaves.len() - 1) / 2);
    info!("Characters per leaf: {}", tree_data.n_characters);
    
    let mut total_homoplasies = 0;
    let mut pair_count = 0;
    
    // Python approach: iterate through all leaf pairs exactly once
    for i in 0..all_leaves.len() {
        for j in i+1..all_leaves.len() {
            let leaf1 = &all_leaves[i];
            let leaf2 = &all_leaves[j];
            
            // Find the LCA of this pair
            let lca_id = tree.get_lca_id(&[leaf1.get_id(), leaf2.get_id()]);
            let lca_states = tree_data.internal_character_states.get(&lca_id);
            
            // Calculate distance from root by traversing up the tree and summing branch lengths
            // This matches tree.get_time(lca) in Python (distance from root in Cassiopeia)
            let lca_height = calculate_distance_from_root(tree, lca_id);
            
            if let (Some(leaf1_taxa), Some(leaf2_taxa)) = (leaf1.get_taxa(), leaf2.get_taxa()) {
                let leaf1_states = tree_data.leaf_character_states.get(leaf1_taxa);
                let leaf2_states = tree_data.leaf_character_states.get(leaf2_taxa);
                
                let mut phs = 0;
                let mut detailed_homoplasies = Vec::new();
                
                if let (Some(lca_chars), Some(leaf1_chars), Some(leaf2_chars)) = 
                    (lca_states, leaf1_states, leaf2_states) {
                    
                    // Character states loaded successfully
                    
                    // Python: for character_id in range(k)
                    for char_idx in 0..tree_data.n_characters {
                        let lca_char = lca_chars.get(char_idx).copied().unwrap_or(tree_data.unedited_state);
                        let leaf1_char = leaf1_chars.get(char_idx).copied().unwrap_or(tree_data.missing_state);
                        let leaf2_char = leaf2_chars.get(char_idx).copied().unwrap_or(tree_data.missing_state);
                        
                        // Python homoplasy condition:
                        // lca_states[character_id] == 0
                        // and leaf_1_states[character_id] > 0
                        // and leaf_1_states[character_id] == leaf_2_states[character_id]
                        if lca_char == tree_data.unedited_state
                            && leaf1_char > tree_data.unedited_state
                            && leaf1_char == leaf2_char {
                            phs += 1;
                            detailed_homoplasies.push((char_idx, lca_char, leaf1_char, leaf2_char));
                        }
                    }
                    
                    total_homoplasies += phs;
                    
                    // Log detailed information for first 10 pairs
                    if pair_count < 10 {
                        info!("Pair {}: {:?}-{:?}", pair_count, leaf1_taxa, leaf2_taxa);
                        info!("  LCA node ID: {}, height: {:.6}", lca_id, lca_height);
                        info!("  Homoplasy count: {}", phs);
                        if !detailed_homoplasies.is_empty() {
                            info!("  Homoplasies at chars: {:?}", detailed_homoplasies.iter().map(|(idx, lca, l1, l2)| 
                                format!("{}(lca:{},l1:{},l2:{})", idx, lca, l1, l2)).collect::<Vec<_>>());
                        }
                        if phs == 0 {
                            // Show a few non-homoplasic characters for debugging
                            let sample_chars: Vec<_> = (0..tree_data.n_characters.min(5)).map(|idx| {
                                let lca_char = lca_chars.get(idx).copied().unwrap_or(tree_data.unedited_state);
                                let leaf1_char = leaf1_chars.get(idx).copied().unwrap_or(tree_data.missing_state);
                                let leaf2_char = leaf2_chars.get(idx).copied().unwrap_or(tree_data.missing_state);
                                format!("{}(lca:{},l1:{},l2:{})", idx, lca_char, leaf1_char, leaf2_char)
                            }).collect();
                            info!("  Sample chars (no homoplasies): {:?}", sample_chars);
                        }
                    }
                }
                
                homoplasy_counts.push(phs);
                lca_heights.push(lca_height);
                pair_count += 1;
            }
        }
    }
    
    info!("STELLARS Homoplasy Summary:");
    info!("  Total pairs analyzed: {}", pair_count);
    info!("  Total homoplasies found: {}", total_homoplasies);
    info!("  Average homoplasies per pair: {:.3}", if pair_count > 0 { total_homoplasies as f64 / pair_count as f64 } else { 0.0 });
    
    // Calculate homoplasy distribution for comparison with Python
    let mut homoplasy_dist = std::collections::HashMap::new();
    for &count in &homoplasy_counts {
        *homoplasy_dist.entry(count).or_insert(0) += 1;
    }
    let mut sorted_dist: Vec<_> = homoplasy_dist.iter().collect();
    sorted_dist.sort_by_key(|(k, _)| *k);
    info!("  Homoplasy distribution: {:?}", sorted_dist);
    
    (homoplasy_counts, lca_heights)
}

/// Calculate PHS score using binomial p-values and Benjamini-Hochberg correction
fn calculate_phs_pvalues(
    homoplasy_counts: &[usize],
    lca_heights: &[f64],
    k: usize, // number of characters
    mutation_rate: f64,
    collision_probability: f64,
) -> f64 {
    let n = homoplasy_counts.len();
    let mut pvalues = Vec::with_capacity(n);
    
    info!("=== STELLARS DETAILED P-VALUE CALCULATION ===");
    info!("Starting PHS p-value calculation with n={} pairs, k={} characters", n, k);
    info!("Parameters: mutation_rate={:.6}, collision_probability={:.6}", mutation_rate, collision_probability);
    
    // Log sample of inputs for comparison
    let sample_size = 10.min(n);
    info!("Sample homoplasy counts (first {}): {:?}", sample_size, &homoplasy_counts[..sample_size]);
    info!("Sample LCA heights (first {}): {:?}", sample_size, 
        &lca_heights[..sample_size].iter().map(|h| format!("{:.6}", h)).collect::<Vec<_>>());
    
    let mut ultra_low_count = 0;
    let mut zero_homoplasy_count = 0;
    let mut probability_sum = 0.0;
    
    // Calculate p-values for each pair
    for i in 0..n {
        let lca_height = lca_heights[i];
        let homoplasy_count = homoplasy_counts[i];
        
        if homoplasy_count == 0 {
            zero_homoplasy_count += 1;
        }
        
        // Calculate probability of homoplasy under null model
        let alpha = (-mutation_rate * lca_height).exp();
        let beta = 1.0 - (-mutation_rate * (1.0 - lca_height)).exp();
        let prob = if lca_height >= 1.0 {
            1.0 // Special case when LCA is at the root - matches Python implementation
        } else {
            alpha * beta * beta * collision_probability
        };
        
        probability_sum += prob;
        
        // Calculate binomial p-value
        // p-value = 1 - P(X <= homoplasy_count - 1) where X ~ Binomial(k, prob)
        let pvalue = if homoplasy_count == 0 {
            1.0
        } else {
            1.0 - binomial_cdf(homoplasy_count - 1, k, prob)
        };
        
        // Avoid exactly zero p-values
        let adjusted_pvalue = if pvalue == 0.0 {
            f64::EPSILON
        } else {
            pvalue
        };
        
        if adjusted_pvalue < 1e-10 {
            ultra_low_count += 1;
        }
        
        // Log detailed information for first 10 pairs
        if i < 10 {
            info!("Pair {}: PHS={}, height={:.6}, alpha={:.6}, beta={:.6}, prob={:.6}, p-value={:.6e}", 
                   i, homoplasy_count, lca_height, alpha, beta, prob, adjusted_pvalue);
            
            // Additional detail for debugging probability calculation
            let exp_term1 = -mutation_rate * lca_height;
            let exp_term2 = -mutation_rate * (1.0 - lca_height);
            info!("    exp({:.6}) = {:.6}, exp({:.6}) = {:.6}", 
                   exp_term1, exp_term1.exp(), exp_term2, exp_term2.exp());
            
            if homoplasy_count > 0 {
                let binomial_cdf_val = binomial_cdf(homoplasy_count - 1, k, prob);
                info!("    binomial_cdf({}, {}, {:.6}) = {:.6e}", homoplasy_count - 1, k, prob, binomial_cdf_val);
            }
        }
        
        pvalues.push(adjusted_pvalue);
    }
    
    info!("P-value calculation summary:");
    info!("  Zero homoplasy pairs: {}/{} ({:.1}%)", zero_homoplasy_count, n, 100.0 * zero_homoplasy_count as f64 / n as f64);
    info!("  Ultra-low p-values (<1e-10): {}/{} ({:.1}%)", ultra_low_count, n, 100.0 * ultra_low_count as f64 / n as f64);
    info!("  Average probability per pair: {:.6e}", probability_sum / n as f64);
    
    let min_pvalue = pvalues.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_pvalue = pvalues.iter().cloned().fold(0.0, f64::max);
    info!("  Raw p-value range: {:.2e} to {:.2e}", min_pvalue, max_pvalue);
    
    // Benjamini-Hochberg correction for multiple testing
    pvalues.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    info!("Sorted p-values (first 10): {:?}", 
        pvalues.iter().take(10).map(|p| format!("{:.2e}", p)).collect::<Vec<_>>());
    
    let mut min_adjusted = f64::INFINITY;
    for (i, &pvalue) in pvalues.iter().enumerate() {
        let adjusted = pvalue * (n as f64) / ((i + 1) as f64);
        
        if i < 10 {
            info!("BH step {}: pvalue={:.8e} * {} / {} = {:.8e}", i, pvalue, n, i+1, adjusted);
        }
        
        if adjusted < min_adjusted {
            min_adjusted = adjusted;
            if i < 10 {
                info!("  -> New minimum: {:.8e}", min_adjusted);
            }
        }
    }
    
    info!("Final Benjamini-Hochberg corrected result: {:.10e}", min_adjusted);
    
    if min_adjusted < 1e-10 {
        info!("🚨 ULTRA-LOW STELLARS RESULT: {:.2e} - This matches the simulation issue!", min_adjusted);
        info!("This suggests the issue is in either:");
        info!("  1. Homoplasy overcounting (too many homoplasies detected)");
        info!("  2. Probability miscalculation (prob values too extreme)");
        info!("  3. Binomial CDF numerical issues");
        info!("  4. Benjamini-Hochberg correction amplifying small errors");
    }
    
    min_adjusted
}

/// Binomial CDF calculation using statrs for numerical precision
fn binomial_cdf(k: usize, n: usize, p: f64) -> f64 {
    if p <= 0.0 {
        return 1.0; // All values are <= k
    }
    if p >= 1.0 {
        return if k >= n { 1.0 } else { 0.0 };
    }
    
    // Use statrs for accurate binomial CDF calculation
    match Binomial::new(p, n as u64) {
        Ok(binomial) => binomial.cdf(k as u64),
        Err(_) => {
            // Fallback for edge cases
            if k >= n { 1.0 } else { 0.0 }
        }
    }
}


/// Main PHS calculation function following Python implementation
fn calculate_phs_direct(
    tree: &PhyloTree,
    tree_data: &PHSTreeData,
    mutation_rate: f64,
    collision_probability: f64,
) -> (f64, usize, String, usize) {
    // Calculate homoplasies and LCA heights directly
    let (homoplasy_counts, lca_heights) = calculate_phs_and_lca_heights_direct(tree, tree_data);
    
    let total_pairs = homoplasy_counts.len();
    
    // Remove debug output for production
    
    if total_pairs == 0 {
        return (1.0, 0, "no_pairs".to_string(), 1);
    }
    
    // Calculate PHS score using same logic as Python
    let phs_score = calculate_phs_pvalues(
        &homoplasy_counts, 
        &lca_heights, 
        tree_data.n_characters, 
        mutation_rate, 
        collision_probability
    );
    
    (phs_score, total_pairs, "direct_calculation".to_string(), 1)
}

/// Main optimized PHS function
pub fn optimized_phs(
    tree: &mut PhyloTree,
    character_matrix: Vec<Vec<i32>>, // character_matrix[leaf_idx][character_idx]
    internal_character_states: HashMap<String, Vec<i32>>, // node_name -> character states
    mutation_rate: Option<f64>,
    collision_probability: Option<f64>,
    missing_state: i32,
    unedited_state: i32,
    _max_threads: Option<usize>, // ignored for now
    use_provided_internal_states: Option<bool>, // NEW: use provided states instead of Fitch
    leaf_names: Option<Vec<String>>, // NEW: leaf names corresponding to character_matrix rows
) -> PHSResult {
    let start_time = std::time::Instant::now();
    
    // Build tree data structure
    let use_provided = use_provided_internal_states.unwrap_or(false);
    let tree_data = PHSTreeData::new(
        tree,
        character_matrix.clone(),
        internal_character_states,
        missing_state,
        unedited_state,
        use_provided,
        leaf_names,
    );
    
    // Estimate mutation rate if not provided
    let mut_rate = mutation_rate.unwrap_or_else(|| {
        // Simple estimation: proportion of mutated characters
        let total_chars: usize = character_matrix.iter()
            .map(|states| states.iter().filter(|&&s| s != missing_state).count())
            .sum();
        let mutated_chars: usize = character_matrix.iter()
            .map(|states| states.iter().filter(|&&s| s != missing_state && s != unedited_state).count())
            .sum();
        
        if total_chars > 0 {
            let proportion = mutated_chars as f64 / total_chars as f64;
            -(1.0 - proportion).ln()
        } else {
            0.1 // Default fallback
        }
    });
    
    // Use provided collision probability or default
    let coll_prob = collision_probability.unwrap_or(0.1); // Default collision probability
    
    // Calculate using direct method
    let (phs_score, total_pairs, method, chunks) = 
        calculate_phs_direct(tree, &tree_data, mut_rate, coll_prob);
    
    let computation_time = start_time.elapsed().as_secs_f64() * 1000.0; // Convert to milliseconds
    
    PHSResult {
        phs_score,
        total_pairs,
        computation_time_ms: computation_time,
        method_used: method,
        parallel_chunks_used: chunks,
    }
}
