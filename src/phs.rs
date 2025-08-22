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
    ) -> Self {
        tree.precompute_constant_time_lca();
        
        // Build leaf character states map
        let mut leaf_character_states = FxHashMap::default();
        for (idx, leaf) in tree.get_leaves().enumerate() {
            if let Some(taxa) = leaf.get_taxa() {
                if idx < character_matrix.len() {
                    leaf_character_states.insert(taxa.clone(), character_matrix[idx].clone());
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
            debug!("Node {}: {:?}", node_id, states);
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

/// Calculate homoplasies following the Python implementation exactly
/// Python counts each leaf pair exactly once with their LCA
fn calculate_phs_and_lca_heights_direct(
    tree: &PhyloTree,
    tree_data: &PHSTreeData,
) -> (Vec<usize>, Vec<f64>) {
    let mut homoplasy_counts = Vec::new();
    let mut lca_heights = Vec::new();
    
    // Get maximum depth for normalization
    let max_depth = tree.get_leaves()
        .map(|leaf| tree.depth(leaf.get_id()))
        .max()
        .unwrap_or(0) as f64;
    
    // Get all leaves once
    let all_leaves: Vec<_> = tree.get_leaves().collect();
    
    // Python approach: iterate through all leaf pairs exactly once
    for i in 0..all_leaves.len() {
        for j in i+1..all_leaves.len() {
            let leaf1 = &all_leaves[i];
            let leaf2 = &all_leaves[j];
            
            // Find the LCA of this pair
            let lca_id = tree.get_lca_id(&[leaf1.get_id(), leaf2.get_id()]);
            let lca_states = tree_data.internal_character_states.get(&lca_id);
            
            // Get normalized height for this LCA
            let lca_depth = tree.depth(lca_id) as f64;
            let lca_height = if max_depth > 0.0 { 
                (max_depth - lca_depth) / max_depth 
            } else { 
                0.0 
            };
            
            if let (Some(leaf1_taxa), Some(leaf2_taxa)) = (leaf1.get_taxa(), leaf2.get_taxa()) {
                let leaf1_states = tree_data.leaf_character_states.get(leaf1_taxa);
                let leaf2_states = tree_data.leaf_character_states.get(leaf2_taxa);
                
                let mut phs = 0;
                
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
                        }
                    }
                }
                
                homoplasy_counts.push(phs);
                lca_heights.push(lca_height);
            }
        }
    }
    
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
    
    // Calculate p-values for each pair
    for i in 0..n {
        let lca_height = lca_heights[i];
        let homoplasy_count = homoplasy_counts[i];
        
        // Calculate probability of homoplasy under null model
        let alpha = (-mutation_rate * lca_height).exp();
        let beta = 1.0 - (-mutation_rate * (1.0 - lca_height)).exp();
        let prob = if lca_height >= 1.0 {
            1.0 // Special case when LCA is at the root - matches Python implementation
        } else {
            alpha * beta * beta * collision_probability
        };
        
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
        
        pvalues.push(adjusted_pvalue);
    }
    
    // Benjamini-Hochberg correction for multiple testing
    pvalues.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    let mut min_adjusted = f64::INFINITY;
    for (i, &pvalue) in pvalues.iter().enumerate() {
        let adjusted = pvalue * (n as f64) / ((i + 1) as f64);
        if adjusted < min_adjusted {
            min_adjusted = adjusted;
        }
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
