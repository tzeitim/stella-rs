use phylo::prelude::*;
use phylo::node::NodeID;
use std::collections::HashMap;

/// Likelihood calculation result
#[derive(Debug, Clone)]
pub struct LikelihoodResult {
    pub log_likelihood: f64,
    pub likelihood: f64,
    pub computation_time_ms: f64,
    pub method_used: String,
    pub n_characters: usize,
    pub n_leaves: usize,
}

/// Calculate likelihood for a tree under the CRISPR-Cas9 non-modifiable model
/// Following the likelihood approach mentioned in the PHS paper
pub fn calculate_likelihood(
    tree: &mut PhyloTree,
    character_matrix: Vec<Vec<i32>>, // character_matrix[leaf_idx][character_idx]
    internal_character_states: Option<HashMap<String, Vec<i32>>>, // node_name -> character states
    mutation_rate: Option<f64>,
    collision_probability: Option<f64>,
    missing_state: i32,
    unedited_state: i32,
) -> LikelihoodResult {
    let start_time = std::time::Instant::now();
    
    let n_leaves = tree.get_leaves().count();
    let n_characters = if character_matrix.is_empty() { 0 } else { character_matrix[0].len() };
    
    if n_characters == 0 || n_leaves == 0 {
        return LikelihoodResult {
            log_likelihood: 0.0,
            likelihood: 1.0,
            computation_time_ms: 0.0,
            method_used: "empty_data".to_string(),
            n_characters: 0,
            n_leaves,
        };
    }
    
    // Use default parameters if not provided
    let lambda = mutation_rate.unwrap_or(0.1);
    let q = collision_probability.unwrap_or(0.1);
    
    // If internal states are not provided, use small MP inference
    let final_internal_states = if let Some(states) = internal_character_states {
        states
    } else {
        infer_internal_states_for_likelihood(tree, &character_matrix, missing_state, unedited_state)
    };
    
    // Calculate likelihood using the CRISPR-Cas9 model
    let log_likelihood = calculate_tree_log_likelihood(
        tree, 
        &character_matrix, 
        &final_internal_states, 
        lambda, 
        q,
        missing_state, 
        unedited_state
    );
    
    let computation_time = start_time.elapsed().as_secs_f64() * 1000.0;
    
    LikelihoodResult {
        log_likelihood,
        likelihood: log_likelihood.exp(),
        computation_time_ms: computation_time,
        method_used: "crispr_cas9_model".to_string(),
        n_characters,
        n_leaves,
    }
}

/// Infer internal states for likelihood calculation
/// Similar to parsimony but optimized for likelihood computation
fn infer_internal_states_for_likelihood(
    tree: &mut PhyloTree,
    character_matrix: &[Vec<i32>],
    missing_state: i32,
    unedited_state: i32,
) -> HashMap<String, Vec<i32>> {
    let mut internal_states = HashMap::new();
    
    let k = if character_matrix.is_empty() { 
        return internal_states; 
    } else { 
        character_matrix[0].len() 
    };
    
    // Create taxa to index mapping
    let mut taxa_to_idx = HashMap::new();
    for (idx, leaf) in tree.get_leaves().enumerate() {
        if let Some(taxa) = leaf.get_taxa() {
            taxa_to_idx.insert(taxa.clone(), idx);
        }
    }
    
    // Use maximum likelihood estimation for internal states
    let root = tree.get_root();
    infer_ml_states_recursive(
        tree, 
        root.get_id(), 
        &taxa_to_idx, 
        character_matrix, 
        &mut internal_states, 
        k, 
        missing_state, 
        unedited_state
    );
    
    // Root starts unedited
    if let Some(root_taxa) = root.get_taxa() {
        internal_states.insert(root_taxa.clone(), vec![unedited_state; k]);
    }
    
    internal_states
}

/// Recursive ML state inference
fn infer_ml_states_recursive(
    tree: &PhyloTree,
    node_id: NodeID,
    taxa_to_idx: &HashMap<String, usize>,
    character_matrix: &[Vec<i32>],
    internal_states: &mut HashMap<String, Vec<i32>>,
    k: usize,
    missing_state: i32,
    unedited_state: i32,
) -> Vec<i32> {
    let node = tree.get_node(node_id).unwrap();
    
    // If it's a leaf, return observed states
    if node.is_leaf() {
        if let Some(taxa) = node.get_taxa() {
            if let Some(&leaf_idx) = taxa_to_idx.get(taxa) {
                if leaf_idx < character_matrix.len() {
                    return character_matrix[leaf_idx].clone();
                }
            }
        }
        return vec![missing_state; k];
    }
    
    // For internal nodes, get states from children
    let children: Vec<NodeID> = node.get_children().collect();
    
    if children.len() == 2 {
        let left_states = infer_ml_states_recursive(
            tree, children[0], taxa_to_idx, character_matrix, internal_states, k, missing_state, unedited_state
        );
        let right_states = infer_ml_states_recursive(
            tree, children[1], taxa_to_idx, character_matrix, internal_states, k, missing_state, unedited_state
        );
        
        // Use ML estimation for internal states under non-modifiability
        let mut node_states = vec![unedited_state; k];
        
        for i in 0..k {
            let left_char = if i < left_states.len() { left_states[i] } else { missing_state };
            let right_char = if i < right_states.len() { right_states[i] } else { missing_state };
            
            // Under non-modifiability, if both children have the same mutated state,
            // the parent most likely had that state too (or was unedited)
            // We use the parsimony principle as a proxy for ML estimation
            if left_char == right_char && left_char != missing_state && left_char != unedited_state {
                node_states[i] = left_char; // Parent likely had this mutation
            } else {
                node_states[i] = unedited_state; // Parent was unedited
            }
        }
        
        if let Some(taxa) = node.get_taxa() {
            internal_states.insert(taxa.clone(), node_states.clone());
        }
        
        node_states
    } else {
        vec![unedited_state; k]
    }
}

/// Calculate log-likelihood for the entire tree under CRISPR-Cas9 model
fn calculate_tree_log_likelihood(
    tree: &PhyloTree,
    character_matrix: &[Vec<i32>],
    internal_character_states: &HashMap<String, Vec<i32>>,
    mutation_rate: f64,
    collision_probability: f64,
    missing_state: i32,
    unedited_state: i32,
) -> f64 {
    let mut total_log_likelihood = 0.0;
    
    // Create taxa to index mapping
    let mut taxa_to_idx = HashMap::new();
    for (idx, leaf) in tree.get_leaves().enumerate() {
        if let Some(taxa) = leaf.get_taxa() {
            taxa_to_idx.insert(taxa.clone(), idx);
        }
    }
    
    // Calculate likelihood for each edge in the tree
    for node in tree.get_nodes() {
        if let Some(parent_id) = node.get_parent() {
            let parent_node = tree.get_node(parent_id).unwrap();
            
            // Get branch length (time)
            let branch_length = get_branch_length(&node, parent_node);
            
            // Get states for parent and child
            let parent_states = get_node_states_likelihood(parent_node, &taxa_to_idx, character_matrix, internal_character_states);
            let child_states = get_node_states_likelihood(&node, &taxa_to_idx, character_matrix, internal_character_states);
            
            // Calculate likelihood for this edge
            let edge_log_likelihood = calculate_edge_log_likelihood(
                &parent_states, 
                &child_states, 
                branch_length, 
                mutation_rate, 
                collision_probability,
                missing_state,
                unedited_state
            );
            
            total_log_likelihood += edge_log_likelihood;
        }
    }
    
    total_log_likelihood
}

/// Get branch length between child and parent nodes
fn get_branch_length(_child_node: &Node<String, f32, f32>, _parent_node: &Node<String, f32, f32>) -> f64 {
    // In the absence of explicit branch lengths, use unit length
    // In practice, this would come from the tree structure or be estimated
    1.0
}

/// Get character states for likelihood calculation
fn get_node_states_likelihood(
    node: &Node<String, f32, f32>,
    taxa_to_idx: &HashMap<String, usize>,
    character_matrix: &[Vec<i32>],
    internal_character_states: &HashMap<String, Vec<i32>>,
) -> Vec<i32> {
    if let Some(taxa) = node.get_taxa() {
        // Check if it's a leaf
        if let Some(&leaf_idx) = taxa_to_idx.get(taxa) {
            if leaf_idx < character_matrix.len() {
                return character_matrix[leaf_idx].clone();
            }
        }
        
        // Check if it's an internal node
        if let Some(states) = internal_character_states.get(taxa) {
            return states.clone();
        }
    }
    
    Vec::new()
}

/// Calculate log-likelihood for a single edge under CRISPR-Cas9 model
fn calculate_edge_log_likelihood(
    parent_states: &[i32],
    child_states: &[i32],
    branch_length: f64,
    mutation_rate: f64,
    collision_probability: f64,
    missing_state: i32,
    unedited_state: i32,
) -> f64 {
    let mut log_likelihood = 0.0;
    let k = parent_states.len().min(child_states.len());
    
    for i in 0..k {
        let parent_char = parent_states[i];
        let child_char = child_states[i];
        
        // Skip missing data
        if parent_char == missing_state || child_char == missing_state {
            continue;
        }
        
        // Calculate likelihood for this character position
        let char_log_likelihood = calculate_character_log_likelihood(
            parent_char,
            child_char,
            branch_length,
            mutation_rate,
            collision_probability,
            unedited_state,
        );
        
        log_likelihood += char_log_likelihood;
    }
    
    log_likelihood
}

/// Calculate log-likelihood for a single character under CRISPR-Cas9 model
fn calculate_character_log_likelihood(
    parent_state: i32,
    child_state: i32,
    branch_length: f64,
    mutation_rate: f64,
    collision_probability: f64,
    unedited_state: i32,
) -> f64 {
    // Probability of mutation occurring along this branch
    let mutation_prob = 1.0 - (-mutation_rate * branch_length).exp();
    
    if parent_state == unedited_state {
        if child_state == unedited_state {
            // No mutation occurred
            (1.0 - mutation_prob).ln()
        } else {
            // Mutation occurred to specific state
            (mutation_prob * collision_probability).ln()
        }
    } else {
        // Parent already mutated - under non-modifiability, no further change
        if child_state == parent_state {
            0.0 // Probability 1 (log = 0)
        } else {
            // This shouldn't happen under strict non-modifiability
            f64::NEG_INFINITY
        }
    }
}

/// Calculate likelihood-based distance as mentioned in the paper
/// LD(T,S) = max{0, log(L(T_GT)) / log(L(T,S)) - 1}
pub fn calculate_likelihood_distance(
    reconstructed_tree: &mut PhyloTree,
    ground_truth_tree: &mut PhyloTree,
    character_matrix: Vec<Vec<i32>>,
    internal_states_reconstructed: Option<HashMap<String, Vec<i32>>>,
    internal_states_ground_truth: Option<HashMap<String, Vec<i32>>>,
    mutation_rate: Option<f64>,
    collision_probability: Option<f64>,
    missing_state: i32,
    unedited_state: i32,
) -> f64 {
    // Calculate likelihood for reconstructed tree
    let reconstructed_likelihood = calculate_likelihood(
        reconstructed_tree,
        character_matrix.clone(),
        internal_states_reconstructed,
        mutation_rate,
        collision_probability,
        missing_state,
        unedited_state,
    );
    
    // Calculate likelihood for ground truth tree  
    let ground_truth_likelihood = calculate_likelihood(
        ground_truth_tree,
        character_matrix,
        internal_states_ground_truth,
        mutation_rate,
        collision_probability,
        missing_state,
        unedited_state,
    );
    
    // Calculate likelihood distance: max{0, log(L_GT) / log(L_recon) - 1}
    if reconstructed_likelihood.log_likelihood == 0.0 {
        return f64::INFINITY;
    }
    
    let ratio = ground_truth_likelihood.log_likelihood / reconstructed_likelihood.log_likelihood;
    (ratio - 1.0).max(0.0)
}