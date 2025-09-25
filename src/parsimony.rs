use phylo::prelude::*;
use phylo::node::NodeID;
use std::collections::{HashMap, HashSet};
use fxhash::FxHashMap;

/// Parsimony calculation result
#[derive(Debug, Clone)]
pub struct ParsimonyResult {
    pub parsimony_score: usize,
    pub total_mutations: usize,
    pub computation_time_ms: f64,
    pub method_used: String,
    pub internal_states_inferred: bool,
}

/// Result for parsimony p-value calculation
#[derive(Debug, Clone)]
pub struct ParsimonyPValueResult {
    pub p_value: f64,
    pub parsimony_score: usize,
    pub computation_time_ms: f64,
    pub method_used: String,
}

/// Calculate parsimony score for a given tree and character matrix
/// Following the small MP problem from the PHS paper (Algorithm 1)
pub fn calculate_parsimony(
    tree: &mut PhyloTree,
    character_matrix: Vec<Vec<i32>>, // character_matrix[leaf_idx][character_idx]
    internal_character_states: Option<HashMap<String, Vec<i32>>>, // node_name -> character states
    missing_state: i32,
    unedited_state: i32,
) -> ParsimonyResult {
    let start_time = std::time::Instant::now();
    
    // If internal states are not provided, infer them using the small MP algorithm
    let (final_internal_states, inferred) = if let Some(states) = internal_character_states {
        (states, false)
    } else {
        let inferred_states = infer_internal_states_small_mp(tree, &character_matrix, missing_state, unedited_state);
        (inferred_states, true)
    };
    
    // Calculate total mutations across the tree
    let total_mutations = count_total_mutations(tree, &character_matrix, &final_internal_states, missing_state, unedited_state);
    
    let computation_time = start_time.elapsed().as_secs_f64() * 1000.0;
    
    ParsimonyResult {
        parsimony_score: total_mutations,
        total_mutations,
        computation_time_ms: computation_time,
        method_used: if inferred { "small_mp_inferred".to_string() } else { "provided_states".to_string() },
        internal_states_inferred: inferred,
    }
}

/// Implementation of the small MP algorithm from the PHS paper (Algorithm 1)
/// This reconstructs ancestral states to minimize mutations under non-modifiability
pub fn infer_internal_states_small_mp(
    tree: &mut PhyloTree,
    character_matrix: &[Vec<i32>],
    missing_state: i32,
    unedited_state: i32,
) -> HashMap<String, Vec<i32>> {
    let mut internal_states = HashMap::new();
    
    // Get the number of characters
    let k = if character_matrix.is_empty() { 
        return internal_states; 
    } else { 
        character_matrix[0].len() 
    };
    
    // Create taxa to index mapping for leaf states
    let mut taxa_to_idx = HashMap::new();
    for (idx, leaf) in tree.get_leaves().enumerate() {
        if let Some(taxa) = leaf.get_taxa() {
            taxa_to_idx.insert(taxa.clone(), idx);
        }
    }
    
    // Perform post-order traversal to infer internal states
    let root = tree.get_root();
    infer_states_recursive(
        tree, 
        root.get_id(), 
        &taxa_to_idx, 
        character_matrix, 
        &mut internal_states, 
        k, 
        missing_state, 
        unedited_state
    );
    
    // Set root to all unedited (as per paper)
    if let Some(root_taxa) = root.get_taxa() {
        internal_states.insert(root_taxa.clone(), vec![unedited_state; k]);
    }
    
    internal_states
}

/// Recursive function for post-order state inference (Algorithm 1 from paper)
fn infer_states_recursive(
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
    
    // If it's a leaf, return its observed states
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
        // Process left and right children
        let left_states = infer_states_recursive(
            tree, children[0], taxa_to_idx, character_matrix, internal_states, k, missing_state, unedited_state
        );
        let right_states = infer_states_recursive(
            tree, children[1], taxa_to_idx, character_matrix, internal_states, k, missing_state, unedited_state
        );
        
        // Infer this node's states according to Algorithm 1
        let mut node_states = vec![unedited_state; k];
        
        for i in 0..k {
            let left_char = if i < left_states.len() { left_states[i] } else { missing_state };
            let right_char = if i < right_states.len() { right_states[i] } else { missing_state };
            
            // Non-modifiability logic from Algorithm 1:
            // If left == 0 or right == 0 or left != right, then parent = 0
            // Else parent = left (== right)
            if left_char == unedited_state || right_char == unedited_state || left_char != right_char {
                node_states[i] = unedited_state;
            } else if left_char != missing_state && right_char != missing_state {
                // Both children have the same non-missing, non-unedited state
                node_states[i] = left_char;
            } else {
                node_states[i] = unedited_state;
            }
        }
        
        // Store the inferred states
        if let Some(taxa) = node.get_taxa() {
            internal_states.insert(taxa.clone(), node_states.clone());
        }
        
        node_states
    } else {
        // Handle non-binary trees (default to unedited)
        vec![unedited_state; k]
    }
}

/// Count total mutations in the tree given character states
fn count_total_mutations(
    tree: &PhyloTree,
    character_matrix: &[Vec<i32>],
    internal_character_states: &HashMap<String, Vec<i32>>,
    missing_state: i32,
    unedited_state: i32,
) -> usize {
    let mut total_mutations = 0;
    
    // Create taxa to index mapping for leaves
    let mut taxa_to_idx = HashMap::new();
    for (idx, leaf) in tree.get_leaves().enumerate() {
        if let Some(taxa) = leaf.get_taxa() {
            taxa_to_idx.insert(taxa.clone(), idx);
        }
    }
    
    // Count mutations along each edge
    for node in tree.get_nodes() {
        if let Some(parent_id) = node.get_parent() {
            let parent_node = tree.get_node(parent_id).unwrap();
            
            // Get states for parent and child
            let parent_states = get_node_states(parent_node, &taxa_to_idx, character_matrix, internal_character_states);
            let child_states = get_node_states(&node, &taxa_to_idx, character_matrix, internal_character_states);
            
            // Count mutations between parent and child
            let k = parent_states.len().min(child_states.len());
            for i in 0..k {
                let parent_char = parent_states[i];
                let child_char = child_states[i];
                
                // A mutation occurs when:
                // 1. Parent is unedited (0) and child is not unedited and not missing
                // 2. Under non-modifiability, once mutated, sites cannot change again
                if parent_char == unedited_state && child_char != unedited_state && child_char != missing_state {
                    total_mutations += 1;
                }
            }
        }
    }
    
    total_mutations
}

/// Get character states for a node (leaf or internal)
fn get_node_states(
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
        
        // Check if it's an internal node with known states
        if let Some(states) = internal_character_states.get(taxa) {
            return states.clone();
        }
    }
    
    // Default: return empty vector
    Vec::new()
}

/// Calculate parsimony p-value according to the paper's conditional distribution
/// This implements the parsimony tail probability mentioned in the paper
pub fn calculate_parsimony_p_value(
    tree: &mut PhyloTree,
    character_matrix: Vec<Vec<i32>>,
    internal_character_states: Option<HashMap<String, Vec<i32>>>,
    mutation_rate: Option<f64>,
    missing_state: i32,
    unedited_state: i32,
) -> ParsimonyPValueResult {
    let start_time = std::time::Instant::now();
    
    // First calculate the parsimony score
    let parsimony_result = calculate_parsimony(
        tree, 
        character_matrix.clone(), 
        internal_character_states, 
        missing_state, 
        unedited_state
    );
    
    // Calculate p-value based on conditional parsimony distribution
    // This is a simplified implementation - the paper refers to analytical calculation
    // For now, we use an approximation based on Poisson distribution
    let lambda = mutation_rate.unwrap_or_else(|| {
        eprintln!("WARNING: Using default mutation_rate=0.1 in parsimony! Pass actual calculated rate!");
        0.1
    });
    let n_leaves = tree.get_leaves().count();
    let expected_mutations = lambda * (n_leaves as f64 - 1.0); // Expected mutations in tree
    
    // Calculate p-value using Poisson approximation
    // P(M >= observed_parsimony | tree_topology)
    let p_value = poisson_tail_probability(expected_mutations, parsimony_result.parsimony_score);
    
    let computation_time = start_time.elapsed().as_secs_f64() * 1000.0;
    
    ParsimonyPValueResult {
        p_value,
        parsimony_score: parsimony_result.parsimony_score,
        computation_time_ms: computation_time,
        method_used: "poisson_approximation".to_string(),
    }
}

/// Calculate Poisson tail probability P(X >= k) where X ~ Poisson(lambda)
fn poisson_tail_probability(lambda: f64, k: usize) -> f64 {
    if lambda <= 0.0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    
    // For large k, use normal approximation
    if k > 100 {
        let z = (k as f64 - lambda) / lambda.sqrt();
        return 0.5 * (1.0 - erf(z / 2_f64.sqrt()));
    }
    
    // Calculate exact Poisson tail probability
    let _tail_prob = 0.0;
    let mut poisson_pmf = (-lambda).exp(); // P(X = 0)
    
    // Sum from 0 to k-1, then subtract from 1
    let mut cumulative_prob = poisson_pmf; // P(X = 0)
    
    for i in 1..k {
        poisson_pmf *= lambda / (i as f64);
        cumulative_prob += poisson_pmf;
    }
    
    1.0 - cumulative_prob
}

/// Error function approximation for normal distribution
fn erf(x: f64) -> f64 {
    // Abramowitz and Stegun approximation
    let a1 =  0.254829592;
    let a2 = -0.284496736;
    let a3 =  1.421413741;
    let a4 = -1.453152027;
    let a5 =  1.061405429;
    let p  =  0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();

    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

    sign * y
}

/// Fitch parsimony algorithm for reconstructing ancestral character states
/// This matches the algorithm used by Cassiopeia and other phylogenetic libraries
/// Returns node ID -> character states mapping for internal nodes
pub fn fitch_parsimony(
    tree: &mut PhyloTree,
    character_matrix: &[Vec<i32>],
    missing_state: i32,
    unedited_state: i32,
) -> FxHashMap<NodeID, Vec<i32>> {
    let mut ancestral_states = FxHashMap::default();
    
    // Get number of characters
    let k = if character_matrix.is_empty() { 
        return ancestral_states; 
    } else { 
        character_matrix[0].len() 
    };
    
    // Create taxa to leaf index mapping
    let mut taxa_to_idx = HashMap::new();
    for (idx, leaf) in tree.get_leaves().enumerate() {
        if let Some(taxa) = leaf.get_taxa() {
            taxa_to_idx.insert(taxa.clone(), idx);
        }
    }
    
    // Phase 1: Postorder pass - compute possible states for each internal node
    let root_id = tree.get_root_id();
    let mut possible_states: HashMap<NodeID, Vec<HashSet<i32>>> = HashMap::new();
    
    fitch_postorder_pass(
        tree, 
        root_id, 
        &taxa_to_idx, 
        character_matrix, 
        &mut possible_states, 
        k, 
        missing_state, 
        unedited_state
    );
    
    // Phase 2: Preorder pass - choose optimal states for each internal node
    fitch_preorder_pass(
        tree,
        root_id,
        None, // No parent for root
        &possible_states,
        &mut ancestral_states,
        k,
        unedited_state
    );
    
    ancestral_states
}

/// Postorder pass of Fitch algorithm - compute possible states
fn fitch_postorder_pass(
    tree: &PhyloTree,
    node_id: NodeID,
    taxa_to_idx: &HashMap<String, usize>,
    character_matrix: &[Vec<i32>],
    possible_states: &mut HashMap<NodeID, Vec<HashSet<i32>>>,
    k: usize,
    missing_state: i32,
    unedited_state: i32,
) -> Vec<HashSet<i32>> {
    let node = tree.get_node(node_id).unwrap();
    
    // If it's a leaf, return its observed states
    if node.is_leaf() {
        let mut node_possible_states = vec![HashSet::new(); k];
        
        if let Some(taxa) = node.get_taxa() {
            if let Some(&leaf_idx) = taxa_to_idx.get(taxa) {
                if leaf_idx < character_matrix.len() {
                    let leaf_states = &character_matrix[leaf_idx];
                    for (char_idx, &state) in leaf_states.iter().enumerate() {
                        if char_idx < k {
                            if state == missing_state {
                                // Missing state: all possible states are allowed
                                node_possible_states[char_idx].insert(unedited_state);
                                // Add some common non-zero states
                                for s in 1..=20 {
                                    node_possible_states[char_idx].insert(s);
                                }
                            } else {
                                // Observed state
                                node_possible_states[char_idx].insert(state);
                            }
                        }
                    }
                }
            }
        }
        
        // Ensure each character has at least the unedited state
        for char_states in &mut node_possible_states {
            if char_states.is_empty() {
                char_states.insert(unedited_state);
            }
        }
        
        return node_possible_states;
    }
    
    // For internal nodes, process children first (postorder)
    let children: Vec<NodeID> = node.get_children().collect();
    let mut child_possible_states = Vec::new();
    
    for child_id in children {
        let child_states = fitch_postorder_pass(
            tree,
            child_id,
            taxa_to_idx,
            character_matrix,
            possible_states,
            k,
            missing_state,
            unedited_state,
        );
        child_possible_states.push(child_states);
    }
    
    // Compute possible states for this internal node
    let mut node_possible_states = vec![HashSet::new(); k];
    
    for char_idx in 0..k {
        if child_possible_states.is_empty() {
            // No children - should not happen in binary tree
            node_possible_states[char_idx].insert(unedited_state);
        } else if child_possible_states.len() == 1 {
            // Single child - inherit its states
            node_possible_states[char_idx] = child_possible_states[0][char_idx].clone();
        } else {
            // Multiple children - find intersection and union
            let mut intersection = child_possible_states[0][char_idx].clone();
            let mut union = child_possible_states[0][char_idx].clone();
            
            for child_states in &child_possible_states[1..] {
                let child_char_states = &child_states[char_idx];
                intersection = intersection.intersection(child_char_states).copied().collect();
                union = union.union(child_char_states).copied().collect();
            }
            
            // Fitch rule: if intersection is non-empty, use intersection; otherwise use union
            if !intersection.is_empty() {
                node_possible_states[char_idx] = intersection;
            } else {
                node_possible_states[char_idx] = union;
            }
        }
        
        // Ensure we have at least one state
        if node_possible_states[char_idx].is_empty() {
            node_possible_states[char_idx].insert(unedited_state);
        }
    }
    
    // Store the possible states for this node
    possible_states.insert(node_id, node_possible_states.clone());
    
    node_possible_states
}

/// Preorder pass of Fitch algorithm - choose optimal states
fn fitch_preorder_pass(
    tree: &PhyloTree,
    node_id: NodeID,
    parent_chosen_states: Option<&Vec<i32>>,
    possible_states: &HashMap<NodeID, Vec<HashSet<i32>>>,
    ancestral_states: &mut FxHashMap<NodeID, Vec<i32>>,
    k: usize,
    unedited_state: i32,
) {
    let node = tree.get_node(node_id).unwrap();
    
    // Skip leaves - they already have their states
    if node.is_leaf() {
        return;
    }
    
    // Get possible states for this node
    let node_possible_states = possible_states.get(&node_id);
    if node_possible_states.is_none() {
        // Fallback to unedited states
        ancestral_states.insert(node_id, vec![unedited_state; k]);
        return;
    }
    
    let node_possible_states = node_possible_states.unwrap();
    let mut chosen_states = vec![unedited_state; k];
    
    for char_idx in 0..k {
        let possible_char_states = &node_possible_states[char_idx];
        
        if possible_char_states.is_empty() {
            chosen_states[char_idx] = unedited_state;
        } else if possible_char_states.len() == 1 {
            // Only one possible state
            chosen_states[char_idx] = *possible_char_states.iter().next().unwrap();
        } else {
            // Multiple possible states - choose based on parent if available
            if let Some(parent_states) = parent_chosen_states {
                if char_idx < parent_states.len() {
                    let parent_state = parent_states[char_idx];
                    if possible_char_states.contains(&parent_state) {
                        // Choose parent's state if possible
                        chosen_states[char_idx] = parent_state;
                    } else {
                        // Choose the unedited state if possible, otherwise any state
                        if possible_char_states.contains(&unedited_state) {
                            chosen_states[char_idx] = unedited_state;
                        } else {
                            chosen_states[char_idx] = *possible_char_states.iter().next().unwrap();
                        }
                    }
                } else {
                    // Choose unedited state if possible
                    if possible_char_states.contains(&unedited_state) {
                        chosen_states[char_idx] = unedited_state;
                    } else {
                        chosen_states[char_idx] = *possible_char_states.iter().next().unwrap();
                    }
                }
            } else {
                // No parent (this is root) - prefer unedited state
                if possible_char_states.contains(&unedited_state) {
                    chosen_states[char_idx] = unedited_state;
                } else {
                    chosen_states[char_idx] = *possible_char_states.iter().next().unwrap();
                }
            }
        }
    }
    
    // Store the chosen states for this internal node
    ancestral_states.insert(node_id, chosen_states.clone());
    
    // Recursively process children
    let children: Vec<NodeID> = node.get_children().collect();
    for child_id in children {
        fitch_preorder_pass(
            tree,
            child_id,
            Some(&chosen_states),
            possible_states,
            ancestral_states,
            k,
            unedited_state,
        );
    }
}