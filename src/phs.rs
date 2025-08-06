use phylo::prelude::*;
use phylo::node::NodeID;
use rayon::prelude::*;
use fxhash::FxHashMap;
use std::collections::HashMap;

/// PHS calculation result using optimized techniques
#[derive(Debug, Clone)]
pub struct PHSResult {
    pub phs_score: f64,
    pub total_pairs: usize,
    pub computation_time_ms: f64,
    pub method_used: String,
    pub parallel_chunks_used: usize,
}

/// Optimized tree data structure for PHS calculation
/// Inspired by phylo-rs and optimized triplets implementation
struct OptimizedPHSTreeData {
    // Fast taxa lookups using FxHashMap (20-30% faster than std HashMap)
    taxa_to_idx: FxHashMap<String, usize>,
    idx_to_taxa: Vec<String>,
    n_leaves: usize,
    
    // Precomputed LCA node IDs in flat 2D array for O(1) access
    // This is the biggest performance win from phylo-rs analysis
    lca_node_ids: Vec<Vec<NodeID>>,
    
    // Precomputed depths for LCA nodes
    lca_depths: Vec<Vec<f64>>,
    
    // Character states for each leaf (indexed by taxa index)
    character_states: Vec<Vec<i32>>,
    
    // Character states for internal nodes (indexed by NodeID)
    internal_character_states: FxHashMap<NodeID, Vec<i32>>,
    
    // Tree reference for getting node character states
    leaf_node_ids: Vec<NodeID>,
    
    // Parameters
    missing_state: i32,
    unedited_state: i32,
}

impl OptimizedPHSTreeData {
    fn new(
        tree: &mut PhyloTree,
        character_states: Vec<Vec<i32>>,
        internal_character_states: HashMap<String, Vec<i32>>,
        missing_state: i32,
        unedited_state: i32,
    ) -> Self {
        tree.precompute_constant_time_lca();
        
        // Use FxHashMap for faster hash operations (phylo-rs technique)
        let mut taxa_to_idx = FxHashMap::default();
        let mut idx_to_taxa = Vec::new();
        let mut leaf_node_ids = Vec::new();
        
        // Build taxa mappings and collect leaf node IDs
        for (idx, leaf) in tree.get_leaves().enumerate() {
            if let Some(taxa) = leaf.get_taxa() {
                taxa_to_idx.insert(taxa.clone(), idx);
                idx_to_taxa.push(taxa.clone());
                leaf_node_ids.push(leaf.get_id());
            }
        }
        
        let n_leaves = idx_to_taxa.len();
        
        // Precompute ALL pairwise LCA node IDs and depths for O(1) access
        // This is inspired by phylo-rs's bipartition preprocessing
        let mut lca_node_ids = vec![vec![leaf_node_ids[0]; n_leaves]; n_leaves];
        let mut lca_depths = vec![vec![0.0; n_leaves]; n_leaves];
        
        // Batch compute all LCA node IDs and depths using efficient tree traversal
        Self::precompute_all_lca_data(&mut lca_node_ids, &mut lca_depths, tree, &leaf_node_ids, n_leaves);
        
        // Convert internal character states to use NodeID keys
        let mut internal_char_states = FxHashMap::default();
        for (node_name, states) in internal_character_states {
            // Find node by name (this could be optimized further if needed)
            for node in tree.get_nodes() {
                if let Some(taxa) = node.get_taxa() {
                    if taxa == &node_name {
                        internal_char_states.insert(node.get_id(), states);
                        break;
                    }
                }
            }
        }
        
        OptimizedPHSTreeData {
            taxa_to_idx,
            idx_to_taxa,
            n_leaves,
            lca_node_ids,
            lca_depths,
            character_states,
            internal_character_states: internal_char_states,
            leaf_node_ids,
            missing_state,
            unedited_state,
        }
    }
    
    /// Optimized batch LCA computation inspired by phylo-rs RF implementation
    fn precompute_all_lca_data(
        lca_node_ids: &mut Vec<Vec<NodeID>>,
        lca_depths: &mut Vec<Vec<f64>>,
        tree: &PhyloTree,
        leaf_ids: &[NodeID],
        n_leaves: usize,
    ) {
        // Diagonal elements (self-to-self) are the leaves themselves
        for i in 0..n_leaves {
            lca_node_ids[i][i] = leaf_ids[i];
            lca_depths[i][i] = tree.depth(leaf_ids[i]) as f64;
        }
        
        // Compute upper triangle, then mirror to lower triangle
        // This reduces computation by half
        for i in 0..n_leaves {
            for j in i+1..n_leaves {
                let lca_id = tree.get_lca_id(&[leaf_ids[i], leaf_ids[j]]);
                let depth = tree.depth(lca_id) as f64;
                lca_node_ids[i][j] = lca_id;
                lca_node_ids[j][i] = lca_id; // Mirror for symmetry
                lca_depths[i][j] = depth;
                lca_depths[j][i] = depth; // Mirror for symmetry
            }
        }
    }
    
    /// Ultra-fast homoplasy counting using precomputed data
    #[inline(always)]
    fn count_homoplasies_optimized(&self, i_idx: usize, j_idx: usize) -> usize {
        // Direct array access - no function calls, maximum performance
        let lca_id = unsafe { *self.lca_node_ids.get_unchecked(i_idx).get_unchecked(j_idx) };
        
        // Get character states
        let leaf1_states = unsafe { self.character_states.get_unchecked(i_idx) };
        let leaf2_states = unsafe { self.character_states.get_unchecked(j_idx) };
        
        // Get LCA character states (if available)
        let lca_states = self.internal_character_states.get(&lca_id);
        
        let mut homoplasy_count = 0;
        
        if let Some(lca_chars) = lca_states {
            // Hot loop: minimize memory accesses and branches
            let k = leaf1_states.len().min(leaf2_states.len()).min(lca_chars.len());
            
            for idx in 0..k {
                let lca_char = unsafe { *lca_chars.get_unchecked(idx) };
                let leaf1_char = unsafe { *leaf1_states.get_unchecked(idx) };
                let leaf2_char = unsafe { *leaf2_states.get_unchecked(idx) };
                
                // Branchless homoplasy detection
                // Homoplasy occurs when:
                // 1. LCA has unedited state
                // 2. Both leaves have same non-missing, non-unedited state
                // 3. The leaves' states are equal
                let lca_is_unedited = lca_char == self.unedited_state;
                let leaf1_is_valid = leaf1_char != self.missing_state && leaf1_char != self.unedited_state;
                let leaf2_is_valid = leaf2_char != self.missing_state && leaf2_char != self.unedited_state;
                let states_equal = leaf1_char == leaf2_char;
                
                // Use arithmetic instead of branches for better CPU prediction
                homoplasy_count += (lca_is_unedited && leaf1_is_valid && leaf2_is_valid && states_equal) as usize;
            }
        }
        
        homoplasy_count
    }
    
    /// Get LCA depth for pair (optimized)
    #[inline(always)]
    fn get_lca_depth_optimized(&self, i_idx: usize, j_idx: usize) -> f64 {
        unsafe { *self.lca_depths.get_unchecked(i_idx).get_unchecked(j_idx) }
    }
}

/// Generate all leaf pairs efficiently using arena-style allocation
fn generate_leaf_pairs_optimized(n: usize) -> Vec<(usize, usize)> {
    let total_pairs = if n < 2 { 0 } else { n * (n - 1) / 2 };
    let mut pairs = Vec::with_capacity(total_pairs); // Pre-allocate exact size
    
    // Generate upper triangle pairs
    for i in 0..n {
        for j in i+1..n {
            pairs.push((i, j));
        }
    }
    pairs
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
        let prob = if lca_height == 1.0 {
            1.0 // Special case when LCA is at the root
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

/// Simple binomial CDF calculation (could be replaced with a more optimized version)
fn binomial_cdf(k: usize, n: usize, p: f64) -> f64 {
    if p == 0.0 {
        return 1.0; // For usize, k is always >= 0
    }
    if p == 1.0 {
        return if k >= n { 1.0 } else { 0.0 };
    }
    
    let mut cdf = 0.0;
    for i in 0..=k {
        cdf += binomial_pmf(i, n, p);
    }
    cdf
}

/// Simple binomial PMF calculation
fn binomial_pmf(k: usize, n: usize, p: f64) -> f64 {
    if k > n {
        return 0.0;
    }
    
    let log_binom_coeff = log_binomial_coefficient(n, k);
    let log_prob = (k as f64) * p.ln() + ((n - k) as f64) * (1.0 - p).ln();
    (log_binom_coeff + log_prob).exp()
}

/// Log binomial coefficient calculation
fn log_binomial_coefficient(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }
    
    let k = k.min(n - k); // Take advantage of symmetry
    let mut result = 0.0;
    for i in 0..k {
        result += (n - i) as f64 / (i + 1) as f64;
    }
    result.ln()
}

/// Optimized parallel PHS calculation using phylo-rs techniques
fn calculate_phs_optimized(
    tree_data: &OptimizedPHSTreeData,
    mutation_rate: f64,
    collision_probability: f64,
    max_threads: Option<usize>,
) -> (f64, usize, String, usize) {
    let n_leaves = tree_data.n_leaves;
    let total_pairs = if n_leaves < 2 { 0 } else { n_leaves * (n_leaves - 1) / 2 };
    
    if total_pairs == 0 {
        return (1.0, 0, "no_pairs".to_string(), 1);
    }
    
    // Use optimized thread calculation
    let target_threads = calculate_optimal_threads(total_pairs, max_threads);
    
    // For small workloads, use sequential with optimized data structures
    if total_pairs < 1000 || target_threads == 1 {
        return calculate_phs_sequential_optimized(tree_data, mutation_rate, collision_probability);
    }
    
    // Generate pairs using optimized allocation
    let all_pairs = generate_leaf_pairs_optimized(n_leaves);
    
    // Advanced chunking strategy inspired by phylo-rs work distribution
    let chunk_size = calculate_optimal_chunk_size(all_pairs.len(), target_threads);
    
    // Parallel processing with optimized data access
    let (homoplasy_counts, lca_heights): (Vec<_>, Vec<_>) = all_pairs
        .par_chunks(chunk_size)
        .flat_map(|pair_chunk| {
            // Hot loop optimization - minimize function calls
            pair_chunk.iter().map(|&(i, j)| {
                let homoplasy_count = tree_data.count_homoplasies_optimized(i, j);
                let lca_height = tree_data.get_lca_depth_optimized(i, j);
                (homoplasy_count, lca_height)
            }).collect::<Vec<_>>()
        })
        .unzip();
    
    // Calculate PHS score
    let k = if tree_data.character_states.is_empty() { 
        0 
    } else { 
        tree_data.character_states[0].len() 
    };
    
    let phs_score = calculate_phs_pvalues(&homoplasy_counts, &lca_heights, k, mutation_rate, collision_probability);
    
    (phs_score, total_pairs, "optimized_parallel".to_string(), target_threads)
}

/// Sequential optimized calculation for small workloads
fn calculate_phs_sequential_optimized(
    tree_data: &OptimizedPHSTreeData,
    mutation_rate: f64,
    collision_probability: f64,
) -> (f64, usize, String, usize) {
    let n_leaves = tree_data.n_leaves;
    let total_pairs = if n_leaves < 2 { 0 } else { n_leaves * (n_leaves - 1) / 2 };
    
    if total_pairs == 0 {
        return (1.0, 0, "no_pairs".to_string(), 1);
    }
    
    let mut homoplasy_counts = Vec::with_capacity(total_pairs);
    let mut lca_heights = Vec::with_capacity(total_pairs);
    
    // Double nested loop with maximum optimization
    for i in 0..n_leaves {
        for j in i+1..n_leaves {
            let homoplasy_count = tree_data.count_homoplasies_optimized(i, j);
            let lca_height = tree_data.get_lca_depth_optimized(i, j);
            homoplasy_counts.push(homoplasy_count);
            lca_heights.push(lca_height);
        }
    }
    
    // Calculate PHS score
    let k = if tree_data.character_states.is_empty() { 
        0 
    } else { 
        tree_data.character_states[0].len() 
    };
    
    let phs_score = calculate_phs_pvalues(&homoplasy_counts, &lca_heights, k, mutation_rate, collision_probability);
    
    (phs_score, total_pairs, "optimized_sequential".to_string(), 1)
}

/// Optimized thread calculation based on workload characteristics
fn calculate_optimal_threads(total_work: usize, max_threads: Option<usize>) -> usize {
    let mut available_threads = rayon::current_num_threads();
    
    if let Some(max) = max_threads {
        available_threads = std::cmp::min(available_threads, max);
    }
    
    // Refined work-per-thread based on phylo-rs analysis
    let optimal_work_per_thread = 2000..8000; // Adjusted for PHS computation characteristics
    
    if total_work < optimal_work_per_thread.start {
        return 1;
    }
    
    let conservative_threads = total_work / optimal_work_per_thread.end;
    let aggressive_threads = (total_work / optimal_work_per_thread.start).max(1);
    
    // Geometric mean for balanced approach
    let target_threads = ((conservative_threads * aggressive_threads) as f64).sqrt() as usize;
    
    std::cmp::min(std::cmp::max(target_threads, 1), available_threads)
}

/// Calculate optimal chunk size for load balancing
fn calculate_optimal_chunk_size(total_items: usize, target_threads: usize) -> usize {
    let base_chunk_size = total_items / target_threads;
    
    // Adaptive chunking based on thread count (phylo-rs inspired)
    let chunk_multiplier = if target_threads >= 32 {
        4 // Many threads - create many small chunks
    } else if target_threads >= 8 {
        3 // Moderate threads
    } else {
        2 // Few threads
    };
    
    std::cmp::max(1, base_chunk_size / chunk_multiplier)
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
    max_threads: Option<usize>,
) -> PHSResult {
    let start_time = std::time::Instant::now();
    
    // Build optimized tree data structure
    let tree_data = OptimizedPHSTreeData::new(
        tree,
        character_matrix,
        internal_character_states,
        missing_state,
        unedited_state,
    );
    
    // Estimate mutation rate if not provided
    let mut_rate = mutation_rate.unwrap_or_else(|| {
        // Simple estimation: proportion of mutated characters
        let total_chars: usize = tree_data.character_states.iter()
            .map(|states| states.iter().filter(|&&s| s != missing_state).count())
            .sum();
        let mutated_chars: usize = tree_data.character_states.iter()
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
    
    // Calculate using optimized functions
    let (phs_score, total_pairs, method, chunks) = 
        calculate_phs_optimized(&tree_data, mut_rate, coll_prob, max_threads);
    
    let computation_time = start_time.elapsed().as_secs_f64() * 1000.0; // Convert to milliseconds
    
    PHSResult {
        phs_score,
        total_pairs,
        computation_time_ms: computation_time,
        method_used: method,
        parallel_chunks_used: chunks,
    }
}