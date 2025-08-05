use phylo::prelude::*;
use rand::prelude::*;
use std::collections::HashMap;

/// Ultra-fast triplet result with additional optimizations
#[derive(Debug, Clone)]
pub struct UltraFastTripletResult {
    pub all_triplets_correct: HashMap<usize, f64>,
    pub resolvable_triplets_correct: HashMap<usize, f64>,
    pub unresolved_triplets_correct: HashMap<usize, f64>,
    pub proportion_unresolvable: HashMap<usize, f64>,
    pub method_used: String, // "exhaustive" or "sampling"
    pub total_triplets_checked: usize,
}

/// Ultra-optimized tree data with flat arrays for maximum speed
struct UltraFastTreeData {
    taxa_to_idx: HashMap<String, usize>,
    idx_to_taxa: Vec<String>,
    n_taxa: usize,
    // Flat 2D array for O(1) LCA depth access: [i][j] = LCA depth of taxa i and j
    lca_depths: Vec<Vec<usize>>,
}

impl UltraFastTreeData {
    fn new(tree: &mut PhyloTree) -> Self {
        tree.precompute_constant_time_lca();
        
        let mut taxa_to_idx = HashMap::new();
        let mut idx_to_taxa = Vec::new();
        
        // Build flat indexing
        for (idx, leaf) in tree.get_leaves().enumerate() {
            if let Some(taxa) = leaf.get_taxa() {
                taxa_to_idx.insert(taxa.clone(), idx);
                idx_to_taxa.push(taxa.clone());
            }
        }
        
        let n_taxa = idx_to_taxa.len();
        let mut lca_depths = vec![vec![0; n_taxa]; n_taxa];
        
        // Precompute all pairwise LCA depths with flat array access
        let leaf_ids: Vec<_> = tree.get_leaves().map(|l| l.get_id()).collect();
        
        for i in 0..n_taxa {
            lca_depths[i][i] = tree.depth(leaf_ids[i]); // Distance to self
            for j in i+1..n_taxa {
                let lca = tree.get_lca_id(&[leaf_ids[i], leaf_ids[j]]);
                let depth = tree.depth(lca);
                lca_depths[i][j] = depth;
                lca_depths[j][i] = depth; // Symmetric
            }
        }
        
        UltraFastTreeData {
            taxa_to_idx,
            idx_to_taxa,
            n_taxa,
            lca_depths,
        }
    }
    
    /// Ultra-fast outgroup detection with direct array access
    #[inline]
    fn get_outgroup_ultra_fast(&self, i_idx: usize, j_idx: usize, k_idx: usize) -> Option<usize> {
        let depth_ij = self.lca_depths[i_idx][j_idx];
        let depth_ik = self.lca_depths[i_idx][k_idx];
        let depth_jk = self.lca_depths[j_idx][k_idx];
        
        if depth_ij > depth_ik && depth_ij > depth_jk {
            Some(k_idx) // k is outgroup
        } else if depth_ik > depth_ij && depth_ik > depth_jk {
            Some(j_idx) // j is outgroup  
        } else if depth_jk > depth_ij && depth_jk > depth_ik {
            Some(i_idx) // i is outgroup
        } else {
            None // Unresolvable
        }
    }
    
    fn taxa_to_idx(&self, taxa: &str) -> Option<usize> {
        self.taxa_to_idx.get(taxa).copied()
    }
    
    fn idx_to_taxa(&self, idx: usize) -> &str {
        &self.idx_to_taxa[idx]
    }
}

/// Calculate total number of possible triplets
fn total_triplets(n: usize) -> usize {
    if n < 3 { 0 } else { n * (n - 1) * (n - 2) / 6 }
}

/// Determine if exhaustive calculation is better than sampling
fn should_use_exhaustive(n_taxa: usize, requested_trials: usize) -> bool {
    let total = total_triplets(n_taxa);
    // Use exhaustive if:
    // 1. Total triplets <= requested trials (exhaustive is more accurate)
    // 2. Total triplets is reasonably small (exhaustive is faster)
    // 3. For very small trees, always use exhaustive
    if n_taxa <= 20 {
        true  // Always use exhaustive for small trees
    } else {
        total <= requested_trials || total <= 100000
    }
}

/// Exhaustive triplet calculation (for small-medium trees)
fn calculate_exhaustive(
    tree1_data: &UltraFastTreeData,
    tree2_data: &UltraFastTreeData,
) -> (usize, usize, usize, usize) {
    let mut all_correct = 0;
    let mut resolvable_correct = 0;
    let mut unresolvable_correct = 0;
    let mut unresolvable_count = 0;
    
    let n = tree1_data.n_taxa;
    
    // Generate all possible triplets
    for i in 0..n {
        for j in i+1..n {
            for k in j+1..n {
                // Get outgroups from both trees using indices
                let outgroup1 = tree1_data.get_outgroup_ultra_fast(i, j, k);
                let outgroup2 = tree2_data.get_outgroup_ultra_fast(i, j, k);
                
                let is_resolvable = outgroup1.is_some();
                if !is_resolvable {
                    unresolvable_count += 1;
                }
                
                let correct = outgroup1 == outgroup2;
                if correct {
                    all_correct += 1;
                    if is_resolvable {
                        resolvable_correct += 1;
                    } else {
                        unresolvable_correct += 1;
                    }
                }
            }
        }
    }
    
    (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count)
}

/// Sampling-based calculation (for large trees)
fn calculate_sampling(
    tree1_data: &UltraFastTreeData,
    tree2_data: &UltraFastTreeData,
    number_of_trials: usize,
    seed: Option<u64>,
) -> (usize, usize, usize, usize) {
    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };
    
    let mut all_correct = 0;
    let mut resolvable_correct = 0;
    let mut unresolvable_correct = 0;
    let mut unresolvable_count = 0;
    
    let n = tree1_data.n_taxa;
    let indices: Vec<usize> = (0..n).collect();
    
    for _ in 0..number_of_trials {
        // Sample 3 distinct indices
        let selected: Vec<_> = indices.choose_multiple(&mut rng, 3).collect();
        let (i, j, k) = (*selected[0], *selected[1], *selected[2]);
        
        let outgroup1 = tree1_data.get_outgroup_ultra_fast(i, j, k);
        let outgroup2 = tree2_data.get_outgroup_ultra_fast(i, j, k);
        
        let is_resolvable = outgroup1.is_some();
        if !is_resolvable {
            unresolvable_count += 1;
        }
        
        let correct = outgroup1 == outgroup2;
        if correct {
            all_correct += 1;
            if is_resolvable {
                resolvable_correct += 1;
            } else {
                unresolvable_correct += 1;
            }
        }
    }
    
    (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count)
}

/// Ultra-fast triplets correct with automatic exhaustive/sampling selection
pub fn ultra_fast_triplets_correct(
    tree1: &mut PhyloTree,
    tree2: &mut PhyloTree,
    number_of_trials: usize,
    _min_triplets_at_depth: usize,
    seed: Option<u64>,
) -> UltraFastTripletResult {
    // Precompute data for both trees
    let tree1_data = UltraFastTreeData::new(tree1);
    let tree2_data = UltraFastTreeData::new(tree2);
    
    // Verify both trees have same taxa
    if tree1_data.n_taxa != tree2_data.n_taxa {
        panic!("Trees must have the same number of taxa");
    }
    
    let n_taxa = tree1_data.n_taxa;
    let use_exhaustive = should_use_exhaustive(n_taxa, number_of_trials);
    
    let (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count, total_checked, method) = 
        if use_exhaustive {
            let (a, r, u, uc) = calculate_exhaustive(&tree1_data, &tree2_data);
            (a, r, u, uc, total_triplets(n_taxa), "exhaustive".to_string())
        } else {
            let (a, r, u, uc) = calculate_sampling(&tree1_data, &tree2_data, number_of_trials, seed);
            (a, r, u, uc, number_of_trials, "sampling".to_string())
        };
    
    // Calculate results
    let mut all_triplets_correct = HashMap::new();
    let mut resolvable_triplets_correct = HashMap::new();
    let mut unresolved_triplets_correct = HashMap::new();
    let mut proportion_unresolvable = HashMap::new();
    
    all_triplets_correct.insert(0, all_correct as f64 / total_checked as f64);
    proportion_unresolvable.insert(0, unresolvable_count as f64 / total_checked as f64);
    
    if unresolvable_count == 0 {
        unresolved_triplets_correct.insert(0, 1.0);
    } else {
        unresolved_triplets_correct.insert(0, unresolvable_correct as f64 / unresolvable_count as f64);
    }
    
    let resolvable_trials = total_checked - unresolvable_count;
    if resolvable_trials > 0 {
        resolvable_triplets_correct.insert(0, resolvable_correct as f64 / resolvable_trials as f64);
    } else {
        resolvable_triplets_correct.insert(0, 1.0);
    }
    
    UltraFastTripletResult {
        all_triplets_correct,
        resolvable_triplets_correct,
        unresolved_triplets_correct,
        proportion_unresolvable,
        method_used: method,
        total_triplets_checked: total_checked,
    }
}