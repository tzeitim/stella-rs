use phylo::prelude::*;
use rand::prelude::*;
use std::collections::HashMap;

/// Fast triplet result using precomputed data
#[derive(Debug, Clone)]
pub struct FastTripletResult {
    pub all_triplets_correct: HashMap<usize, f64>,
    pub resolvable_triplets_correct: HashMap<usize, f64>,
    pub unresolved_triplets_correct: HashMap<usize, f64>,
    pub proportion_unresolvable: HashMap<usize, f64>,
}

/// Precomputed tree data for fast triplet queries
struct FastTreeData {
    taxa_to_id: HashMap<String, usize>,
    id_to_taxa: HashMap<usize, String>,
    leaf_ids: Vec<usize>,
    // Cache all pairwise LCA depths
    lca_depths: HashMap<(usize, usize), usize>,
}

impl FastTreeData {
    fn new(tree: &mut PhyloTree) -> Self {
        // Precompute LCA for O(1) queries
        tree.precompute_constant_time_lca();
        
        let mut taxa_to_id = HashMap::new();
        let mut id_to_taxa = HashMap::new();
        let mut leaf_ids = Vec::new();
        
        // Build taxa mappings
        for leaf in tree.get_leaves() {
            if let Some(taxa) = leaf.get_taxa() {
                let id = leaf.get_id();
                taxa_to_id.insert(taxa.clone(), id);
                id_to_taxa.insert(id, taxa.clone());
                leaf_ids.push(id);
            }
        }
        
        // Precompute all pairwise LCA depths
        let mut lca_depths = HashMap::new();
        for i in 0..leaf_ids.len() {
            for j in i+1..leaf_ids.len() {
                let id1 = leaf_ids[i];
                let id2 = leaf_ids[j];
                let lca = tree.get_lca_id(&[id1, id2]);
                let depth = tree.depth(lca);
                lca_depths.insert((id1, id2), depth);
                lca_depths.insert((id2, id1), depth); // Symmetric
            }
        }
        
        FastTreeData {
            taxa_to_id,
            id_to_taxa,
            leaf_ids,
            lca_depths,
        }
    }
    
    fn get_outgroup_fast(&self, triplet: (&str, &str, &str)) -> Option<String> {
        let (i, j, k) = triplet;
        
        // Fast taxa lookup
        let i_id = *self.taxa_to_id.get(i)?;
        let j_id = *self.taxa_to_id.get(j)?;
        let k_id = *self.taxa_to_id.get(k)?;
        
        // Fast LCA depth lookup
        let depth_ij = *self.lca_depths.get(&(i_id, j_id))?;
        let depth_ik = *self.lca_depths.get(&(i_id, k_id))?;
        let depth_jk = *self.lca_depths.get(&(j_id, k_id))?;
        
        // Same logic as before but with cached data
        if depth_ij > depth_ik && depth_ij > depth_jk {
            Some(k.to_string())
        } else if depth_ik > depth_ij && depth_ik > depth_jk {
            Some(j.to_string())
        } else if depth_jk > depth_ij && depth_jk > depth_ik {
            Some(i.to_string())
        } else {
            None
        }
    }
    
    fn sample_triplet_fast(&self, rng: &mut impl Rng) -> Option<(String, String, String)> {
        if self.leaf_ids.len() < 3 {
            return None;
        }
        
        let selected_ids: Vec<_> = self.leaf_ids.choose_multiple(rng, 3).collect();
        let taxa1 = self.id_to_taxa.get(selected_ids[0])?.clone();
        let taxa2 = self.id_to_taxa.get(selected_ids[1])?.clone();
        let taxa3 = self.id_to_taxa.get(selected_ids[2])?.clone();
        
        Some((taxa1, taxa2, taxa3))
    }
}

/// Fast triplets correct implementation with precomputed data
pub fn fast_triplets_correct(
    tree1: &mut PhyloTree,
    tree2: &mut PhyloTree,
    number_of_trials: usize,
    _min_triplets_at_depth: usize,
    seed: Option<u64>,
) -> FastTripletResult {
    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };
    
    // Precompute data for both trees
    let tree1_data = FastTreeData::new(tree1);
    let tree2_data = FastTreeData::new(tree2);
    
    let mut all_correct = 0;
    let mut resolvable_correct = 0;
    let mut unresolvable_correct = 0;
    let mut unresolvable_count = 0;
    
    // Fast sampling and comparison
    for _ in 0..number_of_trials {
        if let Some((i, j, k)) = tree1_data.sample_triplet_fast(&mut rng) {
            // Fast outgroup detection
            let truth_outgroup = tree1_data.get_outgroup_fast((&i, &j, &k));
            let reconstructed_outgroup = tree2_data.get_outgroup_fast((&i, &j, &k));
            
            let is_resolvable = truth_outgroup.is_some();
            
            if !is_resolvable {
                unresolvable_count += 1;
            }
            
            let correct = match (&truth_outgroup, &reconstructed_outgroup) {
                (Some(t), Some(r)) => t == r,
                (None, None) => true,
                _ => false,
            };
            
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
    
    // Calculate results
    let mut all_triplets_correct = HashMap::new();
    let mut resolvable_triplets_correct = HashMap::new();
    let mut unresolved_triplets_correct = HashMap::new();
    let mut proportion_unresolvable = HashMap::new();
    
    all_triplets_correct.insert(0, all_correct as f64 / number_of_trials as f64);
    proportion_unresolvable.insert(0, unresolvable_count as f64 / number_of_trials as f64);
    
    if unresolvable_count == 0 {
        unresolved_triplets_correct.insert(0, 1.0);
    } else {
        unresolved_triplets_correct.insert(0, unresolvable_correct as f64 / unresolvable_count as f64);
    }
    
    let resolvable_trials = number_of_trials - unresolvable_count;
    if resolvable_trials > 0 {
        resolvable_triplets_correct.insert(0, resolvable_correct as f64 / resolvable_trials as f64);
    } else {
        resolvable_triplets_correct.insert(0, 1.0);
    }
    
    FastTripletResult {
        all_triplets_correct,
        resolvable_triplets_correct,
        unresolved_triplets_correct,
        proportion_unresolvable,
    }
}