use phylo::prelude::*;
use rand::prelude::*;
use std::collections::HashMap;

/// Represents the result of triplet analysis at different depths
#[derive(Debug, Clone)]
pub struct TripletResult {
    pub all_triplets_correct: HashMap<usize, f64>,
    pub resolvable_triplets_correct: HashMap<usize, f64>,
    pub unresolved_triplets_correct: HashMap<usize, f64>,
    pub proportion_unresolvable: HashMap<usize, f64>,
}


/// Get the outgroup of a triplet based on LCA depths
pub fn get_outgroup(tree: &PhyloTree, triplet: (&str, &str, &str)) -> Option<String> {
    let (i, j, k) = triplet;
    
    // Get node IDs for the leaves by taxa
    let i_id = tree.get_taxa_node_id(&i.to_string())
        .or_else(|| {
            // Fallback: find leaf with matching taxa
            tree.get_leaves().find(|leaf| {
                leaf.get_taxa().map_or(false, |taxa| taxa == i)
            }).map(|leaf| leaf.get_id())
        })?;
    
    let j_id = tree.get_taxa_node_id(&j.to_string())
        .or_else(|| {
            tree.get_leaves().find(|leaf| {
                leaf.get_taxa().map_or(false, |taxa| taxa == j)
            }).map(|leaf| leaf.get_id())
        })?;
    
    let k_id = tree.get_taxa_node_id(&k.to_string())
        .or_else(|| {
            tree.get_leaves().find(|leaf| {
                leaf.get_taxa().map_or(false, |taxa| taxa == k)
            }).map(|leaf| leaf.get_id())
        })?;
    
    // Get LCA for each pair
    let lca_ij = tree.get_lca_id(&[i_id, j_id]);
    let lca_ik = tree.get_lca_id(&[i_id, k_id]);
    let lca_jk = tree.get_lca_id(&[j_id, k_id]);
    
    // Get depths of LCAs
    let depth_ij = tree.depth(lca_ij);
    let depth_ik = tree.depth(lca_ik);
    let depth_jk = tree.depth(lca_jk);
    
    
    // The outgroup is the one whose LCA with the ingroup pair is shallowest
    // (the ingroup pair has the deepest/highest LCA depth)
    if depth_ij > depth_ik && depth_ij > depth_jk {
        // i and j are the ingroup, k is outgroup
        Some(k.to_string())
    } else if depth_ik > depth_ij && depth_ik > depth_jk {
        // i and k are the ingroup, j is outgroup
        Some(j.to_string())
    } else if depth_jk > depth_ij && depth_jk > depth_ik {
        // j and k are the ingroup, i is outgroup
        Some(i.to_string())
    } else {
        None // Unresolvable triplet (all LCAs at same depth)
    }
}

/// Sample a random triplet from the tree leaves
pub fn sample_random_triplet(tree: &PhyloTree, rng: &mut impl Rng) -> Option<(String, String, String)> {
    let leaves: Vec<_> = tree.get_leaves().collect();
    
    if leaves.len() < 3 {
        return None;
    }
    
    let selected_leaves: Vec<_> = leaves.choose_multiple(rng, 3).collect();
    
    let taxa1 = selected_leaves[0].get_taxa()?.clone();
    let taxa2 = selected_leaves[1].get_taxa()?.clone();
    let taxa3 = selected_leaves[2].get_taxa()?.clone();
    
    Some((taxa1, taxa2, taxa3))
}

/// Calculate triplets correct between two trees using simple random sampling
pub fn triplets_correct_impl(
    tree1: &PhyloTree,
    tree2: &PhyloTree,
    number_of_trials: usize,
    _min_triplets_at_depth: usize, // Simplified: ignore depth-based sampling for now
    seed: Option<u64>,
) -> TripletResult {
    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };
    
    let mut all_correct = 0;
    let mut resolvable_correct = 0;
    let mut unresolvable_correct = 0;
    let mut unresolvable_count = 0;
    
    // Sample triplets randomly from tree1
    for _ in 0..number_of_trials {
        if let Some((i, j, k)) = sample_random_triplet(tree1, &mut rng) {
            // Get true outgroup from tree1
            let truth_outgroup = get_outgroup(tree1, (&i, &j, &k));
            
            // Get reconstructed outgroup from tree2
            let reconstructed_outgroup = get_outgroup(tree2, (&i, &j, &k));
            
            let is_resolvable = truth_outgroup.is_some();
            
            if !is_resolvable {
                unresolvable_count += 1;
            }
            
            // Check if outgroups match
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
    
    // Calculate proportions (using depth 0 as a single result)
    let mut all_triplets_correct = HashMap::new();
    let mut resolvable_triplets_correct = HashMap::new();
    let mut unresolved_triplets_correct = HashMap::new();
    let mut proportion_unresolvable = HashMap::new();
    
    all_triplets_correct.insert(0, all_correct as f64 / number_of_trials as f64);
    
    let prop_unresolvable = unresolvable_count as f64 / number_of_trials as f64;
    proportion_unresolvable.insert(0, prop_unresolvable);
    
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
    
    TripletResult {
        all_triplets_correct,
        resolvable_triplets_correct,
        unresolved_triplets_correct,
        proportion_unresolvable,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    
    #[test]
    fn test_triplet_outgroup() {
        // Create a simple tree: ((A,B),C);
        let newick = "((A,B),C);";
        let mut tree = PhyloTree::from_newick(newick.as_bytes()).unwrap();
        
        // Precompute LCA for efficient queries
        tree.precompute_constant_time_lca();
        
        // A and B should be closer to each other than to C
        assert_eq!(get_outgroup(&tree, ("A", "B", "C")), Some("C".to_string()));
        assert_eq!(get_outgroup(&tree, ("A", "C", "B")), Some("C".to_string()));
        assert_eq!(get_outgroup(&tree, ("B", "C", "A")), Some("C".to_string()));
    }
}