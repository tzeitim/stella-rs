use phylo::prelude::*;
use std::time::Instant;

/// Result structure for RF distance calculation
#[derive(Debug, Clone)]
pub struct RFDistanceResult {
    pub rf_distance: usize,
    pub computation_time_ms: f64,
    pub method_used: String,
}

/// Calculate Robinson-Foulds distance between two phylogenetic trees
pub fn calculate_rf_distance(
    tree1: &PhyloTree,
    tree2: &PhyloTree,
) -> RFDistanceResult {
    let start_time = Instant::now();
    
    // Use phylo-rs RobinsonFoulds trait to calculate RF distance
    let rf_distance = tree1.rf(tree2);
    
    let computation_time_ms = start_time.elapsed().as_secs_f64() * 1000.0;
    
    RFDistanceResult {
        rf_distance,
        computation_time_ms,
        method_used: "phylo-rs RobinsonFoulds".to_string(),
    }
}

/// Calculate normalized Robinson-Foulds distance (RF distance / max possible RF distance)
pub fn calculate_normalized_rf_distance(
    tree1: &PhyloTree,
    tree2: &PhyloTree,
) -> (RFDistanceResult, f64) {
    let result = calculate_rf_distance(tree1, tree2);
    
    // Calculate maximum possible RF distance
    // For two binary trees with n leaves, max RF distance is 2(n-3)
    // We use the minimum number of leaves between the two trees
    let n_leaves_tree1 = tree1.get_taxa_space().count();
    let n_leaves_tree2 = tree2.get_taxa_space().count();
    let min_leaves = n_leaves_tree1.min(n_leaves_tree2);
    
    let max_rf_distance = if min_leaves > 2 {
        2 * (min_leaves - 3)
    } else {
        1 // Avoid division by zero for very small trees
    };
    
    let normalized_rf = if max_rf_distance > 0 {
        result.rf_distance as f64 / max_rf_distance as f64
    } else {
        0.0
    };
    
    (result, normalized_rf)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_rf_distance_identical_trees() {
        let newick = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);";
        
        let tree1 = PhyloTree::from_newick(newick.as_bytes()).unwrap();
        let tree2 = PhyloTree::from_newick(newick.as_bytes()).unwrap();
        
        let result = calculate_rf_distance(&tree1, &tree2);
        
        assert_eq!(result.rf_distance, 0);
        assert!(result.computation_time_ms >= 0.0);
        assert_eq!(result.method_used, "phylo-rs RobinsonFoulds");
    }
    
    #[test]
    fn test_rf_distance_different_trees() {
        let newick1 = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);";
        let newick2 = "(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);";
        
        let tree1 = PhyloTree::from_newick(newick1.as_bytes()).unwrap();
        let tree2 = PhyloTree::from_newick(newick2.as_bytes()).unwrap();
        
        let result = calculate_rf_distance(&tree1, &tree2);
        
        assert!(result.rf_distance > 0);
        assert!(result.computation_time_ms >= 0.0);
    }
    
    #[test]
    fn test_normalized_rf_distance() {
        let newick1 = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);";
        let newick2 = "(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);";
        
        let tree1 = PhyloTree::from_newick(newick1.as_bytes()).unwrap();
        let tree2 = PhyloTree::from_newick(newick2.as_bytes()).unwrap();
        
        let (result, normalized) = calculate_normalized_rf_distance(&tree1, &tree2);
        
        assert!(result.rf_distance > 0);
        assert!(normalized >= 0.0 && normalized <= 1.0);
    }
}