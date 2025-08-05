use phylo::prelude::*;
use rand::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

/// Parallel triplet result with performance metadata
#[derive(Debug, Clone)]
pub struct ParallelTripletResult {
    pub all_triplets_correct: HashMap<usize, f64>,
    pub resolvable_triplets_correct: HashMap<usize, f64>,
    pub unresolved_triplets_correct: HashMap<usize, f64>,
    pub proportion_unresolvable: HashMap<usize, f64>,
    pub method_used: String,
    pub total_triplets_checked: usize,
    pub parallel_chunks_used: usize,
}

/// Thread-safe tree data for parallel processing
struct ParallelTreeData {
    taxa_to_idx: HashMap<String, usize>,
    idx_to_taxa: Vec<String>,
    n_taxa: usize,
    // Flat 2D array for O(1) LCA depth access - read-only after construction
    lca_depths: Vec<Vec<usize>>,
}

impl ParallelTreeData {
    fn new(tree: &mut PhyloTree) -> Self {
        tree.precompute_constant_time_lca();
        
        let mut taxa_to_idx = HashMap::new();
        let mut idx_to_taxa = Vec::new();
        
        for (idx, leaf) in tree.get_leaves().enumerate() {
            if let Some(taxa) = leaf.get_taxa() {
                taxa_to_idx.insert(taxa.clone(), idx);
                idx_to_taxa.push(taxa.clone());
            }
        }
        
        let n_taxa = idx_to_taxa.len();
        let mut lca_depths = vec![vec![0; n_taxa]; n_taxa];
        
        let leaf_ids: Vec<_> = tree.get_leaves().map(|l| l.get_id()).collect();
        
        for i in 0..n_taxa {
            lca_depths[i][i] = tree.depth(leaf_ids[i]);
            for j in i+1..n_taxa {
                let lca = tree.get_lca_id(&[leaf_ids[i], leaf_ids[j]]);
                let depth = tree.depth(lca);
                lca_depths[i][j] = depth;
                lca_depths[j][i] = depth;
            }
        }
        
        ParallelTreeData {
            taxa_to_idx,
            idx_to_taxa,
            n_taxa,
            lca_depths,
        }
    }
    
    #[inline]
    fn get_outgroup_fast(&self, i_idx: usize, j_idx: usize, k_idx: usize) -> Option<usize> {
        let depth_ij = self.lca_depths[i_idx][j_idx];
        let depth_ik = self.lca_depths[i_idx][k_idx];
        let depth_jk = self.lca_depths[j_idx][k_idx];
        
        if depth_ij > depth_ik && depth_ij > depth_jk {
            Some(k_idx)
        } else if depth_ik > depth_ij && depth_ik > depth_jk {
            Some(j_idx)
        } else if depth_jk > depth_ij && depth_jk > depth_ik {
            Some(i_idx)
        } else {
            None
        }
    }
}

/// Generate all triplets as a flat vector for balanced work distribution
fn generate_all_triplets(n: usize) -> Vec<(usize, usize, usize)> {
    let mut triplets = Vec::with_capacity(total_triplets(n));
    for i in 0..n {
        for j in i+1..n {
            for k in j+1..n {
                triplets.push((i, j, k));
            }
        }
    }
    triplets
}

/// Dynamically determine optimal thread count based on workload and available cores
fn optimal_thread_count(total_work: usize, max_threads: Option<usize>) -> usize {
    let mut available_threads = rayon::current_num_threads();
    
    // Apply user-specified thread limit if provided
    if let Some(max) = max_threads {
        available_threads = std::cmp::min(available_threads, max);
    }
    
    // Base workload per thread - empirically determined optimal range
    let work_per_thread = 2000..8000;
    
    if total_work < work_per_thread.start {
        return 1;  // Sequential for very small workloads
    }
    
    // Calculate optimal threads based on work distribution
    let ideal_threads = total_work / work_per_thread.end;  // Conservative estimate
    let max_beneficial_threads = (total_work / work_per_thread.start).max(1);  // Aggressive estimate
    
    // Use geometric mean between conservative and aggressive, bounded by available cores
    let target_threads = ((ideal_threads * max_beneficial_threads) as f64).sqrt() as usize;
    
    // Ensure we don't exceed available threads and have at least 1
    std::cmp::min(std::cmp::max(target_threads, 1), available_threads)
}

/// Optimized parallel exhaustive calculation with balanced work distribution
fn calculate_exhaustive_parallel(
    tree1_data: &ParallelTreeData,
    tree2_data: &ParallelTreeData,
    max_threads: Option<usize>,
) -> (usize, usize, usize, usize, usize) {
    let n = tree1_data.n_taxa;
    let total_triplet_count = total_triplets(n);
    let target_threads = optimal_thread_count(total_triplet_count, max_threads);
    
    // For small workloads, use sequential processing to avoid parallelization overhead
    if total_triplet_count < 10000 || target_threads == 1 {
        let mut all_correct = 0;
        let mut resolvable_correct = 0;
        let mut unresolvable_correct = 0;
        let mut unresolvable_count = 0;
        
        for i in 0..n {
            for j in i+1..n {
                for k in j+1..n {
                    let outgroup1 = tree1_data.get_outgroup_fast(i, j, k);
                    let outgroup2 = tree2_data.get_outgroup_fast(i, j, k);
                    
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
        
        return (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count, 1);
    }
    
    // Generate all triplets and distribute with adaptive chunk size for optimal load balancing
    let all_triplets = generate_all_triplets(n);
    let base_chunk_size = all_triplets.len() / target_threads;
    
    // Adaptive chunk sizing: more chunks for better load balancing with many threads
    // Target 2-4x more chunks than threads for good work stealing
    let chunk_multiplier = if target_threads >= 50 {
        4  // Many threads - create many small chunks
    } else if target_threads >= 16 {
        3  // Moderate threads - moderate chunking
    } else if target_threads >= 4 {
        2  // Few threads - some extra chunks
    } else {
        1  // Very few threads - minimal chunking
    };
    
    let chunk_size = std::cmp::max(1, base_chunk_size / chunk_multiplier);
    
    let results: Vec<(usize, usize, usize, usize)> = all_triplets
        .par_chunks(chunk_size)
        .map(|triplet_chunk| {
            let mut local_all_correct = 0;
            let mut local_resolvable_correct = 0;
            let mut local_unresolvable_correct = 0;
            let mut local_unresolvable_count = 0;
            
            for &(i, j, k) in triplet_chunk {
                let outgroup1 = tree1_data.get_outgroup_fast(i, j, k);
                let outgroup2 = tree2_data.get_outgroup_fast(i, j, k);
                
                let is_resolvable = outgroup1.is_some();
                if !is_resolvable {
                    local_unresolvable_count += 1;
                }
                
                let correct = outgroup1 == outgroup2;
                if correct {
                    local_all_correct += 1;
                    if is_resolvable {
                        local_resolvable_correct += 1;
                    } else {
                        local_unresolvable_correct += 1;
                    }
                }
            }
            
            (local_all_correct, local_resolvable_correct, local_unresolvable_correct, local_unresolvable_count)
        })
        .collect();
    
    // Combine results from all threads
    let actual_chunks = std::cmp::min(results.len(), target_threads);
    let (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count) = 
        results.into_iter().reduce(
            |acc, item| (acc.0 + item.0, acc.1 + item.1, acc.2 + item.2, acc.3 + item.3)
        ).unwrap_or((0, 0, 0, 0));
    (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count, actual_chunks)
}

/// Optimized parallel sampling calculation
fn calculate_sampling_parallel(
    tree1_data: &ParallelTreeData,
    tree2_data: &ParallelTreeData,
    number_of_trials: usize,
    seed: Option<u64>,
    max_threads: Option<usize>,
) -> (usize, usize, usize, usize, usize) {
    let n = tree1_data.n_taxa;
    let target_threads = optimal_thread_count(number_of_trials, max_threads);
    
    // For small trial counts, use sequential processing to avoid overhead
    if number_of_trials < 5000 || target_threads == 1 {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };
        
        let mut all_correct = 0;
        let mut resolvable_correct = 0;
        let mut unresolvable_correct = 0;
        let mut unresolvable_count = 0;
        
        let indices: Vec<usize> = (0..n).collect();
        
        for _ in 0..number_of_trials {
            let selected: Vec<_> = indices.choose_multiple(&mut rng, 3).collect();
            let (i, j, k) = (*selected[0], *selected[1], *selected[2]);
            
            let outgroup1 = tree1_data.get_outgroup_fast(i, j, k);
            let outgroup2 = tree2_data.get_outgroup_fast(i, j, k);
            
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
        
        return (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count, 1);
    }
    
    let trials_per_thread = number_of_trials / target_threads;
    let extra_trials = number_of_trials % target_threads;
    
    // Generate seeds for each thread to ensure reproducibility
    let thread_seeds: Vec<u64> = {
        let mut base_rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };
        (0..target_threads).map(|_| base_rng.gen()).collect()
    };
    
    let results: Vec<(usize, usize, usize, usize)> = (0..target_threads)
        .into_par_iter()
        .map(|thread_id| {
            let mut thread_rng = StdRng::seed_from_u64(thread_seeds[thread_id]);
            let thread_trials = trials_per_thread + if thread_id < extra_trials { 1 } else { 0 };
            
            let mut local_all_correct = 0;
            let mut local_resolvable_correct = 0;
            let mut local_unresolvable_correct = 0;
            let mut local_unresolvable_count = 0;
            
            let indices: Vec<usize> = (0..n).collect();
            
            for _ in 0..thread_trials {
                let selected: Vec<_> = indices.choose_multiple(&mut thread_rng, 3).collect();
                let (i, j, k) = (*selected[0], *selected[1], *selected[2]);
                
                let outgroup1 = tree1_data.get_outgroup_fast(i, j, k);
                let outgroup2 = tree2_data.get_outgroup_fast(i, j, k);
                
                let is_resolvable = outgroup1.is_some();
                if !is_resolvable {
                    local_unresolvable_count += 1;
                }
                
                let correct = outgroup1 == outgroup2;
                if correct {
                    local_all_correct += 1;
                    if is_resolvable {
                        local_resolvable_correct += 1;
                    } else {
                        local_unresolvable_correct += 1;
                    }
                }
            }
            
            (local_all_correct, local_resolvable_correct, local_unresolvable_correct, local_unresolvable_count)
        })
        .collect();
    
    // Combine results from all threads
    let (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count) = 
        results.into_iter().reduce(
            |acc, item| (acc.0 + item.0, acc.1 + item.1, acc.2 + item.2, acc.3 + item.3)
        ).unwrap_or((0, 0, 0, 0));
    
    (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count, target_threads)
}

/// Calculate total number of possible triplets
fn total_triplets(n: usize) -> usize {
    if n < 3 { 0 } else { n * (n - 1) * (n - 2) / 6 }
}

/// Dynamically determine if exhaustive calculation is better than sampling
fn should_use_exhaustive(n_taxa: usize, requested_trials: usize, max_threads: Option<usize>) -> bool {
    let total = total_triplets(n_taxa);
    let mut available_threads = rayon::current_num_threads();
    
    // Apply user-specified thread limit for consistency with optimal_thread_count
    if let Some(max) = max_threads {
        available_threads = std::cmp::min(available_threads, max);
    }
    
    if n_taxa <= 25 {
        return true;  // Always use exhaustive for small trees
    }
    
    // Dynamic threshold based on available parallelization
    // Scale threshold roughly linearly with available cores, with diminishing returns
    let base_threshold = 100_000;  // Baseline for single core
    let thread_multiplier = (available_threads as f64).sqrt();  // Diminishing returns
    let parallel_threshold = (base_threshold as f64 * thread_multiplier) as usize;
    
    // Use exhaustive if total work is less than requested trials or our parallel threshold
    total <= requested_trials || total <= parallel_threshold
}

/// Parallel triplets correct with automatic exhaustive/sampling selection
pub fn parallel_triplets_correct(
    tree1: &mut PhyloTree,
    tree2: &mut PhyloTree,
    number_of_trials: usize,
    _min_triplets_at_depth: usize,
    seed: Option<u64>,
    max_threads: Option<usize>,
) -> ParallelTripletResult {
    // Precompute data for both trees
    let tree1_data = ParallelTreeData::new(tree1);
    let tree2_data = ParallelTreeData::new(tree2);
    
    if tree1_data.n_taxa != tree2_data.n_taxa {
        panic!("Trees must have the same number of taxa");
    }
    
    let n_taxa = tree1_data.n_taxa;
    let use_exhaustive = should_use_exhaustive(n_taxa, number_of_trials, max_threads);
    
    let (all_correct, resolvable_correct, unresolvable_correct, unresolvable_count, total_checked, method, chunks) = 
        if use_exhaustive {
            let (a, r, u, uc, ch) = calculate_exhaustive_parallel(&tree1_data, &tree2_data, max_threads);
            (a, r, u, uc, total_triplets(n_taxa), "parallel_exhaustive".to_string(), ch)
        } else {
            let (a, r, u, uc, ch) = calculate_sampling_parallel(&tree1_data, &tree2_data, number_of_trials, seed, max_threads);
            (a, r, u, uc, number_of_trials, "parallel_sampling".to_string(), ch)
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
    
    ParallelTripletResult {
        all_triplets_correct,
        resolvable_triplets_correct,
        unresolved_triplets_correct,
        proportion_unresolvable,
        method_used: method,
        total_triplets_checked: total_checked,
        parallel_chunks_used: chunks,
    }
}