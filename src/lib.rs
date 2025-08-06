use pyo3::prelude::*;
use pyo3::types::PyDict;
use phylo::prelude::*;

mod triplets;
mod fast_triplets;
mod ultra_fast_triplets;
mod parallel_triplets;
mod optimized_triplets;
use triplets::triplets_correct_impl;
use fast_triplets::fast_triplets_correct;
use ultra_fast_triplets::ultra_fast_triplets_correct;
use parallel_triplets::parallel_triplets_correct;
use optimized_triplets::optimized_triplets_correct;

/// Python wrapper for triplets_correct function
#[pyfunction]
#[pyo3(signature = (tree1_newick, tree2_newick, number_of_trials=1000, min_triplets_at_depth=1, seed=None))]
fn triplets_correct(
    py: Python,
    tree1_newick: &str,
    tree2_newick: &str,
    number_of_trials: usize,
    min_triplets_at_depth: usize,
    seed: Option<u64>,
) -> PyResult<PyObject> {
    // Parse trees from Newick strings
    let tree1 = PhyloTree::from_newick(tree1_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree1: {:?}", e)))?;
    
    let tree2 = PhyloTree::from_newick(tree2_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree2: {:?}", e)))?;
    
    // Calculate triplets correct
    let result = triplets_correct_impl(&tree1, &tree2, number_of_trials, min_triplets_at_depth, seed);
    
    // Convert result to Python dictionaries
    let py_result = PyDict::new_bound(py);
    
    // Convert HashMaps to Python dicts
    let all_dict = PyDict::new_bound(py);
    for (k, v) in result.all_triplets_correct {
        all_dict.set_item(k, v)?;
    }
    
    let resolvable_dict = PyDict::new_bound(py);
    for (k, v) in result.resolvable_triplets_correct {
        resolvable_dict.set_item(k, v)?;
    }
    
    let unresolved_dict = PyDict::new_bound(py);
    for (k, v) in result.unresolved_triplets_correct {
        unresolved_dict.set_item(k, v)?;
    }
    
    let proportion_dict = PyDict::new_bound(py);
    for (k, v) in result.proportion_unresolvable {
        proportion_dict.set_item(k, v)?;
    }
    
    py_result.set_item("all_triplets_correct", all_dict)?;
    py_result.set_item("resolvable_triplets_correct", resolvable_dict)?;
    py_result.set_item("unresolved_triplets_correct", unresolved_dict)?;
    py_result.set_item("proportion_unresolvable", proportion_dict)?;
    
    Ok(py_result.into())
}

/// Fast Python wrapper for triplets_correct function
#[pyfunction]
#[pyo3(signature = (tree1_newick, tree2_newick, number_of_trials=1000, min_triplets_at_depth=1, seed=None))]
fn triplets_correct_fast(
    py: Python,
    tree1_newick: &str,
    tree2_newick: &str,
    number_of_trials: usize,
    min_triplets_at_depth: usize,
    seed: Option<u64>,
) -> PyResult<PyObject> {
    // Parse trees from Newick strings
    let mut tree1 = PhyloTree::from_newick(tree1_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree1: {:?}", e)))?;
    
    let mut tree2 = PhyloTree::from_newick(tree2_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree2: {:?}", e)))?;
    
    // Calculate triplets correct using fast implementation
    let result = fast_triplets_correct(&mut tree1, &mut tree2, number_of_trials, min_triplets_at_depth, seed);
    
    // Convert result to Python dictionaries (same as before)
    let py_result = PyDict::new_bound(py);
    
    let all_dict = PyDict::new_bound(py);
    for (k, v) in result.all_triplets_correct {
        all_dict.set_item(k, v)?;
    }
    
    let resolvable_dict = PyDict::new_bound(py);
    for (k, v) in result.resolvable_triplets_correct {
        resolvable_dict.set_item(k, v)?;
    }
    
    let unresolved_dict = PyDict::new_bound(py);
    for (k, v) in result.unresolved_triplets_correct {
        unresolved_dict.set_item(k, v)?;
    }
    
    let proportion_dict = PyDict::new_bound(py);
    for (k, v) in result.proportion_unresolvable {
        proportion_dict.set_item(k, v)?;
    }
    
    py_result.set_item("all_triplets_correct", all_dict)?;
    py_result.set_item("resolvable_triplets_correct", resolvable_dict)?;
    py_result.set_item("unresolved_triplets_correct", unresolved_dict)?;
    py_result.set_item("proportion_unresolvable", proportion_dict)?;
    
    Ok(py_result.into())
}

/// Ultra-fast Python wrapper with automatic exhaustive/sampling selection
#[pyfunction]
#[pyo3(signature = (tree1_newick, tree2_newick, number_of_trials=5000, min_triplets_at_depth=1, seed=None))]
fn triplets_correct_ultra(
    py: Python,
    tree1_newick: &str,
    tree2_newick: &str,
    number_of_trials: usize,
    min_triplets_at_depth: usize,
    seed: Option<u64>,
) -> PyResult<PyObject> {
    // Parse trees from Newick strings
    let mut tree1 = PhyloTree::from_newick(tree1_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree1: {:?}", e)))?;
    
    let mut tree2 = PhyloTree::from_newick(tree2_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree2: {:?}", e)))?;
    
    // Calculate triplets correct using ultra-fast implementation
    let result = ultra_fast_triplets_correct(&mut tree1, &mut tree2, number_of_trials, min_triplets_at_depth, seed);
    
    // Convert result to Python dictionaries
    let py_result = PyDict::new_bound(py);
    
    let all_dict = PyDict::new_bound(py);
    for (k, v) in result.all_triplets_correct {
        all_dict.set_item(k, v)?;
    }
    
    let resolvable_dict = PyDict::new_bound(py);
    for (k, v) in result.resolvable_triplets_correct {
        resolvable_dict.set_item(k, v)?;
    }
    
    let unresolved_dict = PyDict::new_bound(py);
    for (k, v) in result.unresolved_triplets_correct {
        unresolved_dict.set_item(k, v)?;
    }
    
    let proportion_dict = PyDict::new_bound(py);
    for (k, v) in result.proportion_unresolvable {
        proportion_dict.set_item(k, v)?;
    }
    
    py_result.set_item("all_triplets_correct", all_dict)?;
    py_result.set_item("resolvable_triplets_correct", resolvable_dict)?;
    py_result.set_item("unresolved_triplets_correct", unresolved_dict)?;
    py_result.set_item("proportion_unresolvable", proportion_dict)?;
    py_result.set_item("method_used", result.method_used)?;
    py_result.set_item("total_triplets_checked", result.total_triplets_checked)?;
    
    Ok(py_result.into())
}

/// Parallel Python wrapper with rayon-based parallelization
#[pyfunction]
#[pyo3(signature = (tree1_newick, tree2_newick, number_of_trials=10000, min_triplets_at_depth=1, seed=None, max_threads=None))]
fn triplets_correct_parallel(
    py: Python,
    tree1_newick: &str,
    tree2_newick: &str,
    number_of_trials: usize,
    min_triplets_at_depth: usize,
    seed: Option<u64>,
    max_threads: Option<usize>,
) -> PyResult<PyObject> {
    // Parse trees from Newick strings
    let mut tree1 = PhyloTree::from_newick(tree1_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree1: {:?}", e)))?;
    
    let mut tree2 = PhyloTree::from_newick(tree2_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree2: {:?}", e)))?;
    
    // Calculate triplets correct using parallel implementation
    let result = parallel_triplets_correct(&mut tree1, &mut tree2, number_of_trials, min_triplets_at_depth, seed, max_threads);
    
    // Convert result to Python dictionaries
    let py_result = PyDict::new_bound(py);
    
    let all_dict = PyDict::new_bound(py);
    for (k, v) in result.all_triplets_correct {
        all_dict.set_item(k, v)?;
    }
    
    let resolvable_dict = PyDict::new_bound(py);
    for (k, v) in result.resolvable_triplets_correct {
        resolvable_dict.set_item(k, v)?;
    }
    
    let unresolved_dict = PyDict::new_bound(py);
    for (k, v) in result.unresolved_triplets_correct {
        unresolved_dict.set_item(k, v)?;
    }
    
    let proportion_dict = PyDict::new_bound(py);
    for (k, v) in result.proportion_unresolvable {
        proportion_dict.set_item(k, v)?;
    }
    
    py_result.set_item("all_triplets_correct", all_dict)?;
    py_result.set_item("resolvable_triplets_correct", resolvable_dict)?;
    py_result.set_item("unresolved_triplets_correct", unresolved_dict)?;
    py_result.set_item("proportion_unresolvable", proportion_dict)?;
    py_result.set_item("method_used", result.method_used)?;
    py_result.set_item("total_triplets_checked", result.total_triplets_checked)?;
    py_result.set_item("parallel_chunks_used", result.parallel_chunks_used)?;
    
    Ok(py_result.into())
}

/// Optimized Python wrapper using phylo-rs inspired techniques
#[pyfunction]
#[pyo3(signature = (tree1_newick, tree2_newick, number_of_trials=10000, min_triplets_at_depth=1, seed=None, max_threads=None))]
fn triplets_correct_optimized(
    py: Python,
    tree1_newick: &str,
    tree2_newick: &str,
    number_of_trials: usize,
    min_triplets_at_depth: usize,
    seed: Option<u64>,
    max_threads: Option<usize>,
) -> PyResult<PyObject> {
    // Parse trees from Newick strings
    let mut tree1 = PhyloTree::from_newick(tree1_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree1: {:?}", e)))?;
    
    let mut tree2 = PhyloTree::from_newick(tree2_newick.as_bytes())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to parse tree2: {:?}", e)))?;
    
    // Calculate triplets correct using optimized implementation
    let result = optimized_triplets_correct(&mut tree1, &mut tree2, number_of_trials, min_triplets_at_depth, seed, max_threads);
    
    // Convert result to Python dictionaries
    let py_result = PyDict::new_bound(py);
    
    let all_dict = PyDict::new_bound(py);
    for (k, v) in result.all_triplets_correct {
        all_dict.set_item(k, v)?;
    }
    
    let resolvable_dict = PyDict::new_bound(py);
    for (k, v) in result.resolvable_triplets_correct {
        resolvable_dict.set_item(k, v)?;
    }
    
    let unresolved_dict = PyDict::new_bound(py);
    for (k, v) in result.unresolved_triplets_correct {
        unresolved_dict.set_item(k, v)?;
    }
    
    let proportion_dict = PyDict::new_bound(py);
    for (k, v) in result.proportion_unresolvable {
        proportion_dict.set_item(k, v)?;
    }
    
    py_result.set_item("all_triplets_correct", all_dict)?;
    py_result.set_item("resolvable_triplets_correct", resolvable_dict)?;
    py_result.set_item("unresolved_triplets_correct", unresolved_dict)?;
    py_result.set_item("proportion_unresolvable", proportion_dict)?;
    py_result.set_item("method_used", result.method_used)?;
    py_result.set_item("total_triplets_checked", result.total_triplets_checked)?;
    py_result.set_item("parallel_chunks_used", result.parallel_chunks_used)?;
    
    Ok(py_result.into())
}


/// Python module definition
#[pymodule]
fn stellars(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(triplets_correct, m)?)?;
    m.add_function(wrap_pyfunction!(triplets_correct_fast, m)?)?;
    m.add_function(wrap_pyfunction!(triplets_correct_ultra, m)?)?;
    m.add_function(wrap_pyfunction!(triplets_correct_parallel, m)?)?;
    m.add_function(wrap_pyfunction!(triplets_correct_optimized, m)?)?;
    Ok(())
}
