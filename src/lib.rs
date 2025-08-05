use pyo3::prelude::*;

#[pyfunction]
fn hello_stella() -> PyResult<String> {
    Ok("Hello from stella-rs! This is a stub implementation.".to_string())
}

#[pyfunction]
fn version() -> PyResult<String> {
    Ok("0.1.0".to_string())
}

#[pymodule]
fn stellars(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hello_stella, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    Ok(())
}