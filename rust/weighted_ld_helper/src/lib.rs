use pyo3::prelude::*;
use pyo3::exceptions;

use ndarray::prelude::*;
use ndarray::{s, Array, ArrayViewMut1};
use numpy::{IntoPyArray, PyArray2};

fn seq_to_arr<'a>(seq: &str, mut dest: ArrayViewMut1<'a, u8>) -> Result<(), &'static str> {
    if !seq.is_ascii() {
        return Err("Non-ascii characters in sequence");
    }
    
    for (i, b) in seq.bytes().enumerate() {
        dest[i] = match b {
            b'a' | b'A' => 0,
            b'c' | b'C' => 1,
            b'g' | b'G' => 2,
            b't' | b'T' => 3,
            b'-' => 4,
            _ => 5,
        };
    }
    
    Ok(())
}

#[pymodule]
fn weighted_ld_helper(_py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "multiseq_to_arr")]
    fn multiseq_to_arr_py<'py>(py: Python<'py>, multiseq: PyObject) -> PyResult<&'py PyArray2<u8>> {
        let num_seqs: usize = multiseq
            .call_method0(py, "__len__")?
            .extract(py)?;
        let seq_len: usize = multiseq
            .call_method0(py, "get_alignment_length")?
            .extract(py)?;

        let mut result = Array::<u8, _>::zeros((num_seqs, seq_len).f());

        for seq_index in 0..num_seqs {
            let seq_str_py = multiseq
                .call_method1(py, "__getitem__", (seq_index, ))?
                .getattr(py, "seq")?
                .call_method0(py, "__str__")?;
            let seq_str: &str = seq_str_py.extract(py)?;
            let dest_slice = result.slice_mut(s![seq_index, ..]);

            seq_to_arr(seq_str, dest_slice)
                .map_err(|e| PyErr::new::<exceptions::PyException, _>(e))?
        }
        
        Ok(result.into_pyarray(py))
    }

    Ok(())
}