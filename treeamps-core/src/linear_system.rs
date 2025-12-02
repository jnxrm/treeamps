use nalgebra::{DMatrix, SVD};

/// Wrapper around a dense linear system matrix.
#[derive(Debug, Clone)]
pub struct LinearSystem {
    pub matrix: DMatrix<f64>,
}

impl LinearSystem {
    /// Create a linear system from a matrix (rows x cols).
    pub fn new(rows: Vec<Vec<f64>>) -> Self {
        if rows.is_empty() {
            return Self {
                matrix: DMatrix::from_row_slice(0, 0, &[]),
            };
        }
        let n_rows = rows.len();
        let n_cols = rows[0].len();
        let flat: Vec<f64> = rows.into_iter().flatten().collect();
        Self {
            matrix: DMatrix::from_row_slice(n_rows, n_cols, &flat),
        }
    }

    /// Compute the nullspace (right nullspace) of the matrix using SVD.
    /// Returns vectors that span the nullspace.
    pub fn nullspace(&self, tol: f64) -> Vec<Vec<f64>> {
        if self.matrix.nrows() == 0 || self.matrix.ncols() == 0 {
            return vec![];
        }

        let svd = SVD::new(self.matrix.clone(), true, true);
        
        // V^T is stored; we want the columns of V corresponding to small singular values
        let singular_values = svd.singular_values;
        let v_t = svd.v_t.expect("SVD did not compute V^T");

        // Find indices where singular value < tol
        let mut null_indices = Vec::new();
        for (i, &sigma) in singular_values.iter().enumerate() {
            if sigma < tol {
                null_indices.push(i);
            }
        }

        // Extract corresponding rows of V^T (= columns of V)
        let mut nullspace_vectors = Vec::new();
        for idx in null_indices {
            let row = v_t.row(idx);
            nullspace_vectors.push(row.iter().copied().collect());
        }

        nullspace_vectors
    }

    /// Return the dimensions (rows, cols) of the matrix.
    pub fn dimensions(&self) -> (usize, usize) {
        (self.matrix.nrows(), self.matrix.ncols())
    }

    /// Compute the rank of the matrix (number of singular values above tolerance).
    pub fn rank(&self, tol: f64) -> usize {
        if self.matrix.nrows() == 0 || self.matrix.ncols() == 0 {
            return 0;
        }
        let svd = SVD::new(self.matrix.clone(), false, false);
        svd.singular_values.iter().filter(|&&s| s >= tol).count()
    }
}
