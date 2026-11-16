//! Dense linear algebra for the normal equations: a row-major symmetric matrix the assembly adds into, with
//! Cholesky solve and inverse through `faer` (LU with partial pivoting as the fallback where the matrix is
//! not positive definite, which is what `np.linalg.solve` / `inv` would do anyway), and the small 3 × 3
//! helpers the checks use.
//!
//! Against an LU solve (`np.linalg.solve` / `inv`, LAPACK) Cholesky gives the same math, with
//! last-bit differences only.
//!
//! The Cholesky path runs faer **sequentially** (`Par::Seq`, through faer's in-place routines rather than
//! its solver objects, which take the process-wide parallelism setting): the matrices here are at most a
//! few thousand square, where one core factors them in tens of milliseconds, and a baseline job runs
//! several solves at once on one rayon pool — faer's fork-join storms inside a shared, busy pool cost
//! far more than they gain. It also makes the numbers independent of the thread count (with faer's
//! parallel kernels the last bits depended on the number of threads).

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::cholesky::llt::factor::LltRegularization;
use faer::linalg::cholesky::llt::{factor, inverse as llt_inverse, solve as llt_solve};
use faer::linalg::solvers::{DenseSolveCore, PartialPivLu, Solve};
use faer::{Mat, Par, Spec};

/// Tile edge of the row-major <-> column-major copies.
const TILE: usize = 64;

/// A dense square matrix, row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct DMat {
    /// Order.
    pub n: usize,
    /// Elements, row-major.
    pub data: Vec<f64>,
}

impl DMat {
    /// Zero matrix.
    #[must_use]
    pub fn zeros(n: usize) -> Self {
        Self {
            n,
            data: vec![0.0; n * n],
        }
    }

    /// Identity.
    #[must_use]
    pub fn eye(n: usize) -> Self {
        let mut m = Self::zeros(n);
        for i in 0..n {
            m.data[i * n + i] = 1.0;
        }
        m
    }

    /// Element `(i, j)`.
    #[inline]
    #[must_use]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    /// Mutable element `(i, j)`.
    #[inline]
    pub fn at_mut(&mut self, i: usize, j: usize) -> &mut f64 {
        &mut self.data[i * self.n + j]
    }

    /// Add `v` to `(i, j)`.
    #[inline]
    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] += v;
    }

    /// Row `i`.
    #[inline]
    #[must_use]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.n..(i + 1) * self.n]
    }

    /// The sub-matrix on `idx × idx`.
    #[must_use]
    pub fn sub(&self, idx: &[usize]) -> DMat {
        let k = idx.len();
        let mut m = DMat::zeros(k);
        for (a, &i) in idx.iter().enumerate() {
            for (b, &j) in idx.iter().enumerate() {
                m.data[a * k + b] = self.at(i, j);
            }
        }
        m
    }

    /// `M v`.
    #[must_use]
    pub fn matvec(&self, v: &[f64]) -> Vec<f64> {
        (0..self.n).map(|i| dot(self.row(i), v)).collect()
    }

    /// `self * k` (every element).
    #[must_use]
    pub fn scaled(&self, k: f64) -> DMat {
        DMat {
            n: self.n,
            data: self.data.iter().map(|x| x * k).collect(),
        }
    }

    /// Copy into faer's column-major matrix, by tiles (both sides stay in cache).
    fn to_faer(&self) -> Mat<f64> {
        let n = self.n;
        let mut a = Mat::<f64>::zeros(n, n);
        for jb in (0..n).step_by(TILE) {
            let je = (jb + TILE).min(n);
            for ib in (0..n).step_by(TILE) {
                for i in ib..(ib + TILE).min(n) {
                    let row = &self.data[i * n..(i + 1) * n];
                    for j in jb..je {
                        a[(i, j)] = row[j];
                    }
                }
            }
        }
        a
    }

    /// Copy back from faer's column-major matrix, by tiles.
    fn from_faer(m: &Mat<f64>) -> DMat {
        let n = m.nrows();
        let mut data = vec![0.0; n * n];
        for jb in (0..n).step_by(TILE) {
            let je = (jb + TILE).min(n);
            for ib in (0..n).step_by(TILE) {
                for i in ib..(ib + TILE).min(n) {
                    let row = &mut data[i * n..(i + 1) * n];
                    for j in jb..je {
                        row[j] = m[(i, j)];
                    }
                }
            }
        }
        DMat { n, data }
    }

    /// The lower triangle copied into faer's column-major matrix (the rest zero), by tiles.
    fn to_faer_lower(&self) -> Mat<f64> {
        let n = self.n;
        let mut a = Mat::<f64>::zeros(n, n);
        for jb in (0..n).step_by(TILE) {
            let je = (jb + TILE).min(n);
            for ib in (jb..n).step_by(TILE) {
                for i in ib..(ib + TILE).min(n) {
                    let row = &self.data[i * n..(i + 1) * n];
                    for j in jb..je.min(i + 1) {
                        a[(i, j)] = row[j];
                    }
                }
            }
        }
        a
    }

    /// Symmetric copy back from a faer matrix whose lower triangle holds the result, by tiles.
    fn from_faer_lower(m: &Mat<f64>) -> DMat {
        let n = m.nrows();
        let mut data = vec![0.0; n * n];
        for jb in (0..n).step_by(TILE) {
            let je = (jb + TILE).min(n);
            for ib in (jb..n).step_by(TILE) {
                for i in ib..(ib + TILE).min(n) {
                    for j in jb..je.min(i + 1) {
                        let v = m[(i, j)];
                        data[i * n + j] = v;
                        data[j * n + i] = v;
                    }
                }
            }
        }
        DMat { n, data }
    }

    /// The Cholesky factor `L` (lower triangle; the rest zero), sequential; `None` when the matrix is not
    /// positive definite.
    fn llt(&self) -> Option<Mat<f64>> {
        let n = self.n;
        let mut l = self.to_faer_lower();
        let mut mem = MemBuffer::new(factor::cholesky_in_place_scratch::<f64>(
            n,
            Par::Seq,
            Spec::default(),
        ));
        factor::cholesky_in_place(
            l.as_mut(),
            LltRegularization::default(),
            Par::Seq,
            MemStack::new(&mut mem),
            Spec::default(),
        )
        .ok()?;
        for j in 0..n {
            for i in 0..j {
                l[(i, j)] = 0.0;
            }
        }
        Some(l)
    }

    /// `L L^T X = B` in place (sequential).
    fn llt_solve(l: &Mat<f64>, rhs: &mut Mat<f64>) {
        let mut mem = MemBuffer::new(llt_solve::solve_in_place_scratch::<f64>(
            l.nrows(),
            rhs.ncols(),
            Par::Seq,
        ));
        llt_solve::solve_in_place(l.as_ref(), rhs.as_mut(), Par::Seq, MemStack::new(&mut mem));
    }

    /// Solve `M x = b` (Cholesky; LU fallback). `None` when singular.
    #[must_use]
    pub fn solve(&self, b: &[f64]) -> Option<Vec<f64>> {
        let mut rhs = Mat::from_fn(self.n, 1, |i, _| b[i]);
        match self.llt() {
            Some(l) => Self::llt_solve(&l, &mut rhs),
            None => {
                let a = self.to_faer();
                rhs = PartialPivLu::new(a.as_ref()).solve(&rhs);
            }
        }
        let out: Vec<f64> = (0..self.n).map(|i| rhs[(i, 0)]).collect();
        out.iter().all(|v| v.is_finite()).then_some(out)
    }

    /// Solve `M X = B` for several right-hand sides given as columns (`cols[k]` is column k).
    #[must_use]
    pub fn solve_many(&self, cols: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
        if cols.is_empty() {
            return Some(Vec::new());
        }
        let mut rhs = Mat::from_fn(self.n, cols.len(), |i, k| cols[k][i]);
        match self.llt() {
            Some(l) => Self::llt_solve(&l, &mut rhs),
            None => {
                let a = self.to_faer();
                rhs = PartialPivLu::new(a.as_ref()).solve(&rhs);
            }
        }
        let out: Vec<Vec<f64>> = (0..cols.len())
            .map(|k| (0..self.n).map(|i| rhs[(i, k)]).collect())
            .collect();
        out.iter()
            .all(|c| c.iter().all(|v| v.is_finite()))
            .then_some(out)
    }

    /// The inverse (Cholesky; LU fallback). `None` when singular.
    #[must_use]
    pub fn inverse(&self) -> Option<DMat> {
        let m = match self.llt() {
            Some(l) => {
                let n = self.n;
                let mut out = Mat::<f64>::zeros(n, n);
                let mut mem = MemBuffer::new(llt_inverse::inverse_scratch::<f64>(n, Par::Seq));
                llt_inverse::inverse(out.as_mut(), l.as_ref(), Par::Seq, MemStack::new(&mut mem));
                DMat::from_faer_lower(&out)
            }
            None => DMat::from_faer(&PartialPivLu::new(self.to_faer().as_ref()).inverse()),
        };
        m.data.iter().all(|v| v.is_finite()).then_some(m)
    }

    /// Symmetrise in place: `(M + Mᵀ) / 2` is not what NumPy does; the inverse of a symmetric matrix from
    /// Cholesky is symmetric to the last bit already, so this only copies the lower triangle up where the
    /// assembly filled one triangle.
    pub fn mirror_lower(&mut self) {
        for i in 0..self.n {
            for j in 0..i {
                let v = self.at(i, j);
                *self.at_mut(j, i) = v;
            }
        }
    }
}

/// Dot product (sequential, as a BLAS ddot of short vectors).
#[inline]
#[must_use]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// `a - b` for 3-vectors.
#[inline]
#[must_use]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// `a + b` for 3-vectors.
#[inline]
#[must_use]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// `k a` for 3-vectors.
#[inline]
#[must_use]
pub fn mul3(a: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] * k, a[1] * k, a[2] * k]
}

/// Euclidean norm of a 3-vector as `np.linalg.norm` computes it (`sqrt(x² + y² + z²)`, left to right).
#[inline]
#[must_use]
pub fn norm3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Dot of 3-vectors.
#[inline]
#[must_use]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline]
#[must_use]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `R v` for a 3 × 3 matrix in rows.
#[inline]
#[must_use]
pub fn mat3_vec(r: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [dot3(r[0], v), dot3(r[1], v), dot3(r[2], v)]
}

/// `Rᵀ v`.
#[inline]
#[must_use]
pub fn mat3t_vec(r: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        r[0][0] * v[0] + r[1][0] * v[1] + r[2][0] * v[2],
        r[0][1] * v[0] + r[1][1] * v[1] + r[2][1] * v[2],
        r[0][2] * v[0] + r[1][2] * v[1] + r[2][2] * v[2],
    ]
}

/// `A B` for 3 × 3.
#[must_use]
pub fn mat3_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    out
}

/// `Aᵀ`.
#[must_use]
pub fn mat3_t(a: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [a[0][0], a[1][0], a[2][0]],
        [a[0][1], a[1][1], a[2][1]],
        [a[0][2], a[1][2], a[2][2]],
    ]
}

/// `R C Rᵀ` for 3 × 3 (a covariance rotated).
#[must_use]
pub fn rotate_cov3(r: &[[f64; 3]; 3], c: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    mat3_mul(&mat3_mul(r, c), &mat3_t(r))
}

/// Inverse of a 3 × 3 (`None` when singular).
#[must_use]
pub fn inv3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let d = DMat {
        n: 3,
        data: m.iter().flatten().copied().collect(),
    };
    let inv = d.inverse()?;
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = inv.at(i, j);
        }
    }
    Some(out)
}

/// Solve a 3 × 3 system.
#[must_use]
pub fn solve3(m: &[[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let d = DMat {
        n: 3,
        data: m.iter().flatten().copied().collect(),
    };
    let x = d.solve(&b)?;
    Some([x[0], x[1], x[2]])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_and_solve_over_several_tiles() {
        // n not a multiple of the tile: the tiled copies cover the edges
        let n = 2 * TILE + 7;
        let mut m = DMat::zeros(n);
        for i in 0..n {
            for j in 0..n {
                let d = (i as f64 - j as f64).abs();
                *m.at_mut(i, j) = 1.0 / (1.0 + d) + if i == j { 2.0 } else { 0.0 };
            }
        }
        let inv = m.inverse().unwrap();
        for i in 0..n {
            for j in 0..n {
                let col: Vec<f64> = (0..n).map(|k| inv.at(k, j)).collect();
                let s = dot(m.row(i), &col);
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((s - want).abs() < 1e-10, "({i}, {j}) {s}");
            }
        }
        let x0: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let x = m.solve(&m.matvec(&x0)).unwrap();
        assert!(x.iter().zip(&x0).all(|(a, b)| (a - b).abs() < 1e-8));
    }

    #[test]
    fn cholesky_solve_and_inverse() {
        let m = DMat {
            n: 3,
            data: vec![4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0],
        };
        let x = m.solve(&[1.0, 2.0, 3.0]).unwrap();
        let back = m.matvec(&x);
        assert!((back[0] - 1.0).abs() < 1e-12 && (back[2] - 3.0).abs() < 1e-12);
        let inv = m.inverse().unwrap();
        let id = DMat {
            n: 3,
            data: (0..9)
                .map(|k| {
                    dot(
                        m.row(k / 3),
                        &inv.data.chunks(3).map(|r| r[k % 3]).collect::<Vec<_>>(),
                    )
                })
                .collect(),
        };
        for i in 0..3 {
            for j in 0..3 {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((id.at(i, j) - want).abs() < 1e-12);
            }
        }
        // indefinite: LU fallback
        let m2 = DMat {
            n: 2,
            data: vec![0.0, 1.0, 1.0, 0.0],
        };
        assert_eq!(m2.solve(&[2.0, 3.0]).unwrap(), vec![3.0, 2.0]);
        assert!(inv3(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).is_some());
    }
}
