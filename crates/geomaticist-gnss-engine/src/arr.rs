//! Epoch-major 2-D arrays `(ne, ns)` as the engine's per-epoch, per-satellite tables,
//! for floats, booleans and arc ids.

use serde::{Deserialize, Serialize};

/// A dense `(ne, ns)` array, row-major (epoch-major).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Arr2<T> {
    /// Rows (epochs).
    pub ne: usize,
    /// Columns (satellites).
    pub ns: usize,
    /// Elements, row-major.
    pub data: Vec<T>,
}

impl<T: Copy> Arr2<T> {
    /// An array filled with `v`.
    pub fn filled(ne: usize, ns: usize, v: T) -> Self {
        Self {
            ne,
            ns,
            data: vec![v; ne * ns],
        }
    }

    /// Element `(e, j)`.
    #[inline]
    #[must_use]
    pub fn at(&self, e: usize, j: usize) -> T {
        self.data[e * self.ns + j]
    }

    /// Mutable element `(e, j)`.
    #[inline]
    pub fn set(&mut self, e: usize, j: usize, v: T) {
        self.data[e * self.ns + j] = v;
    }

    /// Row `e` as a slice.
    #[inline]
    #[must_use]
    pub fn row(&self, e: usize) -> &[T] {
        &self.data[e * self.ns..(e + 1) * self.ns]
    }

    /// Row `e` as a mutable slice.
    #[inline]
    pub fn row_mut(&mut self, e: usize) -> &mut [T] {
        &mut self.data[e * self.ns..(e + 1) * self.ns]
    }

    /// Column `j` as a new vector.
    #[must_use]
    pub fn col(&self, j: usize) -> Vec<T> {
        (0..self.ne).map(|e| self.at(e, j)).collect()
    }

    /// Element-wise map.
    pub fn map<U: Copy>(&self, f: impl Fn(T) -> U) -> Arr2<U> {
        Arr2 {
            ne: self.ne,
            ns: self.ns,
            data: self.data.iter().map(|&x| f(x)).collect(),
        }
    }

    /// Element-wise combination with another array of the same shape.
    pub fn zip_map<U: Copy, V: Copy>(&self, o: &Arr2<U>, f: impl Fn(T, U) -> V) -> Arr2<V> {
        debug_assert_eq!((self.ne, self.ns), (o.ne, o.ns));
        Arr2 {
            ne: self.ne,
            ns: self.ns,
            data: self
                .data
                .iter()
                .zip(&o.data)
                .map(|(&a, &b)| f(a, b))
                .collect(),
        }
    }
}

/// Floats, NaN-filled by default.
pub type F2 = Arr2<f64>;
/// Booleans.
pub type B2 = Arr2<bool>;
/// Arc / segment ids, -1 where none.
pub type I2 = Arr2<i64>;

impl F2 {
    /// NaN-filled.
    #[must_use]
    pub fn nans(ne: usize, ns: usize) -> Self {
        Self::filled(ne, ns, f64::NAN)
    }

    /// Zero-filled.
    #[must_use]
    pub fn zeros(ne: usize, ns: usize) -> Self {
        Self::filled(ne, ns, 0.0)
    }
}

impl B2 {
    /// All false.
    #[must_use]
    pub fn falses(ne: usize, ns: usize) -> Self {
        Self::filled(ne, ns, false)
    }

    /// Number of true elements.
    #[must_use]
    pub fn count(&self) -> usize {
        self.data.iter().filter(|&&b| b).count()
    }

    /// Whether any element is true.
    #[must_use]
    pub fn any(&self) -> bool {
        self.data.iter().any(|&b| b)
    }

    /// Whether any element of column `j` is true.
    #[must_use]
    pub fn any_col(&self, j: usize) -> bool {
        (0..self.ne).any(|e| self.at(e, j))
    }

    /// `a & b`.
    #[must_use]
    pub fn and(&self, o: &B2) -> B2 {
        self.zip_map(o, |a, b| a && b)
    }

    /// `a | b`.
    #[must_use]
    pub fn or(&self, o: &B2) -> B2 {
        self.zip_map(o, |a, b| a || b)
    }
}

/// A `(ne, ns, 3)` array of vectors.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Vec3Arr {
    /// Rows.
    pub ne: usize,
    /// Columns.
    pub ns: usize,
    /// Elements.
    pub data: Vec<[f64; 3]>,
}

impl Vec3Arr {
    /// Zero-filled.
    #[must_use]
    pub fn zeros(ne: usize, ns: usize) -> Self {
        Self {
            ne,
            ns,
            data: vec![[0.0; 3]; ne * ns],
        }
    }

    /// Element `(e, j)`.
    #[inline]
    #[must_use]
    pub fn at(&self, e: usize, j: usize) -> [f64; 3] {
        self.data[e * self.ns + j]
    }

    /// Set element `(e, j)`.
    #[inline]
    pub fn set(&mut self, e: usize, j: usize, v: [f64; 3]) {
        self.data[e * self.ns + j] = v;
    }

    /// Flat f64 view (row-major, 3 per element).
    #[must_use]
    pub fn flat(&self) -> Vec<f64> {
        self.data.iter().flatten().copied().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexing() {
        let mut a = F2::zeros(2, 3);
        a.set(1, 2, 5.0);
        assert_eq!(a.at(1, 2), 5.0);
        assert_eq!(a.row(1), &[0.0, 0.0, 5.0]);
        assert_eq!(a.col(2), vec![0.0, 5.0]);
        let b = a.map(|x| x > 1.0);
        assert_eq!(b.count(), 1);
        assert!(b.any_col(2) && !b.any_col(0));
    }
}
