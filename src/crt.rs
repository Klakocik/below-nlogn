//! Good–Thomas CRT layout: one long cyclic convolution becomes a grid of
//! short ones (Section 8's explicit layout).
//!
//! For pairwise-coprime `n_1, ..., n_k` with `N = prod n_i`, the ring `Z/NZ`
//! is `prod Z/n_i Z` under `j -> (j mod n_1, ..., j mod n_k)`, and this map
//! carries cyclic convolution to the multidimensional cyclic convolution on
//! the product of axes. The paper puts its radix-digit coefficient list on
//! prime axes by exactly this map. Residues are free; convolution is carried.

/// A flat index <-> multidimensional residue layout.
#[derive(Clone, Debug)]
pub struct CrtLayout {
    pub dims: Vec<usize>,
    /// Product of the dimensions.
    pub len: usize,
    /// `flat_to_tuple[t]` is the residue tuple of `t`.
    pub flat_to_tuple: Vec<Vec<usize>>,
    /// `tuple_to_flat` maps a lexicographic tuple index to the flat index.
    pub tuple_to_flat: Vec<usize>,
}

impl CrtLayout {
    pub fn new(dims: &[usize]) -> Self {
        assert!(!dims.is_empty());
        for d in dims {
            assert!(*d >= 1);
        }
        let len: usize = dims.iter().product();
        let mut flat_to_tuple = Vec::with_capacity(len);
        let mut tuple_to_flat = vec![0usize; len];
        for t in 0..len {
            let tuple: Vec<usize> = dims.iter().map(|&d| t % d).collect();
            let lex = tuple_index(dims, &tuple);
            tuple_to_flat[lex] = t;
            flat_to_tuple.push(tuple);
        }
        CrtLayout {
            dims: dims.to_vec(),
            len,
            flat_to_tuple,
            tuple_to_flat,
        }
    }

    fn assert_coprime(&self) {
        for i in 0..self.dims.len() {
            for j in i + 1..self.dims.len() {
                assert_eq!(
                    gcd(self.dims[i], self.dims[j]),
                    1,
                    "dimensions must be pairwise coprime"
                );
            }
        }
    }

    /// Lexicographic tuple index with the first coordinate most significant.
    fn lex_index(&self, tuple: &[usize]) -> usize {
        tuple_index(&self.dims, tuple)
    }

    /// Transports a flat array into multidimensional (lexicographic) layout.
    pub fn to_tensor<T: Clone>(&self, flat: &[T]) -> Vec<T> {
        self.assert_coprime();
        assert_eq!(flat.len(), self.len);
        let mut out = vec![flat[0].clone(); self.len];
        for t in 0..self.len {
            out[self.lex_index(&self.flat_to_tuple[t])] = flat[t].clone();
        }
        out
    }

    /// Transports multidimensional layout back to a flat array.
    pub fn to_flat<T: Clone>(&self, tensor: &[T]) -> Vec<T> {
        self.assert_coprime();
        assert_eq!(tensor.len(), self.len);
        let mut out = vec![tensor[0].clone(); self.len];
        for t in 0..self.len {
            out[t] = tensor[self.lex_index(&self.flat_to_tuple[t])].clone();
        }
        out
    }

    /// Reconstructs `x mod N` from residues via the CRT.
    pub fn reconstruct(&self, residues: &[usize]) -> usize {
        self.assert_coprime();
        let n = self.len as u128;
        let mut acc = 0u128;
        for (i, &d) in self.dims.iter().enumerate() {
            let m = n / d as u128;
            let inv = mod_inverse((m % d as u128) as u64, d as u64) as u128;
            acc = (acc + (residues[i] as u128 % d as u128) * m % n * inv) % n;
        }
        acc as usize
    }
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn tuple_index(dims: &[usize], tuple: &[usize]) -> usize {
    let mut acc = 0usize;
    for (t, &d) in tuple.iter().zip(dims) {
        acc = acc * d + (t % d);
    }
    acc
}

fn mod_inverse(a: u64, m: u64) -> u64 {
    // Extended Euclid; m small and prime to a.
    let (mut t, mut new_t) = (0i128, 1i128);
    let (mut r, mut new_r) = (m as i128, a as i128);
    while new_r != 0 {
        let q = r / new_r;
        (t, new_t) = (new_t, t - q * new_t);
        (r, new_r) = (new_r, r - q * new_r);
    }
    assert!(r == 1, "not invertible");
    t.rem_euclid(m as i128) as u64
}

/// Naive cyclic convolution of length `n`.
pub fn cyclic_convolution_naive(a: &[u64], b: &[u64], n: usize) -> Vec<u64> {
    assert_eq!(a.len(), n);
    assert_eq!(b.len(), n);
    let mut out = vec![0u128; n];
    for i in 0..n {
        for j in 0..n {
            out[(i + j) % n] += a[i] as u128 * b[j] as u128;
        }
    }
    out.into_iter().map(|v| v as u64).collect()
}

/// Multidimensional cyclic convolution on the product of the given dims,
/// stored lexicographically (first coordinate most significant).
pub fn multidimensional_cyclic_convolution(dims: &[usize], a: &[u64], b: &[u64]) -> Vec<u64> {
    let len: usize = dims.iter().product();
    assert_eq!(a.len(), len);
    assert_eq!(b.len(), len);
    let mut out = vec![0u128; len];
    // Coordinates of each lexicographic index.
    let coords: Vec<Vec<usize>> = (0..len).map(|i| lex_tuple(dims, i)).collect();
    for i in 0..len {
        for j in 0..len {
            let sum: Vec<usize> = coords[i]
                .iter()
                .zip(&coords[j])
                .zip(dims)
                .map(|((x, y), d)| (x + y) % d)
                .collect();
            let dest = tuple_index(dims, &sum);
            out[dest] += a[i] as u128 * b[j] as u128;
        }
    }
    out.into_iter().map(|v| v as u64).collect()
}

fn lex_tuple(dims: &[usize], mut index: usize) -> Vec<usize> {
    let mut coords = vec![0usize; dims.len()];
    for i in (0..dims.len()).rev() {
        coords[i] = index % dims[i];
        index /= dims[i];
    }
    coords
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crt_map_is_bijective_and_additive() {
        let dims = [3usize, 5, 7];
        let layout = CrtLayout::new(&dims);
        assert_eq!(layout.len, 105);
        let mut seen = [false; 105];
        for t in 0..105 {
            seen[layout.reconstruct(&layout.flat_to_tuple[t])] = true;
            let m = (t + 50) % 105;
            for (i, &d) in dims.iter().enumerate() {
                assert_eq!(
                    layout.flat_to_tuple[m][i],
                    (layout.flat_to_tuple[t][i] + layout.flat_to_tuple[50][i]) % d
                );
            }
        }
        assert!(seen.iter().all(|&s| s));
    }

    #[test]
    fn crt_layout_preserves_cyclic_convolution() {
        let dims = [4usize, 5, 7];
        let n = 140;
        let layout = CrtLayout::new(&dims);
        let mut state = 0xDEADBEEFu64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 33) % 100
        };
        let a: Vec<u64> = (0..n).map(|_| next()).collect();
        let b: Vec<u64> = (0..n).map(|_| next()).collect();
        let expected = cyclic_convolution_naive(&a, &b, n);
        let ta = layout.to_tensor(&a);
        let tb = layout.to_tensor(&b);
        let md = multidimensional_cyclic_convolution(&dims, &ta, &tb);
        let got = layout.to_flat(&md);
        assert_eq!(got, expected);
    }
}
