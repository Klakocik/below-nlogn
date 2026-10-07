//! Exact convolution modulo the Goldilocks prime `p = 2^64 - 2^32 + 1`.
//!
//! The paper does its transforms with complex arithmetic and then spends a
//! whole section proving the error stays below `1/2` (Section 8). A
//! number-theoretic transform is the exact shortcut: no rounding, no error
//! analysis, no third decimal place to argue about. `7` is a primitive root
//! and `p - 1 = 2^32 * (2^32 - 1)`, so power-of-two lengths up to `2^32` are
//! fair game. With 16-bit digits and at most `2^31` coefficients every
//! coefficient stays below `2^63 < p`: what comes out is the true
//! convolution, not an approximation of it.
//!
//! This is what [`crate::mul`] uses. The floating-point material lives in
//! [`crate::fft`], [`crate::resample`] and [`crate::transforms`], where the
//! paper's own machinery is studied on purpose.

/// The Goldilocks prime `2^64 - 2^32 + 1`.
pub const P: u64 = 0xFFFF_FFFF_0000_0001;

/// A primitive root modulo [`P`].
pub const GENERATOR: u64 = 7;

/// Maximum supported transform length (the 2-adic valuation of `P - 1`).
pub const MAX_LEN: usize = 1 << 32;

fn mul_mod(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % P as u128) as u64
}

fn add_mod(a: u64, b: u64) -> u64 {
    let s = a.wrapping_add(b);
    if s < a {
        // Overflow: a + b - 2^64 + (2^64 - P) = a + b - P, and
        // 2^64 - P = 2^32 - 1.
        s.wrapping_add((1u64 << 32) - 1)
    } else if s >= P {
        s - P
    } else {
        s
    }
}

fn sub_mod(a: u64, b: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        // b - a <= P, so P - (b - a) < P and never overflows.
        P - (b - a)
    }
}

/// Modular exponentiation.
pub fn pow_mod(mut base: u64, mut exp: u64) -> u64 {
    let mut acc = 1u64;
    base %= P;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base);
        }
        base = mul_mod(base, base);
        exp >>= 1;
    }
    acc
}

fn bit_reverse_permute(values: &mut [u64]) {
    let n = values.len();
    if n <= 1 {
        return;
    }
    let bits = n.trailing_zeros();
    for i in 0..n {
        let j = (i as u64).reverse_bits() >> (64 - bits);
        if (j as usize) > i {
            values.swap(i, j as usize);
        }
    }
}

/// In-place number-theoretic transform of a power-of-two length array.
///
/// `inverse = false` computes the forward transform; `inverse = true`
/// computes the inverse and includes the `n^{-1}` normalization.
pub fn transform(values: &mut [u64], inverse: bool) {
    let n = values.len();
    assert!(n.is_power_of_two(), "NTT length must be a power of two");
    assert!(n <= MAX_LEN, "NTT length exceeds the supported 2^32");
    if n <= 1 {
        return;
    }
    bit_reverse_permute(values);
    let mut len = 2usize;
    while len <= n {
        let root = pow_mod(GENERATOR, (P - 1) / len as u64);
        let root = if inverse { pow_mod(root, P - 2) } else { root };
        for block in values.chunks_exact_mut(len) {
            let mut w = 1u64;
            let (left, right) = block.split_at_mut(len / 2);
            for j in 0..len / 2 {
                let u = left[j];
                let v = mul_mod(right[j], w);
                left[j] = add_mod(u, v);
                right[j] = sub_mod(u, v);
                w = mul_mod(w, root);
            }
        }
        len <<= 1;
    }
    if inverse {
        let n_inv = pow_mod(n as u64, P - 2);
        for v in values.iter_mut() {
            *v = mul_mod(*v, n_inv);
        }
    }
}

/// Exact cyclic-free (linear) convolution of nonnegative digit arrays.
///
/// The caller must guarantee that every exact convolution coefficient is
/// `< P`; this is checked against the actual digit maxima and lengths and
/// panics otherwise.
pub fn convolve(a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let out_len = a.len() + b.len() - 1;
    let n = out_len.next_power_of_two();
    assert!(
        n <= MAX_LEN,
        "convolution too large for the supported NTT length"
    );
    let max_a = *a.iter().max().unwrap();
    let max_b = *b.iter().max().unwrap();
    let terms = a.len().min(b.len()) as u128;
    let bound = terms * max_a as u128 * max_b as u128;
    assert!(
        bound < P as u128,
        "coefficient bound {bound} exceeds the prime; use narrower digits or multiple primes"
    );

    let mut fa = vec![0u64; n];
    let mut fb = vec![0u64; n];
    for (i, &v) in a.iter().enumerate() {
        fa[i] = v % P;
    }
    for (i, &v) in b.iter().enumerate() {
        fb[i] = v % P;
    }
    transform(&mut fa, false);
    transform(&mut fb, false);
    for i in 0..n {
        fa[i] = mul_mod(fa[i], fb[i]);
    }
    transform(&mut fa, true);
    fa.truncate(out_len);
    fa
}

#[cfg(test)]
mod tests {
    use super::*;

    fn schoolbook(a: &[u64], b: &[u64]) -> Vec<u64> {
        let mut out = vec![0u64; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        out
    }

    #[test]
    fn ntt_round_trip() {
        let mut v: Vec<u64> = (0..16).map(|i| (i * i + 1) % P).collect();
        let orig = v.clone();
        transform(&mut v, false);
        transform(&mut v, true);
        assert_eq!(v, orig);
    }

    #[test]
    fn convolution_matches_schoolbook() {
        let a: Vec<u64> = (0..37).map(|i| (i * 12345 % 65536) as u64).collect();
        let b: Vec<u64> = (0..41).map(|i| (i * 54321 % 65536) as u64).collect();
        assert_eq!(convolve(&a, &b), schoolbook(&a, &b));
    }

    #[test]
    fn empty_convolution() {
        assert!(convolve(&[], &[1, 2, 3]).is_empty());
    }
}
