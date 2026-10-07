//! Section 6's synthetic transforms over `R_r = C[y]/(y^r + 1)`,
//! implemented literally.
//!
//! Roots of unity are powers of `y`, so multiplying by one is a signed
//! permutation of a coefficient list; the DIF recursion uses the butterfly
//! `H_0 = 1/2 [[1, 1], [1, -1]]` and drops its outputs into bit-reversed
//! slots, because that is where the next butterfly wants them. Every step is
//! checked against the naive definition in the tests. Tape-cost accounting
//! is not here; see `docs/paper-map.md`.

use crate::fft::Complex;

/// An element of `C[y]/(y^r + 1)`, stored as coefficients of degree `< r`.
pub type RingElem = Vec<Complex>;

/// Multiplies a ring element by `y^a` (an exact signed permutation of
/// coefficients), which is how the paper's roots of unity act.
pub fn mul_y_pow(a: i64, x: &[Complex], r: usize) -> RingElem {
    assert!(r >= 1);
    let a = a.rem_euclid(2 * r as i64) as usize;
    let mut out = vec![Complex::default(); r];
    for (k, &xk) in x.iter().enumerate() {
        let m = k + a;
        if m < r {
            out[m] += xk;
        } else if m < 2 * r {
            // y^{k + a} = y^r * y^{k + a - r} = -y^{k + a - r}
            out[m - r] -= xk;
        } else {
            // y^{k + a} = y^{2r} * y^{k + a - 2r} = y^{k + a - 2r}
            out[m - 2 * r] += xk;
        }
    }
    out
}

fn add_assign(dst: &mut [Complex], src: &[Complex]) {
    for (d, s) in dst.iter_mut().zip(src) {
        *d += *s;
    }
}

/// Naive `F_t^-` straight from the definition:
/// `(F_t^- z)_j = (1/t) sum_k omega_t^{-jk} z_k` with `omega_t = y^{2r/t}`.
pub fn naive_forward(z: &[RingElem], r: usize) -> Vec<RingElem> {
    let t = z.len();
    assert!(t >= 1 && 2 * r % t == 0, "require t | 2r");
    let step = (2 * r / t) as i64;
    let mut out = vec![vec![Complex::default(); r]; t];
    for (j, slot) in out.iter_mut().enumerate() {
        for (k, zk) in z.iter().enumerate() {
            let exp = -(j as i64) * (k as i64) * step;
            let term = mul_y_pow(exp, zk, r);
            add_assign(slot, &term);
        }
        for c in slot.iter_mut() {
            *c = *c / t as f64;
        }
    }
    out
}

/// Naive `F_t^+`: `(F_t^+ z)_j = (1/t) sum_k omega_t^{jk} z_k`.
pub fn naive_inverse(z: &[RingElem], r: usize) -> Vec<RingElem> {
    let t = z.len();
    assert!(t >= 1 && 2 * r % t == 0, "require t | 2r");
    let step = (2 * r / t) as i64;
    let mut out = vec![vec![Complex::default(); r]; t];
    for (j, slot) in out.iter_mut().enumerate() {
        for (k, zk) in z.iter().enumerate() {
            let exp = (j as i64) * (k as i64) * step;
            let term = mul_y_pow(exp, zk, r);
            add_assign(slot, &term);
        }
        for c in slot.iter_mut() {
            *c = *c / t as f64;
        }
    }
    out
}

/// The decimation-in-frequency recursion of the paper. `inverse = false`
/// computes `F_t^-`; `inverse = true` computes the reverse procedure.
///
/// The length of `z` is the transform length `t` (a power of two dividing
/// `2r`). Outputs are stored with frequency `j` in slot `rev_a(j)`; call
/// [`bit_reverse_permute`] to obtain natural frequency order.
pub fn dif_in_place(z: &mut [RingElem], r: usize, inverse: bool) {
    let t = z.len();
    if t <= 1 {
        return;
    }
    assert!(t.is_power_of_two(), "t must be a power of two");
    assert!(2 * r % t == 0, "require t | 2r");
    let half = t / 2;
    let step = (2 * r / t) as i64;
    let (lo, hi) = z.split_at_mut(half);
    if inverse {
        dif_in_place(lo, r, true);
        dif_in_place(hi, r, true);
        for k in 0..half {
            let exp = (k as i64) * step;
            hi[k] = mul_y_pow(exp, &hi[k], r);
        }
        // Finally apply H_0 (the reverse level, per the paper).
        for k in 0..half {
            let a = lo[k].clone();
            let b = hi[k].clone();
            for c in 0..r {
                lo[k][c] = (a[c] + b[c]) / 2.0;
                hi[k][c] = (a[c] - b[c]) / 2.0;
            }
        }
    } else {
        for k in 0..half {
            let a = lo[k].clone();
            let b = hi[k].clone();
            for c in 0..r {
                lo[k][c] = (a[c] + b[c]) / 2.0;
                hi[k][c] = (a[c] - b[c]) / 2.0;
            }
        }
        for k in 0..half {
            let exp = -(k as i64) * step;
            hi[k] = mul_y_pow(exp, &hi[k], r);
        }
    }
    if !inverse {
        dif_in_place(lo, r, false);
        dif_in_place(hi, r, false);
    }
}

/// Permutes slots by bit reversal (an involution): moves the DIF's stored
/// order into natural frequency order.
pub fn bit_reverse_permute(z: &mut [RingElem]) {
    let n = z.len();
    if n <= 1 {
        return;
    }
    let bits = n.trailing_zeros();
    for i in 0..n {
        let j = (i as u64).reverse_bits() >> (64 - bits);
        if (j as usize) > i {
            z.swap(i, j as usize);
        }
    }
}

/// `F_t^-` in natural frequency order.
pub fn forward(z: &[RingElem], r: usize) -> Vec<RingElem> {
    let mut out = z.to_vec();
    dif_in_place(&mut out, r, false);
    bit_reverse_permute(&mut out);
    out
}

/// `F_t^+` in natural frequency order.
///
/// The reverse procedure consumes the bit-reversed stored order, so we first
/// apply the (involutive) bit-reversal permutation and then run the recursion
/// in reverse mode.
pub fn inverse(z: &[RingElem], r: usize) -> Vec<RingElem> {
    let mut out = z.to_vec();
    bit_reverse_permute(&mut out);
    dif_in_place(&mut out, r, true);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_ring(r: usize, seed: &mut u64) -> RingElem {
        (0..r)
            .map(|_| {
                *seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let x = (*seed >> 11) as f64 / (1u64 << 53) as f64;
                *seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let y = (*seed >> 11) as f64 / (1u64 << 53) as f64;
                Complex::new(2.0 * x - 1.0, 2.0 * y - 1.0)
            })
            .collect()
    }

    fn max_diff(a: &[RingElem], b: &[RingElem]) -> f64 {
        let mut m: f64 = 0.0;
        for (x, y) in a.iter().zip(b) {
            for (p, q) in x.iter().zip(y) {
                m = m.max((*p - *q).abs());
            }
        }
        m
    }

    #[test]
    fn monomial_multiplication_is_a_signed_permutation() {
        let r = 8;
        let mut seed = 7u64;
        let x = random_ring(r, &mut seed);
        let prod = mul_y_pow(3, &x, r);
        // x -> y^3 x, then multiply by y^5 should give y^8 x = -x.
        let round = mul_y_pow(5, &prod, r);
        for (a, b) in round.iter().zip(&x) {
            assert!((*a + *b).abs() < 1e-12);
        }
    }

    #[test]
    fn character_cancellation() {
        // Lemma 6.1: sum_k omega_t^{jk} is t for t | j and 0 otherwise.
        let (r, t) = (8usize, 4usize);
        let step = (2 * r / t) as i64;
        for j in 0..8i64 {
            let mut acc = vec![Complex::default(); r];
            for k in 0..t {
                let term = mul_y_pow(
                    j * k as i64 * step,
                    &{
                        let mut e = vec![Complex::default(); r];
                        e[0] = Complex::new(1.0, 0.0);
                        e
                    },
                    r,
                );
                add_assign(&mut acc, &term);
            }
            if j % t as i64 == 0 {
                assert!((acc[0] - Complex::new(t as f64, 0.0)).abs() < 1e-12);
            } else {
                assert!(acc.iter().all(|c| c.abs() < 1e-12), "j={j}");
            }
        }
    }

    #[test]
    fn dif_matches_naive_and_is_stored_bit_reversed() {
        for (r, t) in [(2usize, 2usize), (4, 2), (4, 4), (8, 4), (8, 8)] {
            let mut seed = 1234u64;
            let z: Vec<RingElem> = (0..t).map(|_| random_ring(r, &mut seed)).collect();
            let expected = naive_forward(&z, r);
            let mut got = z.clone();
            dif_in_place(&mut got, r, false);
            let bits = t.trailing_zeros();
            for j in 0..t {
                let slot = ((j as u64).reverse_bits() >> (64 - bits)) as usize;
                for c in 0..r {
                    assert!(
                        (got[slot][c] - expected[j][c]).abs() < 1e-9,
                        "r={r} t={t} j={j} slot={slot}"
                    );
                }
            }
            // Natural order must match too.
            assert!(max_diff(&forward(&z, r), &expected) < 1e-9);
        }
    }

    #[test]
    fn forward_then_inverse_is_identity_over_t() {
        for (r, t) in [(2usize, 1usize), (2, 2), (4, 4), (8, 4), (8, 8)] {
            let mut seed = 99u64;
            let z: Vec<RingElem> = (0..t).map(|_| random_ring(r, &mut seed)).collect();
            let f = forward(&z, r);
            let back = inverse(&f, r);
            for j in 0..t {
                for c in 0..r {
                    assert!(
                        (back[j][c] - z[j][c] / t as f64).abs() < 1e-9,
                        "F+F- = I/t, r={r} t={t} j={j}"
                    );
                }
            }
        }
    }
}
