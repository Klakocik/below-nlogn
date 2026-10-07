//! XOR swaps, three-shear block interchange, and the exact rational
//! factorization behind Section 4's fast address rearrangement.
//!
//! The move the paper pulls: give every wire an address permutation, keep the
//! edge operators cheap, and count only ranks. The finite algebra is
//!
//! * Lemma 4.3 — every rational matrix `A` of rank `a` factors as
//!   `A = E1 · Pi · E2` with `E1, E2` invertible lower triangular and `Pi` a
//!   partial permutation with exactly `a` ones;
//! * Lemma 4.4 — the shear `(H, D) -> (H + A·D, D)` then costs exactly `a`
//!   field interchanges plus linear work;
//! * the chunk swap itself is three shears wearing a trench coat
//!   (`D -= H`, `H += D`, `D = H - D`).
//!
//! All three are implemented in exact arithmetic and tested. The *sublinear*
//! tape schedule of Proposition 4.10 is not reproduced; see
//! `docs/paper-map.md`.

use std::collections::HashMap;

/// In-place XOR swap of two equal-length word blocks (the classic trick the
/// paper's bit construction rests on).
#[allow(clippy::manual_swap)] // the XOR dance is the subject here, not an accident
pub fn xor_swap(a: &mut [u64], b: &mut [u64]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter_mut().zip(b.iter_mut()) {
        *x ^= *y;
        *y ^= *x;
        *x ^= *y;
    }
}

/// Interchanges the two contiguous digit vectors by the three-shear identity
/// `D -= H; H += D; D = H - D`, all arithmetic modulo `q`.
pub fn shear_interchange(h: &mut [i64], d: &mut [i64], q: i64) {
    assert_eq!(h.len(), d.len());
    let m = |v: i64| v.rem_euclid(q);
    for i in 0..h.len() {
        d[i] = m(d[i] - h[i]);
        h[i] = m(h[i] + d[i]);
        d[i] = m(h[i] - d[i]);
    }
}

/// Exact rational number with `i128` components, for the factorization demo.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Frac {
    pub num: i128,
    pub den: i128,
}

impl Frac {
    pub fn new(num: i128, den: i128) -> Self {
        assert!(den != 0);
        let mut n = num;
        let mut d = den;
        if d < 0 {
            n = -n;
            d = -d;
        }
        let g = gcd(n.unsigned_abs(), d as u128) as i128;
        Frac {
            num: n / g,
            den: d / g,
        }
    }

    pub fn int(v: i128) -> Self {
        Frac { num: v, den: 1 }
    }

    pub fn zero() -> Self {
        Frac { num: 0, den: 1 }
    }

    pub fn one() -> Self {
        Frac { num: 1, den: 1 }
    }

    pub fn is_zero(self) -> bool {
        self.num == 0
    }
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

impl std::ops::Add for Frac {
    type Output = Frac;
    fn add(self, o: Frac) -> Frac {
        Frac::new(self.num * o.den + o.num * self.den, self.den * o.den)
    }
}

impl std::ops::Sub for Frac {
    type Output = Frac;
    fn sub(self, o: Frac) -> Frac {
        Frac::new(self.num * o.den - o.num * self.den, self.den * o.den)
    }
}

impl std::ops::Mul for Frac {
    type Output = Frac;
    fn mul(self, o: Frac) -> Frac {
        Frac::new(self.num * o.num, self.den * o.den)
    }
}

impl std::ops::Div for Frac {
    type Output = Frac;
    fn div(self, o: Frac) -> Frac {
        assert!(!o.is_zero(), "rational division by zero");
        Frac::new(self.num * o.den, self.den * o.num)
    }
}

impl std::ops::Neg for Frac {
    type Output = Frac;
    fn neg(self) -> Frac {
        Frac {
            num: -self.num,
            den: self.den,
        }
    }
}

/// A dense rational matrix.
pub type Mat = Vec<Vec<Frac>>;

pub fn mat_zeros(m: usize) -> Mat {
    vec![vec![Frac::zero(); m]; m]
}

pub fn mat_identity(m: usize) -> Mat {
    let mut a = mat_zeros(m);
    for (i, row) in a.iter_mut().enumerate() {
        row[i] = Frac::one();
    }
    a
}

pub fn mat_mul(a: &Mat, b: &Mat) -> Mat {
    let m = a.len();
    let mut out = mat_zeros(m);
    for i in 0..m {
        for k in 0..m {
            if a[i][k].is_zero() {
                continue;
            }
            for j in 0..m {
                out[i][j] = out[i][j] + a[i][k] * b[k][j];
            }
        }
    }
    out
}

fn mat_inverse(a: &Mat) -> Mat {
    let n = a.len();
    let mut work = a.clone();
    let mut inv = mat_identity(n);
    for col in 0..n {
        let pivot = (col..n)
            .find(|&r| !work[r][col].is_zero())
            .unwrap_or_else(|| panic!("singular matrix at column {col}"));
        work.swap(col, pivot);
        inv.swap(col, pivot);
        let scale = work[col][col];
        for j in 0..n {
            work[col][j] = work[col][j] / scale;
            inv[col][j] = inv[col][j] / scale;
        }
        for r in 0..n {
            if r == col || work[r][col].is_zero() {
                continue;
            }
            let factor = work[r][col];
            for j in 0..n {
                work[r][j] = work[r][j] - factor * work[col][j];
                inv[r][j] = inv[r][j] - factor * inv[col][j];
            }
        }
    }
    inv
}

fn is_lower_triangular(a: &Mat) -> bool {
    for (i, row) in a.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            if j > i && !v.is_zero() {
                return false;
            }
        }
    }
    true
}

/// Lemma 4.3: exact lower-lower factorization `A = E1 · Pi · E2`.
///
/// Returns `(E1, Pi, E2, rank)` where `rank` is the number of ones in `Pi`.
pub fn lower_lower_factorization(a: &Mat) -> (Mat, Mat, Mat, usize) {
    let m = a.len();
    assert!(a.iter().all(|r| r.len() == m), "matrix must be square");
    let mut work = a.clone();
    let mut l = mat_identity(m);
    let mut r = mat_identity(m);
    let mut pi = mat_zeros(m);
    let mut active_rows: Vec<usize> = (0..m).collect();
    let mut active_cols: Vec<usize> = (0..m).collect();
    while !active_rows.is_empty() && !active_cols.is_empty() {
        // Topmost active row with a nonzero active entry; rightmost such entry.
        let mut chosen: Option<(usize, usize)> = None;
        for &ri in &active_rows {
            if let Some(&ci) = active_cols
                .iter()
                .rev()
                .find(|&&ci| !work[ri][ci].is_zero())
            {
                chosen = Some((ri, ci));
                break;
            }
        }
        let Some((pr, pc)) = chosen else { break };
        // Clear entries left of the pivot in row `pr` by column operations.
        let left_cols: Vec<usize> = active_cols.iter().copied().filter(|&c| c < pc).collect();
        for ci in left_cols {
            if work[pr][ci].is_zero() {
                continue;
            }
            let factor = work[pr][ci] / work[pr][pc];
            for i in 0..m {
                work[i][ci] = work[i][ci] - factor * work[i][pc];
                r[i][ci] = r[i][ci] - factor * r[i][pc];
            }
        }
        // Clear entries below the pivot by row operations.
        let below_rows: Vec<usize> = active_rows.iter().copied().filter(|&x| x > pr).collect();
        for ri in below_rows {
            if work[ri][pc].is_zero() {
                continue;
            }
            let factor = work[ri][pc] / work[pr][pc];
            for j in 0..m {
                work[ri][j] = work[ri][j] - factor * work[pr][j];
                l[ri][j] = l[ri][j] - factor * l[pr][j];
            }
        }
        // Scale the pivot column to one.
        let scale = work[pr][pc];
        for i in 0..m {
            work[i][pc] = work[i][pc] / scale;
            r[i][pc] = r[i][pc] / scale;
        }
        debug_assert!(work[pr][pc] == Frac::one());
        pi[pr][pc] = Frac::one();
        active_rows.retain(|&x| x != pr);
        active_cols.retain(|&x| x != pc);
    }
    let rank = pi.iter().flatten().filter(|v| !v.is_zero()).count();
    if cfg!(debug_assertions) {
        let check = mat_mul(&mat_mul(&l, a), &r);
        for i in 0..m {
            for j in 0..m {
                assert_eq!(check[i][j], pi[i][j], "L·A·R != Pi at ({i},{j})");
            }
        }
    }
    let e1 = mat_inverse(&l);
    let e2 = mat_inverse(&r);
    debug_assert!(is_lower_triangular(&e1) && is_lower_triangular(&e2));
    (e1, pi, e2, rank)
}

/// Reduces a rational to `Z/qZ` (requires `q` prime and denominators coprime to `q`).
fn frac_mod(v: Frac, q: u64) -> u64 {
    let den = v.den.rem_euclid(q as i128) as u64;
    let num = v.num.rem_euclid(q as i128) as u64;
    mul_mod(num, pow_mod(den, q - 2, q), q)
}

fn mul_mod(a: u64, b: u64, q: u64) -> u64 {
    ((a as u128 * b as u128) % q as u128) as u64
}

fn pow_mod(mut b: u64, mut e: u64, q: u64) -> u64 {
    let mut acc = 1u64 % q;
    b %= q;
    while e > 0 {
        if e & 1 == 1 {
            acc = mul_mod(acc, b, q);
        }
        b = mul_mod(b, b, q);
        e >>= 1;
    }
    acc
}

fn mat_vec_mul_mod(a: &Mat, v: &[u64], q: u64) -> Vec<u64> {
    let m = a.len();
    let mut out = vec![0u64; m];
    for i in 0..m {
        let mut acc = 0u128;
        for j in 0..m {
            acc += (frac_mod(a[i][j], q) as u128) * (v[j] as u128);
            acc %= q as u128;
        }
        out[i] = acc as u64;
    }
    out
}

/// Lemma 4.4: performs `(H, D) -> (H + A·D, D)` over `Z/qZ` using exactly
/// `rank(A)` field interchanges, via the factorization `A = E1·Pi·E2`.
/// Returns the number of interchanges used.
pub fn shear_via_factorization(h: &mut [u64], d: &mut [u64], a: &Mat, q: u64) -> usize {
    let m = h.len();
    assert_eq!(d.len(), m);
    assert_eq!(a.len(), m);
    let (e1, pi, e2, _rank) = lower_lower_factorization(a);
    let e1_inv = mat_inverse(&e1);
    let e2_inv = mat_inverse(&e2);
    let add = |x: u64, y: u64| {
        let s = x + y;
        if s >= q {
            s - q
        } else {
            s
        }
    };
    // D <- E2 D and H <- E1^{-1} H
    let d_new = mat_vec_mul_mod(&e2, d, q);
    d.copy_from_slice(&d_new);
    let h_new = mat_vec_mul_mod(&e1_inv, h, q);
    h.copy_from_slice(&h_new);
    // H <- H + Pi D, one interchange per one of Pi.
    let mut used = 0;
    for i in 0..m {
        for j in 0..m {
            if !pi[i][j].is_zero() {
                // The paper's single interchange: (h, d) -> (h + d, d).
                h[i] = add(h[i], d[j]);
                used += 1;
            }
        }
    }
    // H <- E1 H and D <- E2^{-1} D
    let h_new = mat_vec_mul_mod(&e1, h, q);
    h.copy_from_slice(&h_new);
    let d_new = mat_vec_mul_mod(&e2_inv, d, q);
    d.copy_from_slice(&d_new);
    used
}

/// A bit-reversal permutation as an index map (used by the transform layouts).
pub fn bit_reversal_map(n: usize) -> Vec<usize> {
    if n <= 1 {
        return vec![0usize; n];
    }
    let bits = n.trailing_zeros();
    let mut map = vec![0usize; n];
    for (i, m) in map.iter_mut().enumerate() {
        *m = ((i as u64).reverse_bits() >> (64 - bits)) as usize;
    }
    map
}

/// Splits `[P] x [R] x [B]` rows cyclically into `w` role streams
/// (Lemma 2.2, item 4). Returns, for each role, the row classes it contains.
pub fn cyclic_row_classes(
    p: usize,
    r: usize,
    b: usize,
    w: usize,
) -> Vec<Vec<(usize, usize, usize)>> {
    assert!(w >= 1);
    let mut roles: Vec<Vec<(usize, usize, usize)>> = vec![Vec::new(); w];
    for pi in 0..p {
        for ri in 0..r {
            for bi in 0..b {
                roles[ri % w].push((pi, ri, bi));
            }
        }
    }
    roles
}

/// Simple occurrence index of `(source, target)` pairs, used by the motif demo.
pub fn pair_index(pairs: &[(usize, usize)]) -> HashMap<(usize, usize), usize> {
    pairs.iter().enumerate().map(|(i, &p)| (p, i)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(rows: &[&[i128]]) -> Mat {
        rows.iter()
            .map(|r| r.iter().map(|&v| Frac::int(v)).collect())
            .collect()
    }

    #[test]
    fn xor_swap_works() {
        let mut a = [0x1234_5678_9ABC_DEF0u64, 42];
        let mut b = [0x0FED_CBA9_8765_4321u64, 7];
        let (oa, ob) = (a, b);
        xor_swap(&mut a, &mut b);
        assert_eq!(a, ob);
        assert_eq!(b, oa);
    }

    #[test]
    fn three_shear_swaps_digits() {
        let mut rng = 123u64;
        let next = |rng: &mut u64| {
            *rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            (*rng >> 33) as i64 % 1000
        };
        let mut h: Vec<i64> = (0..32).map(|_| next(&mut rng)).collect();
        let mut d: Vec<i64> = (0..32).map(|_| next(&mut rng)).collect();
        let (oh, od) = (h.clone(), d.clone());
        shear_interchange(&mut h, &mut d, 1_000_003);
        assert_eq!(h, od);
        assert_eq!(d, oh);
    }

    #[test]
    fn factorization_reconstructs_and_counts_rank() {
        let cases: Vec<Mat> = vec![
            m(&[&[1, 2], &[3, 4]]),
            m(&[&[1, 2, 3], &[4, 5, 6], &[7, 8, 9]]), // rank 2
            m(&[&[0, 0, 0], &[0, 0, 0], &[0, 0, 0]]), // rank 0
            m(&[&[2, 0, 0, 0], &[0, 3, 0, 0], &[0, 0, 0, 5], &[0, 0, 0, 0]]),
            m(&[&[1, 0, -2, 1], &[3, 1, 0, 2], &[0, 0, 0, 0], &[-1, 2, 1, 0]]),
        ];
        for a in &cases {
            let (e1, pi, e2, rank) = lower_lower_factorization(a);
            // A = E1 · Pi · E2
            assert_eq!(mat_mul(&mat_mul(&e1, &pi), &e2), *a);
            assert!(is_lower_triangular(&e1) && is_lower_triangular(&e2));
            assert_eq!(rank, pi.iter().flatten().filter(|v| !v.is_zero()).count());
        }
    }

    #[test]
    fn shear_matches_direct_matrix_vector_product() {
        let q = 1_000_003u64;
        let a = m(&[&[1, 2, -1], &[0, 3, 4], &[2, 0, 1]]);
        let mut h = [11u64, 22, 33];
        let mut d = [5u64, 6, 7];
        let expected = mat_vec_mul_mod(&a, &d, q);
        let used = shear_via_factorization(&mut h, &mut d, &a, q);
        for i in 0..3 {
            assert_eq!(h[i], (expected[i] + [11u64, 22, 33][i]) % q, "row {i}");
        }
        assert_eq!(used, 3); // full rank
        assert_eq!(d, [5u64, 6, 7]);
    }

    #[test]
    fn bit_reversal_map_involutive() {
        assert_eq!(bit_reversal_map(1), vec![0]);
        let map = bit_reversal_map(16);
        let mut seen = [false; 16];
        for (i, &j) in map.iter().enumerate() {
            assert!(!seen[j]);
            seen[j] = true;
            assert_eq!(map[j], i);
        }
    }
}
