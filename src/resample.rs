//! Gaussian resampling, permutations left inside the formula (Section 7).
//!
//! The paper connects coprime transform lengths to power-of-two shapes through
//! Harvey–van der Hoeven's Gaussian resampling identity, keeping the
//! coordinate permutations in the formula instead of sorting twice:
//!
//! ```text
//! P_s F_s = 2^{2 alpha^2} B_0 P_t F_t A,
//! ```
//!
//! with `A = S/2`, `B_0 = D' J' C`, `N = C T D`, `J' = N^{-1}/2`, and
//! `||N - I|| < 2.01 e^{-pi alpha^2 theta/2}`, `theta = t/s - 1`.
//!
//! This module builds `S`, `T`, `C`, `D` and `N` as explicit `f64` matrices
//! for small coprime `s < t`, checks the contraction estimate, and verifies
//! the identity numerically. It is a readable demonstration, not a
//! reimplementation of the cited fixed-point contraction algorithms — do not
//! cite it for error bounds.

use crate::fft::Complex;
use std::f64::consts::PI;

/// Dense complex matrix (row-major).
pub type CxMat = Vec<Vec<Complex>>;

fn zeros(rows: usize, cols: usize) -> CxMat {
    vec![vec![Complex::default(); cols]; rows]
}

fn identity(n: usize) -> CxMat {
    let mut m = zeros(n, n);
    for i in 0..n {
        m[i][i] = Complex::new(1.0, 0.0);
    }
    m
}

fn mat_mul(a: &CxMat, b: &CxMat) -> CxMat {
    let n = a.len();
    let k = b.len();
    let m = b[0].len();
    let mut out = zeros(n, m);
    for i in 0..n {
        for t in 0..k {
            let av = a[i][t];
            if av == Complex::default() {
                continue;
            }
            for j in 0..m {
                out[i][j] += av * b[t][j];
            }
        }
    }
    out
}

fn mat_vec(a: &CxMat, v: &[Complex]) -> Vec<Complex> {
    a.iter()
        .map(|row| {
            row.iter()
                .zip(v)
                .fold(Complex::default(), |acc, (x, y)| acc + *x * *y)
        })
        .collect()
}

fn inverse(a: &CxMat) -> CxMat {
    let n = a.len();
    let mut work = a.clone();
    let mut inv = identity(n);
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&r1, &r2| {
                work[r1][col]
                    .abs()
                    .partial_cmp(&work[r2][col].abs())
                    .unwrap()
            })
            .unwrap();
        work.swap(col, pivot);
        inv.swap(col, pivot);
        let scale = work[col][col];
        for j in 0..n {
            work[col][j] = cdiv(work[col][j], scale);
            inv[col][j] = cdiv(inv[col][j], scale);
        }
        for r in 0..n {
            if r == col || work[r][col] == Complex::default() {
                continue;
            }
            let factor = work[r][col];
            for j in 0..n {
                let w = work[col][j];
                let iv = inv[col][j];
                work[r][j] -= factor * w;
                inv[r][j] -= factor * iv;
            }
        }
    }
    inv
}

fn cdiv(a: Complex, b: Complex) -> Complex {
    let d = b.norm_sq();
    Complex::new(
        (a.re * b.re + a.im * b.im) / d,
        (a.im * b.re - a.re * b.im) / d,
    )
}

/// Normalized DFT matrix with the paper's sign: `F_n[k][j] = e^{-2πi jk/n}/n`.
pub fn dft_matrix(n: usize) -> CxMat {
    let mut m = zeros(n, n);
    for (k, row) in m.iter_mut().enumerate() {
        for (j, x) in row.iter_mut().enumerate() {
            *x = Complex::root_of_unity(-1.0, 2 * (j as u64) * (k as u64), n as u64) / n as f64;
        }
    }
    m
}

/// `(P_s u)_j = u_{tj mod s}` as a matrix.
pub fn perm_s_matrix(s: usize, t: usize) -> CxMat {
    assert_eq!(gcd(s, t), 1);
    let mut m = zeros(s, s);
    for j in 0..s {
        m[j][(t * j) % s] = Complex::new(1.0, 0.0);
    }
    m
}

/// `(P_t v)_k = v_{-sk mod t}` as a matrix.
pub fn perm_t_matrix(s: usize, t: usize) -> CxMat {
    assert_eq!(gcd(s, t), 1);
    let mut m = zeros(t, t);
    for k in 0..t {
        m[k][(s * (t - k % t)) % t] = Complex::new(1.0, 0.0);
    }
    m
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn truncation_for(alpha: f64) -> f64 {
    5.0 * alpha + 4.0
}

/// Gaussian map `S: C^s -> C^t` with periodic extension in the input index.
pub fn gaussian_s(s: usize, t: usize, alpha: f64) -> CxMat {
    let mut m = zeros(t, s);
    let cut = truncation_for(alpha);
    let m_range = (cut / s as f64).ceil() as isize + 2;
    for (k, row) in m.iter_mut().enumerate() {
        let center = s as f64 * k as f64 / t as f64;
        for (j, x) in row.iter_mut().enumerate() {
            let mut acc = 0.0;
            for mm in -m_range..=m_range {
                let d = j as f64 + (mm * s as isize) as f64 - center;
                if d.abs() > cut {
                    continue;
                }
                acc += (-PI * (d / alpha).powi(2)).exp();
            }
            *x = Complex::new(acc / alpha, 0.0);
        }
    }
    m
}

/// Gaussian map `T: C^s -> C^t` with periodic extension in the input index.
pub fn gaussian_t(s: usize, t: usize, alpha: f64) -> CxMat {
    let mut m = zeros(t, s);
    let cut = truncation_for(alpha);
    let m_range = (cut / s as f64).ceil() as isize + 2;
    for (k, row) in m.iter_mut().enumerate() {
        for (j, x) in row.iter_mut().enumerate() {
            let mut acc = 0.0;
            for mm in -m_range..=m_range {
                let d = t as f64 * (j as f64 + (mm * s as isize) as f64) / s as f64 - k as f64;
                if d.abs() > cut {
                    continue;
                }
                acc += (-PI * (alpha * d).powi(2)).exp();
            }
            *x = Complex::new(acc, 0.0);
        }
    }
    m
}

/// Selection `C`: coordinate `j` maps to `q_j = floor(t j / s + 1/2)`.
pub fn selector_c(s: usize, t: usize) -> CxMat {
    let mut m = zeros(s, t);
    for (j, row) in m.iter_mut().enumerate() {
        let q = (t as f64 * j as f64 / s as f64 + 0.5).floor() as usize;
        row[q] = Complex::new(1.0, 0.0);
    }
    m
}

/// Diagonal `D_j = e^{pi alpha^2 beta_j^2}` (real, at least one).
pub fn diagonal_d(s: usize, t: usize, alpha: f64) -> CxMat {
    let mut m = zeros(s, s);
    for j in 0..s {
        let q = (t as f64 * j as f64 / s as f64 + 0.5).floor();
        let beta = t as f64 * j as f64 / s as f64 - q;
        m[j][j] = Complex::new((PI * alpha * alpha * beta * beta).exp(), 0.0);
    }
    m
}

/// The resampling operators `(A, B_0, N)` of Lemma 7.2.
pub struct Resampling {
    pub s: usize,
    pub t: usize,
    pub alpha: f64,
    /// `A = S/2` (t x s).
    pub a: CxMat,
    /// `B_0 = D' J' C` (s x t).
    pub b0: CxMat,
    /// `N = C T D` (s x s).
    pub n: CxMat,
}

impl Resampling {
    pub fn new(s: usize, t: usize, alpha: f64) -> Self {
        assert!(
            s >= 2 && t > s && gcd(s, t) == 1,
            "require 2 <= s < t coprime"
        );
        let sm = gaussian_s(s, t, alpha);
        let tm = gaussian_t(s, t, alpha);
        let c = selector_c(s, t);
        let d = diagonal_d(s, t, alpha);
        let n = mat_mul(&mat_mul(&c, &tm), &d);

        let n_inv = inverse(&n);
        let half = Complex::new(0.5, 0.0);
        let jprime: CxMat = n_inv
            .iter()
            .map(|row| row.iter().map(|x| *x * half).collect())
            .collect();
        let scale = 2f64.powf(-(2.0 * alpha * alpha - 2.0));
        let dprime: CxMat = d
            .iter()
            .map(|row| row.iter().map(|x| *x * scale).collect())
            .collect();
        let b0 = mat_mul(&mat_mul(&dprime, &jprime), &c);

        let mut a: CxMat = Vec::with_capacity(t);
        for row in &sm {
            let mut out = Vec::with_capacity(s);
            for x in row {
                out.push(*x * half);
            }
            a.push(out);
        }
        Resampling {
            s,
            t,
            alpha,
            a,
            b0,
            n,
        }
    }

    /// `theta = t/s - 1`.
    pub fn theta(&self) -> f64 {
        self.t as f64 / self.s as f64 - 1.0
    }

    /// The paper's contraction bound `2.01 e^{-pi alpha^2 theta / 2}`.
    pub fn contraction_bound(&self) -> f64 {
        2.01 * (-PI * self.alpha * self.alpha * self.theta() / 2.0).exp()
    }

    /// Max row sum of `|N - I|`, the induced infinity norm of the error.
    pub fn n_minus_i_norm(&self) -> f64 {
        let s = self.s;
        let mut worst: f64 = 0.0;
        for i in 0..s {
            let mut row = 0.0;
            for j in 0..s {
                let v = if i == j {
                    self.n[i][j] - Complex::new(1.0, 0.0)
                } else {
                    self.n[i][j]
                };
                row += v.abs();
            }
            worst = worst.max(row);
        }
        worst
    }

    /// Numerically checks `P_s F_s = 2^{2 alpha^2} B_0 P_t F_t A` on a
    /// pseudorandom complex vector, returning the maximum absolute error.
    pub fn identity_error(&self, seed: u64) -> f64 {
        let s = self.s;
        let t = self.t;
        let mut state = seed | 1;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let x = ((state >> 11) as f64) / (1u64 << 53) as f64;
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let y = ((state >> 11) as f64) / (1u64 << 53) as f64;
            Complex::new(2.0 * x - 1.0, 2.0 * y - 1.0)
        };
        let u: Vec<Complex> = (0..s).map(|_| next()).collect();

        let fs = dft_matrix(s);
        let ft = dft_matrix(t);
        let ps = perm_s_matrix(s, t);
        let pt = perm_t_matrix(s, t);
        let lhs = mat_vec(&mat_mul(&ps, &fs), &u);

        let a_u = mat_vec(&self.a, &u);
        let ft_a_u = mat_vec(&ft, &a_u);
        let pt_ft_a_u = mat_vec(&pt, &ft_a_u);
        let b0_pt_ft_a_u = mat_vec(&self.b0, &pt_ft_a_u);
        let scale = 2f64.powf(2.0 * self.alpha * self.alpha);
        let rhs: Vec<Complex> = b0_pt_ft_a_u.iter().map(|x| *x * scale).collect();

        lhs.iter()
            .zip(&rhs)
            .map(|(x, y)| (*x - *y).abs())
            .fold(0.0, f64::max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn max_mat_diff(a: &CxMat, b: &CxMat) -> f64 {
        let mut m: f64 = 0.0;
        for (ra, rb) in a.iter().zip(b) {
            for (x, y) in ra.iter().zip(rb) {
                m = m.max((*x - *y).abs());
            }
        }
        m
    }

    #[test]
    fn harvey_van_der_hoeven_resampling_identity() {
        // Theorem 4.2 of Harvey–van der Hoeven: T P_s F_s = P_t F_t S,
        // exactly (up to Gaussian-truncation roundoff).
        for (s, t, alpha) in [(5usize, 7usize, 2.0), (11, 13, 1.5), (17, 33, 2.0)] {
            let sm = gaussian_s(s, t, alpha);
            let tm = gaussian_t(s, t, alpha);
            let fs = dft_matrix(s);
            let ft = dft_matrix(t);
            let ps = perm_s_matrix(s, t);
            let pt = perm_t_matrix(s, t);
            let lhs = mat_mul(&tm, &mat_mul(&ps, &fs));
            let rhs = mat_mul(&pt, &mat_mul(&ft, &sm));
            assert!(
                max_mat_diff(&lhs, &rhs) < 1e-9,
                "s={s} t={t} alpha={alpha} diff={}",
                max_mat_diff(&lhs, &rhs)
            );
        }
    }

    #[test]
    fn resampling_identity_composition() {
        for (s, t, alpha) in [(17usize, 33usize, 2.0), (11, 13, 1.5)] {
            let r = Resampling::new(s, t, alpha);
            // P_s F_s = 2^{2a^2} B0 P_t F_t A must hold up to truncation.
            assert!(
                r.identity_error(12345) < 1e-6,
                "s={s} t={t} alpha={alpha} err={}",
                r.identity_error(12345)
            );
        }
    }

    #[test]
    fn identity_matches_theory() {
        for (s, t, alpha) in [(17usize, 33usize, 2.0), (11, 13, 1.5)] {
            let r = Resampling::new(s, t, alpha);
            // The algebraic identity P_s F_s = 2^{2a^2} B0 P_t F_t A must hold
            // up to Gaussian truncation, independent of the contraction size.
            assert!(
                r.identity_error(12345) < 1e-3,
                "s={s} t={t} alpha={alpha} err={}",
                r.identity_error(12345)
            );
        }
    }

    #[test]
    fn contraction_is_small_and_below_bound_order() {
        let r = Resampling::new(17, 33, 2.0);
        let measured = r.n_minus_i_norm();
        assert!(
            measured < r.contraction_bound() * 2.0,
            "measured {measured} vs bound {}",
            r.contraction_bound()
        );
        assert!(measured < 0.5, "N must be invertible by Neumann series");
    }

    #[test]
    fn theta_and_bound_are_sane() {
        let r = Resampling::new(53, 59, 4.0);
        assert!((r.theta() - (59.0 / 53.0 - 1.0)).abs() < 1e-12);
        // alpha^2 < p and theta > p/alpha^4 for a large enough precision p;
        // here we only record the analytic bound shape.
        assert!(r.contraction_bound() > 0.0);
    }
}
