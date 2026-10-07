//! Complex FFT and Bluestein's chirp transform, in plain `f64`.
//!
//! These are the conventional floating-point cousins of the paper's exact
//! normalized transforms over `C[y]/(y^r + 1)` (Sections 6–8). They exist so
//! [`crate::resample`], [`crate::layers`] and [`crate::transforms`] have
//! something readable to compare against; the exact multiplier never touches
//! this code. If you borrow it for real convolution, mind the floating dust.

use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// A complex number in `f64` components.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub fn new(re: f64, im: f64) -> Self {
        Complex { re, im }
    }

    pub fn from_angle(angle: f64) -> Self {
        Complex {
            re: angle.cos(),
            im: angle.sin(),
        }
    }

    pub fn conj(self) -> Self {
        Complex {
            re: self.re,
            im: -self.im,
        }
    }

    pub fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    pub fn abs(self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// `e^{i * sign * pi * numerator / n}` with `numerator` reduced modulo `2n`.
    pub fn root_of_unity(sign: f64, numerator: u64, n: u64) -> Self {
        let reduced = (numerator % (2 * n)) as f64;
        Self::from_angle(sign * PI * reduced / n as f64)
    }
}

impl Add for Complex {
    type Output = Complex;
    fn add(self, rhs: Complex) -> Complex {
        Complex::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, rhs: Complex) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

impl SubAssign for Complex {
    fn sub_assign(&mut self, rhs: Complex) {
        self.re -= rhs.re;
        self.im -= rhs.im;
    }
}

impl Sub for Complex {
    type Output = Complex;
    fn sub(self, rhs: Complex) -> Complex {
        Complex::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl Mul for Complex {
    type Output = Complex;
    fn mul(self, rhs: Complex) -> Complex {
        Complex::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl MulAssign for Complex {
    fn mul_assign(&mut self, rhs: Complex) {
        *self = *self * rhs;
    }
}

impl Mul<f64> for Complex {
    type Output = Complex;
    fn mul(self, rhs: f64) -> Complex {
        Complex::new(self.re * rhs, self.im * rhs)
    }
}

impl Div<f64> for Complex {
    type Output = Complex;
    fn div(self, rhs: f64) -> Complex {
        Complex::new(self.re / rhs, self.im / rhs)
    }
}

impl Neg for Complex {
    type Output = Complex;
    fn neg(self) -> Complex {
        Complex::new(-self.re, -self.im)
    }
}

/// Naive `O(n^2)` DFT in the standard convention: the forward transform is
/// unnormalized, and the inverse transform carries the `1/n` factor. (The
/// paper's *normalized* transforms, where both directions carry `1/n`, are
/// built explicitly in [`crate::transforms`] and [`crate::resample`].)
pub fn dft_naive(values: &[Complex], inverse: bool) -> Vec<Complex> {
    let n = values.len();
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut out = vec![Complex::default(); n];
    for (k, slot) in out.iter_mut().enumerate() {
        let mut acc = Complex::default();
        for (j, &x) in values.iter().enumerate() {
            acc += x * Complex::root_of_unity(sign, 2 * (j as u64) * (k as u64), n as u64);
        }
        *slot = if inverse { acc / n as f64 } else { acc };
    }
    out
}

fn bit_reverse_permute(values: &mut [Complex]) {
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

/// In-place radix-2 FFT with the normalization and sign convention above.
///
/// The input length must be a power of two (use [`bluestein_dft`] otherwise).
pub fn fft(values: &mut [Complex], inverse: bool) {
    let n = values.len();
    assert!(n.is_power_of_two(), "fft length must be a power of two");
    if n <= 1 {
        return;
    }
    bit_reverse_permute(values);
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2usize;
    while len <= n {
        let root = Complex::from_angle(sign * 2.0 * PI / len as f64);
        for block in values.chunks_exact_mut(len) {
            let mut w = Complex::new(1.0, 0.0);
            let (left, right) = block.split_at_mut(len / 2);
            for j in 0..len / 2 {
                let u = left[j];
                let v = right[j] * w;
                left[j] = u + v;
                right[j] = u - v;
                w *= root;
            }
        }
        len <<= 1;
    }
    if inverse {
        for v in values.iter_mut() {
            *v = *v / n as f64;
        }
    }
}

/// Linear convolution of two complex sequences (zero-padded power-of-two FFT).
pub fn convolve(a: &[Complex], b: &[Complex]) -> Vec<Complex> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let out_len = a.len() + b.len() - 1;
    let n = out_len.next_power_of_two();
    let mut fa = vec![Complex::default(); n];
    let mut fb = vec![Complex::default(); n];
    fa[..a.len()].copy_from_slice(a);
    fb[..b.len()].copy_from_slice(b);
    fft(&mut fa, false);
    fft(&mut fb, false);
    for i in 0..n {
        fa[i] *= fb[i];
    }
    fft(&mut fa, true);
    fa.truncate(out_len);
    fa
}

/// Bluestein's chirp transform: an `O(n log n)` DFT for *arbitrary* lengths,
/// used in the paper for the middle transform (Section 8, "Chirp and
/// source-transform procedures").
pub fn bluestein_dft(values: &[Complex], inverse: bool) -> Vec<Complex> {
    let n = values.len();
    if n == 0 {
        return Vec::new();
    }
    if n.is_power_of_two() {
        let mut v = values.to_vec();
        fft(&mut v, inverse);
        return v;
    }
    let n64 = n as u64;
    let sign = if inverse { 1.0 } else { -1.0 };
    // w(j) = e^{sign * i * pi * j^2 / n}; using
    // e^{∓2πijk/n} = w(j) * w(k) * conj(w(k-j)) with the appropriate sign,
    // the DFT becomes a linear convolution `a ⊛ conj(w)` scaled by w(k).
    let w = |j: u64| Complex::root_of_unity(sign, (j * j) % (2 * n64), n64);
    let mut a = vec![Complex::default(); n];
    let mut b = vec![Complex::default(); n];
    for j in 0..n {
        let wj = w(j as u64);
        a[j] = values[j] * wj;
        b[j] = wj.conj();
    }
    let l = (2 * n - 1).next_power_of_two();
    let mut fa = vec![Complex::default(); l];
    let mut fb = vec![Complex::default(); l];
    fa[..n].copy_from_slice(&a);
    fb[..n].copy_from_slice(&b);
    for m in 1..n {
        fb[l - m] = b[m];
    }
    fft(&mut fa, false);
    fft(&mut fb, false);
    for i in 0..l {
        fa[i] *= fb[i];
    }
    fft(&mut fa, true);
    let mut out = vec![Complex::default(); n];
    for k in 0..n {
        let wk = w(k as u64);
        out[k] = if inverse {
            fa[k] * wk / n as f64
        } else {
            fa[k] * wk
        };
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(n: usize) -> Vec<Complex> {
        (0..n)
            .map(|i| Complex::new((i as f64 * 0.7).sin(), (i as f64 * 1.3).cos()))
            .collect()
    }

    fn max_error(a: &[Complex], b: &[Complex]) -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (*x - *y).abs())
            .fold(0.0, f64::max)
    }

    #[test]
    fn fft_matches_naive() {
        for n in [2usize, 4, 8, 16, 32] {
            let x = sample(n);
            let mut f = x.clone();
            fft(&mut f, false);
            assert!(max_error(&f, &dft_naive(&x, false)) < 1e-9, "n={n}");
            let mut inv = f.clone();
            fft(&mut inv, true);
            assert!(max_error(&inv, &x) < 1e-9, "inverse n={n}");
        }
    }

    #[test]
    fn bluestein_matches_naive() {
        for n in [3usize, 5, 7, 11, 13, 17] {
            let x = sample(n);
            assert!(
                max_error(&bluestein_dft(&x, false), &dft_naive(&x, false)) < 1e-9,
                "forward n={n}"
            );
            assert!(
                max_error(&bluestein_dft(&x, true), &dft_naive(&x, true)) < 1e-9,
                "inverse n={n}"
            );
        }
    }

    #[test]
    fn convolution_matches_naive() {
        let a = sample(9);
        let b = sample(7);
        let expected: Vec<Complex> = (0..a.len() + b.len() - 1)
            .map(|k| {
                let mut acc = Complex::default();
                for i in 0..a.len() {
                    if k >= i && k - i < b.len() {
                        acc += a[i] * b[k - i];
                    }
                }
                acc
            })
            .collect();
        assert!(max_error(&convolve(&a, &b), &expected) < 1e-9);
    }
}
