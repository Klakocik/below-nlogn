//! Section 3's cancellation networks, small enough to actually run.
//!
//! The published instance fixes `h = 100`: `C(100,3) = 161700` triples,
//! roughly `1.77e20` wires, and an eight-row schedule arranged so that
//! arbitrary initial scratch values cancel. A forward invocation is then just
//! `y_S += x_S` (XOR in the bit case), and three invocations swap the two
//! banks. That network is a fixed object, which is lovely for a theorem and
//! impossible for a test harness.
//!
//! So we instantiate it for small `h` and check the algebra by execution —
//! same trick, smaller shoes. Lemma 3.2 and its bit-case corollary are
//! verified for `h` up to 10; the published wire counts and rank deficits live
//! in [`crate::asymptotics`].

use crate::fft::Complex;
use std::collections::HashMap;

/// All three-element subsets of `{0, ..., h-1}`, in lexicographic order.
pub fn triples(h: usize) -> Vec<[usize; 3]> {
    let mut out = Vec::new();
    for a in 0..h {
        for b in a + 1..h {
            for c in b + 1..h {
                out.push([a, b, c]);
            }
        }
    }
    out
}

fn intersection_size(a: &[usize; 3], b: &[usize; 3]) -> usize {
    a.iter().filter(|x| b.contains(x)).count()
}

/// One instantiated motif network over `h`-element ground set.
pub struct Motif {
    pub h: usize,
    /// `true` selects the Gaussian-dyadic (`Z[i, 1/2]`) construction,
    /// `false` the bit (`F_2`) construction.
    pub complex: bool,
    pub triples: Vec<[usize; 3]>,
    index: HashMap<[usize; 3], usize>,
    /// Ordered neighbor pairs `(S, T)` carrying auxiliary wires `A_{S,T}`.
    pub pairs: Vec<(usize, usize)>,
    pair_index: HashMap<(usize, usize), usize>,
}

impl Motif {
    pub fn new(h: usize, complex: bool) -> Self {
        assert!(h >= 3);
        let triples = triples(h);
        let index: HashMap<[usize; 3], usize> = triples
            .iter()
            .copied()
            .enumerate()
            .map(|(i, t)| (t, i))
            .collect();
        let mut pairs = Vec::new();
        for s in 0..triples.len() {
            for t in 0..triples.len() {
                if s == t {
                    continue;
                }
                let inter = intersection_size(&triples[s], &triples[t]);
                let neighbor = if complex { inter % 2 == 0 } else { inter == 1 };
                if neighbor {
                    pairs.push((s, t));
                }
            }
        }
        let pair_index = pairs.iter().enumerate().map(|(i, &p)| (p, i)).collect();
        Motif {
            h,
            complex,
            triples,
            index,
            pairs,
            pair_index,
        }
    }

    pub fn triple_count(&self) -> usize {
        self.triples.len()
    }

    /// Index of a triple in the canonical ordering.
    pub fn triple_position(&self, t: &[usize; 3]) -> Option<usize> {
        self.index.get(t).copied()
    }

    /// Index of an ordered neighbor pair `(S, T)` in the side-wire list.
    pub fn pair_position(&self, s: usize, t: usize) -> Option<usize> {
        self.pair_index.get(&(s, t)).copied()
    }

    pub fn scratch_len(&self) -> usize {
        self.pairs.len()
    }

    fn side_coefficient(&self, s: usize, t: usize) -> f64 {
        let inter = intersection_size(&self.triples[s], &self.triples[t]);
        if self.complex {
            // Side injection coefficient -(|S ∩ T| - 1)/2 on neighbors.
            -((inter as f64) - 1.0) / 2.0
        } else {
            1.0
        }
    }

    /// Bit construction: one forward invocation `X -> Y` over `F_2`.
    /// `y_S` becomes `y_S XOR x_S`; all scratch values are restored.
    pub fn invoke_bit(&self, x: &[u64], y: &mut [u64], a: &mut [u64], c: &mut [u64]) {
        let v = self.triple_count();
        assert_eq!(x.len(), v);
        assert_eq!(y.len(), v);
        assert_eq!(a.len(), self.pairs.len());
        assert_eq!(c.len(), self.h);

        let mut acc = vec![0u64; v];

        // Row 0: subtract side injection.
        acc.fill(0);
        for (i, &(s, _t)) in self.pairs.iter().enumerate() {
            acc[s] ^= a[i];
        }
        for s in 0..v {
            y[s] ^= acc[s];
        }

        // Row 1: subtract scatter.
        for s in 0..v {
            let mut r = 0u64;
            for &i in &self.triples[s] {
                r ^= c[i];
            }
            y[s] ^= r;
        }

        // Row 2: copy x into the side wires.
        for (i, &(_s, t)) in self.pairs.iter().enumerate() {
            a[i] ^= x[t];
        }

        // Row 3: gather x into the centers.
        for i in 0..self.h {
            let mut r = 0u64;
            for (t, tri) in self.triples.iter().enumerate() {
                if tri.contains(&i) {
                    r ^= x[t];
                }
            }
            c[i] ^= r;
        }

        // Row 4: scatter with the updated centers.
        for s in 0..v {
            let mut r = 0u64;
            for &i in &self.triples[s] {
                r ^= c[i];
            }
            y[s] ^= r;
        }

        // Row 5: side injection with the updated side wires.
        acc.iter_mut().for_each(|v| *v = 0);
        for (i, &(s, _t)) in self.pairs.iter().enumerate() {
            acc[s] ^= a[i];
        }
        for s in 0..v {
            y[s] ^= acc[s];
        }

        // Row 6: undo gather.
        for i in 0..self.h {
            let mut r = 0u64;
            for (t, tri) in self.triples.iter().enumerate() {
                if tri.contains(&i) {
                    r ^= x[t];
                }
            }
            c[i] ^= r;
        }

        // Row 7: undo copy.
        for (i, &(_s, t)) in self.pairs.iter().enumerate() {
            a[i] ^= x[t];
        }
    }

    /// Gaussian-dyadic construction: one forward invocation `X -> Y`.
    /// `y_S` becomes `y_S + x_S` (complex addition); scratch is restored.
    pub fn invoke_complex(
        &self,
        x: &[Complex],
        y: &mut [Complex],
        a: &mut [Complex],
        c: &mut [Complex],
        c_star: &mut Complex,
    ) {
        self.invoke_complex_signed(x, y, a, c, c_star, 1.0);
    }

    /// Signed Gaussian-dyadic invocation: `y_S` becomes `y_S + sign * x_S`.
    /// `sign = -1` is the "inverse forward schedule from Y to X" of Lemma 3.2.
    pub fn invoke_complex_signed(
        &self,
        x: &[Complex],
        y: &mut [Complex],
        a: &mut [Complex],
        c: &mut [Complex],
        c_star: &mut Complex,
        sign: f64,
    ) {
        let v = self.triple_count();
        assert_eq!(x.len(), v);
        assert_eq!(y.len(), v);
        assert_eq!(a.len(), self.pairs.len());
        assert_eq!(c.len(), self.h);

        let side = |me: &Self, a: &[Complex], s: usize| -> Complex {
            let mut acc = Complex::default();
            for (i, &(ps, _pt)) in me.pairs.iter().enumerate() {
                if ps == s {
                    let coef = me.side_coefficient(ps, _pt);
                    acc += a[i] * coef;
                }
            }
            acc
        };
        let scatter = |c: &[Complex], c_star: Complex, s: usize| -> Complex {
            let mut acc = Complex::default();
            for &i in &self.triples[s] {
                acc += c[i];
            }
            (acc - c_star) / 2.0
        };

        // Row 0: subtract side injection.
        for s in 0..v {
            y[s] -= side(self, a, s) * sign;
        }
        // Row 1: subtract scatter.
        for s in 0..v {
            y[s] -= scatter(c, *c_star, s) * sign;
        }
        // Row 2: copy.
        for (i, &(_s, t)) in self.pairs.iter().enumerate() {
            a[i] += x[t];
        }
        // Row 3: gather.
        for i in 0..self.h {
            let mut acc = Complex::default();
            for (t, tri) in self.triples.iter().enumerate() {
                if tri.contains(&i) {
                    acc += x[t];
                }
            }
            c[i] += acc;
        }
        {
            let mut acc = Complex::default();
            for &v in x {
                acc += v;
            }
            *c_star += acc;
        }
        // Row 4: scatter.
        for s in 0..v {
            y[s] += scatter(c, *c_star, s) * sign;
        }
        // Row 5: side injection.
        for s in 0..v {
            y[s] += side(self, a, s) * sign;
        }
        // Row 6: undo gather.
        for i in 0..self.h {
            let mut acc = Complex::default();
            for (t, tri) in self.triples.iter().enumerate() {
                if tri.contains(&i) {
                    acc += x[t];
                }
            }
            c[i] -= acc;
        }
        {
            let mut acc = Complex::default();
            for &v in x {
                acc += v;
            }
            *c_star -= acc;
        }
        // Row 7: undo copy.
        for (i, &(_s, t)) in self.pairs.iter().enumerate() {
            a[i] -= x[t];
        }
    }

    /// Bit construction, Lemma 3.2: three invocations exchange the banks.
    pub fn exchange_bits(&self, x: &mut [u64], y: &mut [u64]) {
        let mut a = vec![0u64; self.pairs.len()];
        let mut c = vec![0u64; self.h];
        self.invoke_bit(x, y, &mut a, &mut c); // (x, y+x)
        self.invoke_bit(y, x, &mut a, &mut c); // (x-(y+x), y+x) = (-y, y+x)
        self.invoke_bit(x, y, &mut a, &mut c); // (-y, y+x-y) = (-y, x)
                                               // Over F_2 negation is the identity: (y, x).
    }

    /// Gaussian-dyadic construction, Lemma 3.2: three invocations send
    /// `(X, Y)` to `(-Y, X)`: forward `X -> Y`, inverse forward `Y -> X`
    /// (i.e. `X <- X - Y`), then forward `X -> Y`.
    pub fn exchange_complex(&self, x: &mut [Complex], y: &mut [Complex]) {
        let mut a = vec![Complex::default(); self.pairs.len()];
        let mut c = vec![Complex::default(); self.h];
        let mut c_star = Complex::default();
        self.invoke_complex_signed(x, y, &mut a, &mut c, &mut c_star, 1.0); // (x, y+x)
        self.invoke_complex_signed(y, x, &mut a, &mut c, &mut c_star, -1.0); // (-y, y+x)
        self.invoke_complex_signed(x, y, &mut a, &mut c, &mut c_star, 1.0); // (-y, x)
    }

    /// Paper's neighbor counts per fixed triple: `z_b` and `z_c`.
    /// By symmetry (`S_3`-invariance of the intersection rule) this is the
    /// same for every triple; we compute it for triple 0.
    pub fn neighbors_per_triple(&self) -> usize {
        let v = self.triple_count();
        let t = 0usize;
        (0..v)
            .filter(|&s| {
                s != t && {
                    let inter = intersection_size(&self.triples[s], &self.triples[t]);
                    if self.complex {
                        inter % 2 == 0
                    } else {
                        inter == 1
                    }
                }
            })
            .count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state >> 11
    }

    #[test]
    fn neighbor_counts_match_paper_formulas() {
        // z_b = 3 * C(h-3, 2); z_c = C(h-3, 3) + 3(h-3).
        for h in 3..=10usize {
            let bit = Motif::new(h, false);
            let expected_b = 3 * (h.saturating_sub(3)) * (h.saturating_sub(4)) / 2;
            let per_b = bit
                .triples
                .iter()
                .map(|t| {
                    bit.triples
                        .iter()
                        .filter(|s| *s != t && intersection_size(s, t) == 1)
                        .count()
                })
                .collect::<Vec<_>>();
            assert!(per_b.iter().all(|&x| x == expected_b), "h={h}");

            let cx = Motif::new(h, true);
            let m = h - 3;
            let expected_c = m * m.saturating_sub(1) * m.saturating_sub(2) / 6 + 3 * m;
            let per_c = cx
                .triples
                .iter()
                .map(|t| {
                    cx.triples
                        .iter()
                        .filter(|s| *s != t && intersection_size(s, t) % 2 == 0)
                        .count()
                })
                .collect::<Vec<_>>();
            assert!(per_c.iter().all(|&x| x == expected_c), "h={h}");
        }
    }

    #[test]
    fn bit_invocation_adds_and_restores_scratch() {
        let motif = Motif::new(6, false);
        let v = motif.triple_count();
        let mut state = 7u64;
        for _ in 0..4 {
            let x: Vec<u64> = (0..v).map(|_| rng(&mut state) & 1).collect();
            let y0: Vec<u64> = (0..v).map(|_| rng(&mut state) & 1).collect();
            let a0: Vec<u64> = (0..motif.pairs.len())
                .map(|_| rng(&mut state) & 1)
                .collect();
            let c0: Vec<u64> = (0..motif.h).map(|_| rng(&mut state) & 1).collect();
            let mut y = y0.clone();
            let mut a = a0.clone();
            let mut c = c0.clone();
            motif.invoke_bit(&x, &mut y, &mut a, &mut c);
            for s in 0..v {
                assert_eq!(y[s], y0[s] ^ x[s], "target update");
            }
            assert_eq!(a, a0, "side wires restored");
            assert_eq!(c, c0, "center wires restored");
        }
    }

    #[test]
    fn complex_invocation_adds_and_restores_scratch() {
        let motif = Motif::new(5, true);
        let v = motif.triple_count();
        let mut state = 11u64;
        let rand_c = |state: &mut u64| {
            Complex::new(
                (rng(state) % 21) as f64 - 10.0,
                (rng(state) % 21) as f64 - 10.0,
            )
        };
        for _ in 0..4 {
            let x: Vec<Complex> = (0..v).map(|_| rand_c(&mut state)).collect();
            let y0: Vec<Complex> = (0..v).map(|_| rand_c(&mut state)).collect();
            let a0: Vec<Complex> = (0..motif.pairs.len()).map(|_| rand_c(&mut state)).collect();
            let c0: Vec<Complex> = (0..motif.h).map(|_| rand_c(&mut state)).collect();
            let cs0 = rand_c(&mut state);
            let mut y = y0.clone();
            let mut a = a0.clone();
            let mut c = c0.clone();
            let mut cs = cs0;
            motif.invoke_complex(&x, &mut y, &mut a, &mut c, &mut cs);
            for s in 0..v {
                assert!((y[s] - (y0[s] + x[s])).abs() < 1e-9, "target update at {s}");
            }
            for (i, (g, o)) in a.iter().zip(&a0).enumerate() {
                assert!((*g - *o).abs() < 1e-9, "side wire {i} restored");
            }
            for (i, (g, o)) in c.iter().zip(&c0).enumerate() {
                assert!((*g - *o).abs() < 1e-9, "center wire {i} restored");
            }
            assert!((cs - cs0).abs() < 1e-9, "central star restored");
        }
    }

    #[test]
    fn three_invocations_exchange_bit_banks() {
        let motif = Motif::new(5, false);
        let v = motif.triple_count();
        let mut state = 3u64;
        let x0: Vec<u64> = (0..v).map(|_| rng(&mut state) & 1).collect();
        let y0: Vec<u64> = (0..v).map(|_| rng(&mut state) & 1).collect();
        let mut x = x0.clone();
        let mut y = y0.clone();
        motif.exchange_bits(&mut x, &mut y);
        assert_eq!(x, y0);
        assert_eq!(y, x0);
    }

    #[test]
    fn three_invocations_negate_exchange_complex_banks() {
        let motif = Motif::new(4, true);
        let v = motif.triple_count();
        let mut state = 5u64;
        let x: Vec<Complex> = (0..v)
            .map(|_| Complex::new((rng(&mut state) % 11) as f64, 0.0))
            .collect();
        let y: Vec<Complex> = (0..v)
            .map(|_| Complex::new(0.0, (rng(&mut state) % 11) as f64))
            .collect();
        let (x0, y0) = (x.clone(), y.clone());
        let mut x = x;
        let mut y = y;
        motif.exchange_complex(&mut x, &mut y);
        for i in 0..v {
            assert!((x[i] + y0[i]).abs() < 1e-9, "X becomes -Y");
            assert!((y[i] - x0[i]).abs() < 1e-9, "Y becomes X");
        }
    }
}
