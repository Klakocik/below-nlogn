//! Custom benchmark harness (`cargo bench`), no external dependencies.
//!
//! Measures schoolbook / Karatsuba / NTT multiplication across bit sizes and
//! prints the honest picture: the exact NTT tracks `n log n`-like scaling,
//! and no "below n log n" effect is observable at any realizable size.

use below_nlogn::limbs::BigUint;
use below_nlogn::mul::{self, Algorithm};
use std::time::{Duration, Instant};

fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

fn random_bits(state: &mut u64, bits: usize) -> BigUint {
    let mut v = BigUint::zero();
    for i in 0..bits {
        if splitmix(state) & 1 == 1 {
            v.set_bit(i);
        }
    }
    v
}

fn bench_algo(a: &BigUint, b: &BigUint, alg: Algorithm, repeats: usize) -> Duration {
    let mut best = Duration::MAX;
    for _ in 0..repeats {
        let start = Instant::now();
        let r = mul::mul(a, b, alg);
        let elapsed = start.elapsed();
        std::hint::black_box(r);
        best = best.min(elapsed);
    }
    best
}

fn main() {
    let max_bits: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(1 << 20);
    println!("# below-nlogn multiplication benchmark");
    println!("# algorithms: schoolbook (S), karatsuba (K), ntt (N)");
    println!("# sizes in bits; times in milliseconds (best of repeats)");
    println!(
        "{:>10} {:>14} {:>14} {:>14} {:>16}",
        "bits", "S ms", "K ms", "N ms", "N ns/bit"
    );
    let mut state = 0x5EED_1234_5678_9ABCu64;
    let mut bits = 4096usize;
    while bits <= max_bits {
        let a = random_bits(&mut state, bits);
        let b = random_bits(&mut state, bits);
        let repeats = if bits <= 1 << 16 {
            5
        } else if bits <= 1 << 18 {
            3
        } else {
            1
        };
        let s = if bits <= 1 << 18 {
            Some(bench_algo(&a, &b, Algorithm::Schoolbook, repeats))
        } else {
            None
        };
        let k = if bits <= 1 << 20 {
            Some(bench_algo(&a, &b, Algorithm::Karatsuba, repeats))
        } else {
            None
        };
        let n = bench_algo(&a, &b, Algorithm::Ntt, repeats);
        let fmt = |d: Option<Duration>| match d {
            Some(t) => format!("{:.3}", t.as_secs_f64() * 1e3),
            None => "-".to_string(),
        };
        println!(
            "{bits:>10} {:>14} {:>14} {:>14} {:>16.2}",
            fmt(s),
            fmt(k),
            fmt(Some(n)),
            n.as_secs_f64() * 1e9 / bits as f64
        );
        bits *= 4;
    }
    println!();
    println!("# Reminder: the paper's O(n (lg n)^(1-κ)) with κ = 2^-182 is a");
    println!("# bit-complexity statement. The ratio to n lg n is (lg n)^(-κ);");
    println!("# see `below-nlogn info` for why this is unobservable in practice.");
}
