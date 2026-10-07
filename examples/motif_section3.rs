//! Executes the Section 3 scalar motif network on a small ground set `h`
//! and demonstrates its defining property: a forward invocation updates the
//! target bank while restoring every auxiliary wire, and three invocations
//! exchange the two banks (up to sign in the Gaussian-dyadic construction).
//!
//! The paper fixes `h = 100` (giving `C(100,3) = 161700` triples and roughly
//! `1.77e20` wires); enumerating that network is impossible here, but the
//! algebraic cancellation property is independent of `h` and is verified
//! exactly for the small instantiations below.
//!
//! Run with: `cargo run --release --example motif_section3`

use below_nlogn::fft::Complex;
use below_nlogn::motif::Motif;

fn main() {
    println!("Section 3 motif network (Lemma 3.2) — executable small-h checks");
    println!();
    for h in [4usize, 6, 8] {
        let motif = Motif::new(h, false);
        println!(
            "h = {h}: v = C(h,3) = {}, side wires = {}, z_b = {}",
            motif.triple_count(),
            motif.pairs.len(),
            motif.neighbors_per_triple()
        );
    }
    println!();

    let h = 8;
    let bit = Motif::new(h, false);
    let v = bit.triple_count();
    // "Random" deterministic bits.
    let mut state = 0x1234_5678_9ABC_DEF0u64;
    let mut rng = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state >> 11
    };
    // Lemma 3.2, bit case: three invocations exchange X and Y exactly.
    let x0: Vec<u64> = (0..v).map(|_| rng() & 1).collect();
    let y0: Vec<u64> = (0..v).map(|_| rng() & 1).collect();
    let (mut x, mut y) = (x0.clone(), y0.clone());
    bit.exchange_bits(&mut x, &mut y);
    println!(
        "bit exchange (h={h}): X' == Y and Y' == X : {}",
        if x == y0 && y == x0 { "PASS" } else { "FAIL" }
    );

    // Gaussian-dyadic case: (X, Y) -> (-Y, X).
    let cx = Motif::new(5, true);
    let v = cx.triple_count();
    let x0: Vec<Complex> = (0..v)
        .map(|_| Complex::new((rng() % 21) as f64 - 10.0, (rng() % 21) as f64 - 10.0))
        .collect();
    let y0: Vec<Complex> = (0..v)
        .map(|_| Complex::new((rng() % 21) as f64 - 10.0, (rng() % 21) as f64 - 10.0))
        .collect();
    let (mut x, mut y) = (x0.clone(), y0.clone());
    cx.exchange_complex(&mut x, &mut y);
    let ok = (0..v).all(|i| (x[i] + y0[i]).abs() < 1e-9 && (y[i] - x0[i]).abs() < 1e-9);
    println!(
        "complex exchange (h=5): (X, Y) -> (-Y, X) : {}",
        if ok { "PASS" } else { "FAIL" }
    );
}
