//! Numerical demonstration of Gaussian resampling (Section 7):
//!
//! * the Harvey–van der Hoeven identity `T P_s F_s = P_t F_t S`
//!   (Theorem 4.2 of the cited manuscript, operator form `A = S/2`,
//!   `B_0 = D' J' C` of Lemma 7.2);
//! * the contraction estimate `||N - I|| < 2.01 e^{-pi alpha^2 theta / 2}`.
//!
//! The paper's versions are exact fixed-point algorithms with a rigorous
//! error analysis; this example is a float64 numerical check of the interface.
//!
//! Run with: `cargo run --release --example gaussian_resampling`

use below_nlogn::resample::Resampling;

fn main() {
    println!("Gaussian resampling (Section 7) — numerical checks");
    println!();
    println!(
        "{:>4} {:>4} {:>5} {:>10} {:>14} {:>14} {:>14}",
        "s", "t", "alpha", "theta", "||N-I||", "bound", "identity err"
    );
    for (s, t, alpha) in [(17usize, 33usize, 2.0f64), (53, 59, 4.0), (11, 13, 1.5)] {
        let r = Resampling::new(s, t, alpha);
        println!(
            "{s:>4} {t:>4} {alpha:>5.1} {:>10.4} {:>14.6} {:>14.6} {:>14.2e}",
            r.theta(),
            r.n_minus_i_norm(),
            r.contraction_bound(),
            r.identity_error(0xBEEF)
        );
    }
    println!();
    println!("Interpretation: N = C T D is the normalised version of the Gaussian map T;");
    println!("||N - I|| < 1/2 guarantees invertibility by the Neumann series, which is how");
    println!("the paper synthesises the inverse resampling contraction B_0.");
}
