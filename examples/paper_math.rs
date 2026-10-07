//! The paper's headline constants, the speedup percentages, and the crossover
//! analysis, in one place.
//!
//! Run with: `cargo run --release --example paper_math`

use below_nlogn::asymptotics;

fn main() {
    println!("Theorem (OpenAI, Sept 23 2026):");
    println!("  T(n) = O(n (lg n)^(1-kappa)),  kappa = 2^-182,");
    println!("  on one fixed finite-alphabet multitape Turing machine.");
    println!();
    println!("{}", asymptotics::kappa_summary());
    println!();
    println!("Published parameters:");
    for (name, value) in asymptotics::PAPER_CONSTANTS {
        println!("  {name:<26} {value}");
    }
    println!();

    // Growth-rate advantage over n lg n: 100 * (1 - (lg n)^-kappa).
    println!("Growth-rate advantage over n lg n at named sizes:");
    for (lg_n, label) in [
        (2048.0f64, "RSA-2048"),
        (4096.0, "RSA-4096"),
        (2f64.powi(36), "8 GiB operand"),
        (2f64.powi(64), "2 EiB operand"),
        (1e80, "10^80 bits (1 per atom, observable universe)"),
        (1e122, "10^122 bits (Bekenstein information bound)"),
        (2e185, "~10^185 bits (1 per Planck volume)"),
    ] {
        println!("  {label:<46} {:e} %", asymptotics::speedup_percent(lg_n));
    }
    println!();

    println!("Input size required for a target growth-rate advantage:");
    for (percent, label) in [
        (1.3e-16f64, "0.00000000000000013%"),
        (1e-12, "10^-12%"),
        (1e-4, "10^-4%"),
        (0.01, "0.01%"),
        (1.0, "1%"),
        (10.0, "10%"),
    ] {
        let log10_lg = asymptotics::log10_lg_n_for_percent(percent);
        println!("  {label:<24} lg n = 10^{log10_lg:e},  n = 2^(10^{log10_lg:e}) bits");
    }
    let parity = asymptotics::log10_lg_n_for_percent(50.0);
    println!(
        "  {:<24} lg n = 2^(2^182) = 10^{parity:e},  n = 2^(2^(2^182)) bits",
        "50% (parity)"
    );
    println!();

    println!("{}", asymptotics::crossover_description(1));
    println!();
    println!("For any n you could ever store, the ratio (lg n)^-kappa is 1 to all digits");
    println!("that fit in f64. The theorem's victory is qualitative: it tends to zero.");
}
