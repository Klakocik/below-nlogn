//! The numbers that keep this project honest: κ = 2^-182, the paper's chosen
//! constants, and the speedup arithmetic.
//!
//! The growth-rate advantage over `n lg n` is `1 - (lg n)^{-κ}`. At RSA-2048
//! that is `1.24e-52 %`, and even at the Bekenstein bound of the observable
//! universe it is `4.58e-51 %`. [`speedup_percent`] and
//! [`log10_lg_n_for_percent`] do the arithmetic carefully (`exp_m1` for the
//! tiny values, log space for the huge ones) so the README tables cannot
//! quietly drift from the math.

use crate::limbs::BigUint;

/// `κ = 2^-182`, the exponent improvement of the main theorem.
pub const KAPPA_EXPONENT: u32 = 182;

/// `τ = σ = 1 - 2^-50` (motif exponents, eq. (explicit-motif-exponents)).
pub const TAU_DENOMINATOR_EXPONENT: u32 = 50;

/// `λ = 1 - 2^-52` (layer exponent).
pub const LAMBDA_DENOMINATOR_EXPONENT: u32 = 52;

/// `λ' = 1 - 2^-54` (simultaneous butterfly-layer exponent).
pub const LAMBDA_PRIME_DENOMINATOR_EXPONENT: u32 = 54;

/// `c = 2^-56` (K exponent).
pub const C_EXPONENT: u32 = 56;

/// `ε = 2^-75` (dimension exponent).
pub const EPSILON_EXPONENT: u32 = 75;

/// `C_1 = 20` (layer constant).
pub const C1: u32 = 20;

/// The published constants, as decimal strings where they are large.
pub const PAPER_CONSTANTS: &[(&str, &str)] = &[
    ("κ", "2^-182"),
    ("τ = σ", "1 - 2^-50"),
    ("λ", "1 - 2^-52"),
    ("λ'", "1 - 2^-54"),
    ("c", "2^-56"),
    ("ε", "2^-75"),
    ("C_1", "20"),
    ("β", "1/2"),
    ("δ", "1/16"),
    ("h (motif ground set)", "100"),
    ("v = C(h,3)", "161700"),
    ("N = v^3", "4227952113000000"),
    ("W_b (bit wires)", "177176569091445000000"),
    ("W_c (complex wires)", "1873807244643542670000"),
    ("η_b (bit rank deficit)", "339/22587335000000"),
    ("η_c (complex rank deficit)", "73/19906842167500"),
];

/// Decimal digits of `κ = 5^182 · 10^-182`; returns `(mantissa_digits, exponent)`.
pub fn kappa_decimal() -> (String, i32) {
    let mut pow = BigUint::one();
    for _ in 0..KAPPA_EXPONENT {
        pow = pow.mul_small(5);
    }
    let digits = pow.to_dec();
    (digits, -(KAPPA_EXPONENT as i32))
}

/// Human-readable summary of how small `κ` is.
pub fn kappa_summary() -> String {
    let (digits, _) = kappa_decimal();
    // κ = 5^182 · 10^-182, and 5^182 has len(digits) decimal digits, so
    // κ = 1.[digits] × 10^(-182 + len - 1).
    let leading = &digits[..1];
    let mut frac = digits[1..].to_string();
    if frac.len() > 12 {
        frac.truncate(12);
    }
    let exponent = -(KAPPA_EXPONENT as i32) + (digits.len() as i32 - 1);
    format!("κ = 2^-182 = 5^182 · 10^-182 ≈ {leading}.{frac} × 10^{exponent}")
}

/// For a requested asymptotic saving factor `2^k` in the log exponent,
/// returns a description of the smallest `n` at which the paper's bound
/// improves on `n lg n` by that factor.
///
/// The saving is `(lg n)^κ = 2^k` iff `lg n = 2^{k/κ} = 2^{k·2^182}`.
pub fn crossover_description(saving_log2: u32) -> String {
    if saving_log2 == 0 {
        return "k = 0: no saving (identical log exponent).".to_string();
    }
    let a = saving_log2 as f64 * 2f64.powi(KAPPA_EXPONENT as i32); // k·2^182
    let log10_lg_n = a * 2f64.log10();
    format!(
        "saving 2^-k with k = {saving_log2}: need lg n = 2^(k·2^182) = 2^{a:.3e}, \
         so log10(lg n) ≈ {log10_lg_n:.3e}: the logarithm itself has about \
         {log10_lg_n:.3e} decimal digits, and n has 2^(k·2^182) bits. Far beyond \
         any physical or Gedanken computation."
    )
}

/// The asymptotic time bound `n (lg n)^{1-κ}` evaluated with both the paper's
/// exponent and the classical `n lg n`, as a ratio, for a given bit length.
pub fn bound_ratio(n_bits: f64) -> f64 {
    let lg = n_bits.max(2.0).log2();
    lg.powf(-(2f64.powi(-(KAPPA_EXPONENT as i32))))
}

/// Relative growth-rate advantage of the paper's bound over `n lg n` for a
/// given value of `lg n`: `1 - (lg n)^{-κ}`.
///
/// Computed with `exp_m1` so that tiny values keep their relative precision
/// (the naive expression rounds to exactly zero for every realistic size).
pub fn size_advantage(lg_n: f64) -> f64 {
    let kappa = 2f64.powi(-(KAPPA_EXPONENT as i32));
    -(-kappa * lg_n.ln()).exp_m1()
}

/// [`size_advantage`] expressed as a percentage, e.g. `4.07e-52` (%) for an
/// 8 GiB operand and `1.24e-52` (%) for an RSA-2048 operand.
pub fn speedup_percent(lg_n: f64) -> f64 {
    100.0 * size_advantage(lg_n)
}

/// `log10(lg n)` required for a target percentage saving `target_percent`.
///
/// Inverting `p = 100 (1 - e^{-κ ln lg n})` gives
/// `lg n = exp(-ln(1 - p/100) / κ)`; the result is returned in log space so
/// that targets like `1.3e-16` (%) do not overflow `f64`.
pub fn log10_lg_n_for_percent(target_percent: f64) -> f64 {
    let kappa = 2f64.powi(-(KAPPA_EXPONENT as i32));
    let x = -(-target_percent / 100.0).ln_1p(); // -ln(1 - p/100)
    x / kappa * std::f64::consts::LOG10_E
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kappa_decimal_is_correct_magnitude() {
        let (digits, exp) = kappa_decimal();
        assert_eq!(exp, -182);
        // 5^182 has 128 decimal digits and starts with 1.
        assert_eq!(digits.len(), 128, "digits={digits}");
        assert!(digits.starts_with('1'));
        // κ ≈ 1.63 × 10^-55.
        assert!(kappa_summary().contains("10^-55"));
    }

    #[test]
    fn bound_ratio_is_small_but_indistinguishable() {
        // Even at an absurd n = 2^1e15 bits the saving is negligible.
        let r = bound_ratio(2f64.powi(50));
        assert!(r > 0.999);
        assert!(r <= 1.0);
    }

    #[test]
    fn speedup_percentages_match_the_readme() {
        // Growth-rate advantage over n lg n, as a percentage.
        let close = |got: f64, want: f64| (got - want).abs() <= want * 1e-3;
        assert!(
            close(speedup_percent(2048.0), 1.2438e-52),
            "RSA-2048: {:e}",
            speedup_percent(2048.0)
        );
        assert!(
            close(speedup_percent(2f64.powi(36)), 4.0707e-52),
            "8 GiB: {:e}",
            speedup_percent(2f64.powi(36))
        );
        assert!(
            close(speedup_percent(2f64.powi(64)), 7.2368e-52),
            "2 EiB: {:e}",
            speedup_percent(2f64.powi(64))
        );
        assert!(
            close(speedup_percent(1e80), 3.0050e-51),
            "10^80 bits: {:e}",
            speedup_percent(1e80)
        );
        assert!(
            close(speedup_percent(1e122), 4.5826e-51),
            "10^122 bits: {:e}",
            speedup_percent(1e122)
        );
        // Growth-rate parity with n lg n (advantage 1/2) requires
        // lg n = 2^(1/κ) = 2^(2^182), which overflows f64; verify the
        // log-space inversion instead.
        let parity_log10 = log10_lg_n_for_percent(50.0);
        assert!(
            (parity_log10 - 1.8453085040877162e54).abs() <= 1e48,
            "parity log10: {:e}",
            parity_log10
        );
    }

    #[test]
    fn required_sizes_for_target_percentages() {
        // 0.00000000000000013% needs lg n = 10^(3.46e36).
        let l = log10_lg_n_for_percent(1.3e-16);
        assert!((l - 3.461e36).abs() < 1e34, "got {l:e}");
        // 1% needs lg n = 10^(2.68e52).
        let l = log10_lg_n_for_percent(1.0);
        assert!((l - 2.676e52).abs() < 1e50, "got {l:e}");
        // The formula is monotone.
        assert!(log10_lg_n_for_percent(0.1) < log10_lg_n_for_percent(1.0));
    }
}
