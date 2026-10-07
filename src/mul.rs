//! Three multipliers and a dispatcher: schoolbook, Karatsuba, and a
//! single-prime NTT.
//!
//! The paper's algorithm is asymptotic and unreachable; this is the practical
//! exact multiplier the crate actually ships — and the kind of subroutine the
//! paper assumes for short inner products. `Algorithm::Auto` picks a path by
//! operand size, so nobody has to think about thresholds unless they want to
//! benchmark them.

use crate::limbs::BigUint;
use crate::ntt;

/// Multiplication algorithm selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Algorithm {
    /// Pick by operand size (schoolbook → Karatsuba → NTT).
    Auto,
    /// `O(n^2)` schoolbook on base-2^32 limbs.
    Schoolbook,
    /// `O(n^1.585)` Karatsuba.
    Karatsuba,
    /// Exact NTT over the Goldilocks prime (falls back if out of range).
    Ntt,
}

impl Algorithm {
    pub fn name(self) -> &'static str {
        match self {
            Algorithm::Auto => "auto",
            Algorithm::Schoolbook => "schoolbook",
            Algorithm::Karatsuba => "karatsuba",
            Algorithm::Ntt => "ntt",
        }
    }
}

// Thresholds picked by benchmarking, not by vibes.
/// Schoolbook threshold in limbs (min operand length).
const SCHOOLBOOK_LIMBS: usize = 24;
/// Karatsuba threshold before handing over to the NTT (max operand limbs).
const KARATSUBA_LIMBS: usize = 4096;

/// Multiplies with an explicit algorithm choice.
pub fn mul(a: &BigUint, b: &BigUint, alg: Algorithm) -> BigUint {
    if a.is_zero() || b.is_zero() {
        return BigUint::zero();
    }
    match alg {
        Algorithm::Auto => auto(a, b),
        Algorithm::Schoolbook => schoolbook(a, b),
        Algorithm::Karatsuba => karatsuba_big(a, b),
        Algorithm::Ntt => ntt_mul(a, b).unwrap_or_else(|| karatsuba_big(a, b)),
    }
}

/// Size-based dispatch.
pub fn auto(a: &BigUint, b: &BigUint) -> BigUint {
    if a.is_zero() || b.is_zero() {
        return BigUint::zero();
    }
    let min_limbs = a.limbs().len().min(b.limbs().len());
    let max_limbs = a.limbs().len().max(b.limbs().len());
    if min_limbs <= SCHOOLBOOK_LIMBS {
        schoolbook(a, b)
    } else if max_limbs <= KARATSUBA_LIMBS {
        karatsuba_big(a, b)
    } else {
        ntt_mul(a, b).unwrap_or_else(|| karatsuba_big(a, b))
    }
}

/// Reference `O(n^2)` multiplication.
pub fn schoolbook(a: &BigUint, b: &BigUint) -> BigUint {
    BigUint::from_limbs(schoolbook_slices(a.limbs(), b.limbs()))
}

fn schoolbook_slices(a: &[u32], b: &[u32]) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0u64; a.len() + b.len()];
    for (i, &x) in a.iter().enumerate() {
        let mut carry = 0u64;
        for (j, &y) in b.iter().enumerate() {
            let acc = out[i + j] + x as u64 * y as u64 + carry;
            out[i + j] = acc & 0xFFFF_FFFF;
            carry = acc >> 32;
        }
        let mut k = i + b.len();
        while carry != 0 {
            let acc = out[k] + carry;
            out[k] = acc & 0xFFFF_FFFF;
            carry = acc >> 32;
            k += 1;
        }
    }
    normalize(out.into_iter().map(|v| v as u32).collect())
}

/// Karatsuba multiplication on base-2^32 limbs.
pub fn karatsuba_big(a: &BigUint, b: &BigUint) -> BigUint {
    BigUint::from_limbs(karatsuba(a.limbs(), b.limbs()))
}

fn karatsuba(a: &[u32], b: &[u32]) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    if a.len().min(b.len()) <= SCHOOLBOOK_LIMBS {
        return schoolbook_slices(a, b);
    }
    let m = a.len().max(b.len()).div_ceil(2);
    let (a0, a1) = split_at(a, m);
    let (b0, b1) = split_at(b, m);
    let z0 = karatsuba(a0, b0);
    let z2 = karatsuba(a1, b1);
    let sa = add_slices(a0, a1);
    let sb = add_slices(b0, b1);
    let z1_full = karatsuba(&sa, &sb);
    let z1 = sub_slices(&sub_slices(&z1_full, &z0), &z2);
    let mut out = vec![0u32; a.len() + b.len()];
    add_shifted(&mut out, &z0, 0);
    add_shifted(&mut out, &z1, m);
    add_shifted(&mut out, &z2, 2 * m);
    normalize(out)
}

fn split_at(s: &[u32], at: usize) -> (&[u32], &[u32]) {
    if at >= s.len() {
        (s, &[])
    } else {
        s.split_at(at)
    }
}

fn normalize(mut v: Vec<u32>) -> Vec<u32> {
    while v.last() == Some(&0) {
        v.pop();
    }
    v
}

fn add_slices(a: &[u32], b: &[u32]) -> Vec<u32> {
    let n = a.len().max(b.len()) + 1;
    let mut out = Vec::with_capacity(n);
    let mut carry = 0u64;
    for i in 0..n {
        let x = *a.get(i).unwrap_or(&0) as u64;
        let y = *b.get(i).unwrap_or(&0) as u64;
        let s = x + y + carry;
        out.push(s as u32);
        carry = s >> 32;
    }
    normalize(out)
}

fn sub_slices(a: &[u32], b: &[u32]) -> Vec<u32> {
    debug_assert!(compare_slices(a, b) != std::cmp::Ordering::Less);
    let mut out = Vec::with_capacity(a.len());
    let mut borrow = 0i64;
    for i in 0..a.len() {
        let x = a[i] as i64;
        let y = *b.get(i).unwrap_or(&0) as i64;
        let mut d = x - y - borrow;
        if d < 0 {
            d += 1 << 32;
            borrow = 1;
        } else {
            borrow = 0;
        }
        out.push(d as u32);
    }
    normalize(out)
}

fn compare_slices(a: &[u32], b: &[u32]) -> std::cmp::Ordering {
    let mut i = a.len();
    let mut j = b.len();
    while i > 0 && a[i - 1] == 0 {
        i -= 1;
    }
    while j > 0 && b[j - 1] == 0 {
        j -= 1;
    }
    if i != j {
        return i.cmp(&j);
    }
    for k in (0..i).rev() {
        match a[k].cmp(&b[k]) {
            std::cmp::Ordering::Equal => continue,
            ord => return ord,
        }
    }
    std::cmp::Ordering::Equal
}

fn add_shifted(out: &mut [u32], v: &[u32], shift: usize) {
    let mut carry = 0u64;
    for (i, &x) in v.iter().enumerate() {
        let idx = shift + i;
        if idx >= out.len() {
            debug_assert_eq!(x, 0);
            continue;
        }
        let acc = out[idx] as u64 + x as u64 + carry;
        out[idx] = acc as u32;
        carry = acc >> 32;
    }
    let mut idx = shift + v.len();
    while carry != 0 && idx < out.len() {
        let acc = out[idx] as u64 + carry;
        out[idx] = acc as u32;
        carry = acc >> 32;
        idx += 1;
    }
    debug_assert_eq!(carry, 0, "add_shifted overflow");
}

/// Extracts little-endian fixed-width digits of `width` bits.
pub fn extract_digits(v: &BigUint, width: usize) -> Vec<u64> {
    let mask = (1u64 << width) - 1;
    let count = v.bit_len().div_ceil(width).max(1);
    let limbs = v.limbs();
    let mut out = Vec::with_capacity(count);
    let mut bit = 0usize;
    for _ in 0..count {
        let limb = bit / 32;
        let off = bit % 32;
        let mut window: u128 = 0;
        for k in 0..4 {
            if limb + k < limbs.len() {
                window |= (limbs[limb + k] as u128) << (32 * k);
            }
        }
        out.push(((window >> off) as u64) & mask);
        bit += width;
    }
    while out.last() == Some(&0) {
        out.pop();
    }
    out
}

/// Folds base-`2^width` carries so every output digit is `< 2^width`.
pub fn fold_carries(values: &[u64], width: usize) -> Vec<u64> {
    let mask = (1u64 << width) - 1;
    let mut out = Vec::with_capacity(values.len() + 2);
    let mut carry = 0u64;
    for &v in values {
        // u128 on purpose: a coefficient can sit just below the Goldilocks
        // prime, and v + carry would wrap a u64 for huge (adversarial) inputs.
        let s = v as u128 + carry as u128;
        out.push((s & mask as u128) as u64);
        carry = (s >> width) as u64;
    }
    while carry != 0 {
        out.push(carry & mask);
        carry >>= width;
    }
    normalize_digits(&mut out);
    out
}

fn normalize_digits(v: &mut Vec<u64>) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

/// Re-packs little-endian `width`-bit digits into a [`BigUint`].
pub fn from_digits(digits: &[u64], width: usize) -> BigUint {
    let mut limbs: Vec<u32> = Vec::with_capacity(digits.len() * width / 32 + 1);
    let mut acc: u128 = 0;
    let mut acc_bits = 0usize;
    for &d in digits {
        acc |= (d as u128) << acc_bits;
        acc_bits += width;
        while acc_bits >= 32 {
            limbs.push(acc as u32);
            acc >>= 32;
            acc_bits -= 32;
        }
    }
    if acc_bits > 0 {
        limbs.push(acc as u32);
    }
    BigUint::from_limbs(limbs)
}

/// NTT-based multiplication. Returns `None` when the operands exceed the
/// single-prime range (then the caller should fall back to Karatsuba).
pub fn ntt_mul(a: &BigUint, b: &BigUint) -> Option<BigUint> {
    if a.is_zero() || b.is_zero() {
        return Some(BigUint::zero());
    }
    for width in [16usize, 15usize] {
        let da = extract_digits(a, width);
        let db = extract_digits(b, width);
        let out_len = da.len() + db.len() - 1;
        if out_len.next_power_of_two() > ntt::MAX_LEN {
            continue;
        }
        let max_a = da.iter().copied().max().unwrap_or(0) as u128;
        let max_b = db.iter().copied().max().unwrap_or(0) as u128;
        let bound = da.len().min(db.len()) as u128 * max_a * max_b;
        if bound >= ntt::P as u128 {
            continue;
        }
        let conv = ntt::convolve(&da, &db);
        let folded = fold_carries(&conv, width);
        return Some(from_digits(&folded, width));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::limbs::BigUint;

    fn splitmix(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    fn random_bits(rng: &mut u64, bits: usize) -> BigUint {
        let mut v = BigUint::zero();
        for i in 0..bits {
            if splitmix(rng) & 1 == 1 {
                v.set_bit(i);
            }
        }
        v
    }

    #[test]
    fn known_product() {
        let a = BigUint::from_u64(12345);
        let b = BigUint::from_u64(67890);
        assert_eq!(mul(&a, &b, Algorithm::Auto).to_u64(), Some(838_102_050));
    }

    #[test]
    fn algorithms_agree_across_sizes() {
        let mut rng = 0xC0FFEE_u64;
        for bits in [1usize, 7, 63, 64, 65, 127, 512, 2047, 4096, 20000] {
            let a = random_bits(&mut rng, bits);
            let b = random_bits(&mut rng, bits);
            let reference = schoolbook(&a, &b);
            assert_eq!(karatsuba_big(&a, &b), reference, "karatsuba at {bits} bits");
            assert_eq!(ntt_mul(&a, &b).unwrap(), reference, "ntt at {bits} bits");
            assert_eq!(auto(&a, &b), reference, "auto at {bits} bits");
        }
    }

    #[test]
    fn ntt_threshold_digit_round_trip() {
        let mut rng = 42u64;
        let v = random_bits(&mut rng, 12345);
        for width in [15usize, 16] {
            let d = extract_digits(&v, width);
            let back = from_digits(&d, width);
            assert_eq!(back, v);
        }
    }

    #[test]
    fn mersenne_square() {
        // (2^k - 1)^2 = 2^{2k} - 2^{k+1} + 1
        for k in [64usize, 1000, 5000] {
            let m = BigUint::one().shl(k).sub(&BigUint::one());
            let expected = BigUint::one()
                .shl(2 * k)
                .sub(&BigUint::one().shl(k + 1))
                .add(&BigUint::one());
            assert_eq!(ntt_mul(&m, &m).unwrap(), expected, "k={k}");
        }
    }

    #[test]
    fn fold_carries_survives_coefficients_near_the_prime() {
        // Coefficients can legitimately sit just below the Goldilocks prime;
        // adding the running carry must not wrap a u64.
        let near = crate::ntt::P - 1;
        let folded = fold_carries(&[near, near, near], 15);
        assert_eq!(folded[0], near & ((1u64 << 15) - 1));
        // Round-trip through digit packing must still recover the input sum.
        let digits = fold_carries(&[near, 0, 0, 0], 15);
        let value = from_digits(&digits, 15);
        assert_eq!(value, BigUint::from_u128(near as u128));
    }
}
