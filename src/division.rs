//! Division and integer square root, the paper's Section 9 corollary.
//!
//! Both are Newton iterations on integers. The reciprocal doubles its correct
//! bits until the exact test `2^p - c*x < c` says stop; division then corrects
//! at most a few ulps; `isqrt` is the classic `x <- (x + a/x)/2` with the same
//! belt-and-braces correction at the end. Slow to write, hard to get wrong,
//! boring in the best sense.

use crate::limbs::BigUint;
use crate::mul;

/// Exact `(quotient, remainder)` from the reference long division.
pub fn div_rem(a: &BigUint, b: &BigUint) -> (BigUint, BigUint) {
    a.div_rem(b)
}

/// `floor(2^p / c)` by Newton iteration, with exact final corrections.
///
/// The iteration keeps `x` slightly below `2^scale / c` and squares the
/// relative error; it needs `O(1)` multiplications of geometrically growing
/// precision, so the total is `O(M(n))` with `M` the multiplication cost.
pub fn reciprocal_floor(c: &BigUint, p: usize) -> BigUint {
    assert!(!c.is_zero(), "reciprocal of zero");
    let one = BigUint::one();
    if c.is_one() {
        return one.shl(p);
    }
    let n = c.bit_len();
    // Small divisors: exact u128 arithmetic, no Newton required. Careful with
    // the cast: a quotient can reach 2^125 and does not fit in u64.
    if let Some(cv) = c.to_u64() {
        if p <= 126 {
            let t = 1u128 << p;
            return BigUint::from_u128(t / cv as u128);
        }
    }
    if p < n {
        // c >= 2^(n-1) > 2^p, so the quotient is 0, except in the boundary
        // case where c is exactly 2^p.
        return if *c == one.shl(p) {
            one
        } else {
            BigUint::zero()
        };
    }
    if p == n {
        // c is in [2^(n-1), 2^n): quotient 1, or exactly 2 when c = 2^(n-1).
        return if *c == one.shl(n - 1) {
            BigUint::from_u64(2)
        } else {
            one
        };
    }
    // Seed at scale n: floor(2^n / c) = 1 because c in [2^{n-1}, 2^n);
    // the relative error is below 1/2.
    let mut scale = n;
    let mut x = one.clone();
    while scale < p {
        let new_scale = (2 * scale).min(p);
        // Lift to the new scale (the relative error is preserved)...
        x = x.shl(new_scale - scale);
        // ...then one Newton step at the new scale squares the relative error.
        let target = one.shl(new_scale);
        let big = one.shl(new_scale + 1);
        let mut guard = 0;
        while mul::auto(c, &x) > big {
            x = x.sub(&one);
            guard += 1;
            assert!(guard < 1 << 20, "Newton lift clamp failed");
        }
        let cx = mul::auto(c, &x);
        let numerator = mul::auto(&x, &big.sub(&cx));
        let mut x_next = numerator.shr(new_scale);
        let mut guard = 0;
        while mul::auto(c, &x_next) > target {
            x_next = x_next.sub(&one);
            guard += 1;
            assert!(guard < 1 << 20, "Newton step clamp failed");
        }
        x = x_next;
        scale = new_scale;
    }
    // Refine at the final scale until x is exactly floor(2^p / c):
    // Newton squarings at fixed scale converge in O(log p) steps, and the
    // exact test `2^p - c x < c` terminates the loop at the true floor.
    let target = one.shl(p);
    let big = target.shl(1);
    let mut guard = 0;
    loop {
        let cx = mul::auto(c, &x);
        if cx <= target {
            if target.sub(&cx) < *c {
                break;
            }
        } else {
            x = x.sub(&one);
            guard += 1;
            assert!(guard < 1 << 20, "reciprocal refinement failed");
            continue;
        }
        let numerator = mul::auto(&x, &big.sub(&cx));
        let x_new = numerator.shr(p);
        if x_new > x {
            x = x_new;
        } else {
            // Newton can stall on the last ulp (x = floor - 1); one exact
            // increment reaches the floor and the test above will stop.
            x = x.add(&one);
        }
        guard += 1;
        assert!(guard < 1 << 20, "reciprocal refinement failed");
    }
    x
}

/// Fast exact division using a Newton reciprocal with 64 guard bits.
pub fn fast_div_rem(a: &BigUint, b: &BigUint) -> (BigUint, BigUint) {
    assert!(!b.is_zero(), "division by zero");
    let one = BigUint::one();
    if a < b {
        return (BigUint::zero(), a.clone());
    }
    if b.is_one() {
        return (a.clone(), BigUint::zero());
    }
    if b.bit_len() <= 64 {
        // Fall back on the reference path for small divisors.
        return a.div_rem(b);
    }
    let guard_bits = 64;
    let p = a.bit_len() + guard_bits;
    let recip = reciprocal_floor(b, p);
    let mut q = mul::auto(a, &recip).shr(p);
    // At most a few corrections.
    let mut guard = 0;
    while mul::auto(&q, b) > *a {
        q = q.sub(&one);
        guard += 1;
        assert!(guard < 64, "quotient correction failed");
    }
    guard = 0;
    loop {
        let q1 = q.add(&one);
        if mul::auto(&q1, b) <= *a {
            q = q1;
            guard += 1;
            assert!(guard < 64, "quotient correction failed");
        } else {
            break;
        }
    }
    let r = a.sub(&mul::auto(&q, b));
    (q, r)
}

/// Exact integer square root: `floor(sqrt(a))`.
pub fn isqrt(a: &BigUint) -> BigUint {
    if a.is_zero() {
        return BigUint::zero();
    }
    let one = BigUint::one();
    if *a < one.shl(2) {
        return one; // a in {1, 2, 3}
    }
    // Newton from above: x0 = 2^{ceil(bits/2)} >= sqrt(a).
    let mut x = one.shl(a.bit_len().div_ceil(2));
    loop {
        let (q, _) = fast_div_rem(a, &x);
        let y = x.add(&q).shr(1);
        if y >= x {
            break;
        }
        x = y;
    }
    // Exact corrections.
    let mut guard = 0;
    while mul::auto(&x, &x) > *a {
        x = x.sub(&one);
        guard += 1;
        assert!(guard < 64, "isqrt correction failed");
    }
    guard = 0;
    loop {
        let x1 = x.add(&one);
        if mul::auto(&x1, &x1) <= *a {
            x = x1;
            guard += 1;
            assert!(guard < 64, "isqrt correction failed");
        } else {
            break;
        }
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

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
    #[ignore = "diagnostic"]
    fn debug_reciprocal() {
        let mut rng = 4242u64;
        let c = random_bits(&mut rng, 200).add(&BigUint::one());
        let n = c.bit_len();
        let p = 512usize;
        let one = BigUint::one();
        let mut scale = n;
        let mut x = one.clone();
        println!("n={n} p={p}");
        while scale < p {
            let new_scale = (2 * scale).min(p);
            x = x.shl(new_scale - scale);
            let target = one.shl(new_scale);
            let big = one.shl(new_scale + 1);
            let mut guard = 0;
            while mul::auto(&c, &x) > big {
                x = x.sub(&one);
                guard += 1;
            }
            if guard > 0 {
                println!("clamp at scale {new_scale}: {guard} steps");
            }
            let cx = mul::auto(&c, &x);
            let numerator = mul::auto(&x, &big.sub(&cx));
            let mut x_next = numerator.shr(new_scale);
            let mut guard = 0;
            while mul::auto(&c, &x_next) > target {
                x_next = x_next.sub(&one);
                guard += 1;
            }
            // Relative error in units of 2^-50.
            let err = mul::auto(&c, &x_next);
            let err = if err > target {
                err.sub(&target)
            } else {
                target.sub(&err)
            };
            println!(
                "scale {scale} -> {new_scale}: clamp {guard}, err/2^scale ≈ 2^{:.1}",
                (err.bit_len() as f64) - (new_scale as f64)
            );
            x = x_next;
            scale = new_scale;
        }
        let got = reciprocal_floor(&c, p);
        let expected = one.shl(p).div_rem(&c).0;
        println!("equal = {}", got == expected);
    }

    #[test]
    fn reciprocal_matches_reference() {
        let mut rng = 12345u64;
        for bits in [8usize, 63, 64, 65, 200, 1000, 4096] {
            let c = random_bits(&mut rng, bits);
            if c.is_zero() || c.is_one() {
                continue;
            }
            for p in [10usize, 64, 65, 256, 1000] {
                if p == 0 {
                    continue;
                }
                let got = reciprocal_floor(&c, p);
                let expected = BigUint::one().shl(p).div_rem(&c).0;
                assert_eq!(got, expected, "c bits={bits} p={p}");
            }
        }
    }

    #[test]
    fn reciprocal_edge_cases() {
        let one = BigUint::one();
        // p == bit_len(c) with c a power of two: floor(2^n / 2^(n-1)) = 2.
        for n in [2usize, 5, 64, 200] {
            let c = one.shl(n - 1);
            assert_eq!(reciprocal_floor(&c, n), BigUint::from_u64(2), "n={n}");
        }
        // p below bit_len(c): zero, except when c is exactly 2^p.
        let c = one.shl(100);
        assert_eq!(reciprocal_floor(&c, 99), BigUint::zero());
        assert_eq!(reciprocal_floor(&c, 100), BigUint::from_u64(1));
        // Small divisor, high precision: the quotient must not be truncated
        // to 64 bits by the u128 fast path.
        let c = BigUint::from_u64(3);
        let got = reciprocal_floor(&c, 126);
        assert_eq!(got, one.shl(126).div_rem(&c).0);
        assert!(got.bit_len() > 64);
    }

    #[test]
    fn fast_division_matches_reference() {
        let mut rng = 999u64;
        for bits in [100usize, 300, 1024, 5000] {
            let a = random_bits(&mut rng, bits * 2);
            let b = random_bits(&mut rng, bits);
            if b.is_zero() {
                continue;
            }
            let (q1, r1) = fast_div_rem(&a, &b);
            let (q2, r2) = div_rem(&a, &b);
            assert_eq!(q1, q2, "quotient bits={bits}");
            assert_eq!(r1, r2, "remainder bits={bits}");
        }
    }

    #[test]
    fn isqrt_is_exact() {
        let mut rng = 7u64;
        for bits in [4usize, 16, 64, 129, 1024, 4093] {
            let a = random_bits(&mut rng, bits);
            let s = isqrt(&a);
            assert!(mul::auto(&s, &s) <= a, "s^2 <= a at bits={bits}");
            let s1 = s.add(&BigUint::one());
            assert!(mul::auto(&s1, &s1) > a, "(s+1)^2 > a at bits={bits}");
        }
        for k in [1usize, 2, 3, 100, 10_000] {
            let x = random_bits(&mut rng, k).add(&BigUint::one());
            let sq = mul::auto(&x, &x);
            assert_eq!(isqrt(&sq), x, "perfect square k={k}");
        }
    }
}
