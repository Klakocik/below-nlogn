//! Integration tests over the public API: exactness of the multiplier on
//! random inputs, the corollary's division and square root, the CRT layout,
//! the motif exchange, and the Gaussian resampling identity.

use below_nlogn::crt::CrtLayout;
use below_nlogn::division;
use below_nlogn::limbs::BigUint;
use below_nlogn::motif::Motif;
use below_nlogn::mul::{self, Algorithm};
use below_nlogn::resample::Resampling;

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

#[test]
fn all_algorithms_agree_on_random_inputs() {
    let mut state = 0xABCD_1234_5678_9EF0u64;
    for bits in [1usize, 31, 200, 4096, 40000] {
        let a = random_bits(&mut state, bits);
        let b = random_bits(&mut state, bits);
        let reference = mul::schoolbook(&a, &b);
        assert_eq!(mul::mul(&a, &b, Algorithm::Karatsuba), reference, "{bits}");
        assert_eq!(mul::mul(&a, &b, Algorithm::Ntt), reference, "{bits}");
        assert_eq!(mul::mul(&a, &b, Algorithm::Auto), reference, "{bits}");
    }
}

#[test]
fn hex_round_trip_and_product() {
    let a = BigUint::from_hex("cafebabe1234567890abcdef").unwrap();
    let b = BigUint::from_hex("fedcba98765432100f0f0f0f").unwrap();
    let p = mul::auto(&a, &b);
    let q = mul::auto(&b, &a);
    assert_eq!(p, q);
    assert_eq!(BigUint::from_hex(&p.to_hex()).unwrap(), p);
}

#[test]
fn division_and_sqrt_corollary() {
    let mut state = 0x1111_2222_3333_4444u64;
    let a = random_bits(&mut state, 2048);
    let b = random_bits(&mut state, 900).add(&BigUint::one());
    let (q, r) = division::fast_div_rem(&a, &b);
    assert!(r < b);
    assert_eq!(mul::auto(&q, &b).add(&r), a);
    let s = division::isqrt(&a);
    assert!(mul::auto(&s, &s) <= a);
    assert!(mul::auto(&s.add(&BigUint::one()), &s.add(&BigUint::one())) > a);
}

#[test]
fn crt_layout_round_trip() {
    let dims = [5usize, 7, 9, 11];
    let layout = CrtLayout::new(&dims);
    let data: Vec<u64> = (0..layout.len as u64).map(|i| i * i + 7).collect();
    let tensor = layout.to_tensor(&data);
    assert_eq!(layout.to_flat(&tensor), data);
}

#[test]
fn motif_exchange_public_api() {
    let motif = Motif::new(5, false);
    let v = motif.triple_count();
    let x0: Vec<u64> = (0..v).map(|i| (i % 2) as u64).collect();
    let y0: Vec<u64> = (0..v).map(|i| ((i / 2) % 2) as u64).collect();
    let (mut x, mut y) = (x0.clone(), y0.clone());
    motif.exchange_bits(&mut x, &mut y);
    assert_eq!(x, y0);
    assert_eq!(y, x0);
}

#[test]
fn resampling_identity_small_case() {
    let r = Resampling::new(11, 13, 1.5);
    assert!(r.identity_error(777) < 1e-6);
    assert!(r.n_minus_i_norm() < 0.5);
}
