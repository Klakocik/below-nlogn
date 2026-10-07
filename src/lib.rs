//! # below-nlogn
//!
//! An experimental, dependency-free Rust implementation of the executable
//! gadgets from the OpenAI preprint
//! [*Integer multiplication below n log n*](https://github.com/openai/math/tree/main/preprints/Integer-multiplication-below-n-log-n-September-23-2026)
//! (September 23, 2026).
//!
//! The paper proves that a fixed finite-alphabet multitape Turing machine can
//! multiply two `n`-bit integers in
//! `O(n (lg n)^(1 - κ))` worst-case time with `κ = 2^-182`, refuting the
//! Schönhage–Strassen `n log n` optimality conjecture in that model.
//!
//! **Honest positioning.** The paper is a bit-complexity result: its constants
//! and thresholds are astronomically large and it does not describe a practical
//! multiplier. This crate does not claim to be `n log n`-beating at any
//! realizable size — if you want the fastest multiplier on Earth, GMP is a few
//! doors down and has decades of a head start. What this crate does:
//!
//! * implements and tests the *finite, executable* pieces of the construction
//!   (the Section 3 scalar motif network, XOR/shear interchanges, packed
//!   selected-bit operations, CRT layouts, the synthetic ring transform,
//!   Gaussian resampling numerics, Bluestein's chirp transform);
//! * ships a working exact integer multiplier (schoolbook, Karatsuba, and a
//!   Goldilocks-prime NTT), together with the corollary's division and integer
//!   square root;
//! * records the paper's parameters and exposes the crossover analysis showing
//!   why the asymptotic improvement is unobservable in practice.
//!
//! Every module documents the paper section it corresponds to. See
//! `docs/paper-map.md` for the mapping and `docs/honest-faq.md` for caveats.
//!
//! ## Quickstart
//!
//! ```rust
//! use below_nlogn::{BigUint, mul};
//!
//! // Operands can come from hex, decimal, or raw bits.
//! let a = BigUint::from_hex("cafebabe1234567890abcdef").unwrap();
//! let b = BigUint::from_dec("123456789012345678901234567890").unwrap();
//! let product = mul::mul(&a, &b, mul::Algorithm::Auto);
//!
//! // The Section 9 corollary: exact division and integer square root.
//! let (q, r) = below_nlogn::division::fast_div_rem(&product, &b);
//! assert_eq!(q, a);
//! assert!(r.is_zero());
//! let root = below_nlogn::division::isqrt(&product);
//! assert!(below_nlogn::mul::auto(&root, &root) <= product);
//! ```

pub mod asymptotics;
pub mod crt;
pub mod division;
pub mod fft;
pub mod layers;
pub mod limbs;
pub mod motif;
pub mod mul;
pub mod ntt;
pub mod resample;
pub mod transforms;
pub mod xor_net;

pub use limbs::BigUint;
