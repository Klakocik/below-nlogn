# below-nlogn

[![CI](https://github.com/Klakocik/below-nlogn/actions/workflows/ci.yml/badge.svg)](https://github.com/Klakocik/below-nlogn/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Dependencies: 0](https://img.shields.io/badge/dependencies-0-brightgreen.svg)](Cargo.toml)
[![Speedup vs GMP: eventually](https://img.shields.io/badge/speedup_vs_GMP-eventually-blueviolet.svg)](#could-this-replace-gmp)

A zero-dependency Rust toolbox for the executable machinery of OpenAI's
[**Integer multiplication below n log n**](https://github.com/openai/math/tree/main/preprints/Integer-multiplication-below-n-log-n-September-23-2026)
preprint (September 23, 2026) — *plus* a genuinely working exact integer
multiplier that makes no impossible promises.

> **The result.** One fixed finite-alphabet multitape Turing machine multiplies
> two `n`-bit integers in `O(n (lg n)^(1-κ))` worst-case time with
> **`κ = 2^-182`**, refuting the Schönhage–Strassen `n log n` optimality
> conjecture in that model.
>
> **The fine print.** The paper says its constants and thresholds are
> extremely large and calls the construction a bit-complexity result. This
> README quantifies exactly how large "extremely large" is, in percent signs.

---

## The 60-second version

| Question | Answer |
|---|---|
| Is it asymptotically faster than `n log n`? | Yes, provably (assuming the preprint holds up to review). |
| Is it faster on any computer that exists? | No. See [the compatibility matrix](#could-this-replace-gmp). |
| What does this repo actually run? | The paper's finite gadgets, executed and tested; plus a practical exact multiplier with division and square root. |
| Should I replace GMP? | [Please read the matrix first.](#could-this-replace-gmp) Then keep using GMP. |

## What's in the box

| Paper section | Module | Status |
|---|---|---|
| §1, constants | `asymptotics` | κ = 2⁻¹⁸², all published parameters, speedup/crossover math |
| §2 streams | `layers`, `xor_net` | row splitting, packed selected-bit XOR semantics |
| §3 motif networks (Lemma 3.2) | `motif` | **executed** for small `h`: arbitrary-scratch cancellation + bank exchange |
| §4 shear factorization (Lemmas 4.3, 4.4) | `xor_net` | **exact** rational `A = E1·Π·E2`; shear with `rank(A)` interchanges |
| §4 Prop 4.10 tape schedule | — | documented only (`docs/paper-map.md`) |
| §5 butterfly layers | `layers` | `H0` layers, packed XOR; tape-cost lemma documented |
| §6 synthetic transforms | `transforms` | `C[y]/(y^r+1)`, DIF recursion, bit-reversed storage, `F⁺F⁻ = I/t` |
| §7 Gaussian resampling | `resample` | matrices `S,T,C,D,N`; Harvey–van der Hoeven identity; contraction bound |
| §8 assembly | `crt`, `mul`, `ntt` | CRT layout + practical exact multiplier (NTT); full assembly documented |
| §9 division & square root | `division` | exact Newton reciprocal, `fast_div_rem`, `isqrt` |
| §10 transposition | — | pure corollary, documented only |
| §11 grouped rectangles | `layers` | packed XOR semantics only |

"Documented only" is never hand-waving: `docs/paper-map.md` states what the
section proves, why reproducing it means building a fixed Turing machine with
`1.77×10^20` wires, and where the code stops.

---

## Could this replace GMP?

Short version: **eventually yes, practically no**, and "eventually" is doing
Olympic-level lifting. The paper's growth rate beats `n lg n` by a factor

```text
speedup(n) = 100 · (1 − (lg n)^(−κ))  %,   κ = 2^-182 ≈ 1.63×10^-55,
```

which for every size you can build is, to sufficient precision, **zero**.

### Part 1 — what you get today

Growth-rate advantage over `n lg n` at realistic (and progressively less
realistic) operand sizes. Percent signs are load-bearing.

| Your operand | `lg n` | Advantage over `n lg n` |
|---|---:|---:|
| RSA-2048 key | 2 048 | `1.24×10⁻⁵² %` |
| RSA-4096 key | 4 096 | `1.36×10⁻⁵² %` |
| 8 GiB of operands | 2³⁶ | `4.07×10⁻⁵² %` |
| 2 EiB of operands | 2⁶⁴ | `7.24×10⁻⁵² %` |
| one bit per atom in the observable universe | 10⁸⁰ | `3.01×10⁻⁵¹ %` |
| the Bekenstein information bound of the universe | 10¹²² | `4.58×10⁻⁵¹ %` |
| one bit per Planck volume | ≈10¹⁸⁵ | `6.96×10⁻⁵¹ %` |

So: RSA-2048 gains **0.00000000000000000000000000000000000000000000000000124 %**.
Certificate authorities are advised to update their roadmaps accordingly.

### Part 2 — what you would have to build

| If you want… | …your operands must be | Perspective |
|---|---|---|
| `0.00000000000000013 %` faster | `2^(10^(3.46×10^36))` bits | bit-length itself ≈ 10^(3.46×10^36); the number of bits dwarfs the atom count of the universe |
| `10⁻¹² %` | `2^(10^(2.66×10^40))` bits | still adorable by comparison |
| `0.01 %` | `2^(10^(2.66×10^50))` bits | you will need a second universe for the exponent |
| `1 %` | `2^(10^(2.68×10^52))` bits | we recommend starting smaller and working up |
| `10 %` | `2^(10^(2.80×10^53))` bits | the heat death of the universe is now a build-time dependency |
| `50 %` (growth-rate parity) | `2^(2^(2^182))` bits | the paper's idea of "sufficiently large n" |

The last row is where `(lg n)^κ = 2`, i.e. the paper's growth curve is half of
`n lg n` — *up to the paper's constant*. The paper states that constant is
enormous, so any actual wall-clock crossover is later still. We are comparing
growth rates here, not stopwatches.

### Part 3 — the compatibility matrix

| If you currently use | Asymptotically replaced when… | Recommendation |
|---|---|---|
| GMP `mpz_mul` | `n ≥ 2^(2^(2^182))` bits | keep GMP |
| FLINT `fmpz_mul` | same threshold | keep FLINT |
| Rust `num-bigint` / `rug` | same threshold | keep them (or use this crate for study) |
| Python `int` | same threshold | keep CPython |
| Go `math/big` | same threshold | keep Go |
| Java `BigInteger` | same threshold | keep the JVM |
| OpenSSL `BN_mul` (RSA) | same threshold; RSA-2048 gets 1.24×10⁻⁵² % off | keep OpenSSL |
| JavaScript `BigInt` | same threshold | keep the engine |

**What this does not claim.** GMP, FLINT and friends already implement
near-`n log n` multiplication, so the paper's result is precisely a statement
about eventually beating *them* asymptotically. "Eventually" is the part no
computer will ever reach. This repository is an educational, experimental
implementation; it is not competitive with any of the libraries above, and it
does not try to be.

**Migration guide.** 1. Keep using GMP. 2. Wait. 3. Keep waiting.
4. *(Optional)* wait for proton decay. 5. Update `Cargo.toml`.

---

## Quickstart

Prerequisites: a stable Rust toolchain (`rustup` is enough). Nothing else to
install — the crate has no dependencies.

```console
# 1. Run the whole test suite (46 unit + 6 integration + doctest).
$ cargo test

# 2. Multiply two hex integers. Hex works with or without the 0x prefix, and
#    the CLI cross-checks the result against a second algorithm, so
#    "verified = yes" is a real check rather than decoration.
$ cargo run --release -- mul deadbeefcafebabe1122334455667788 ffffffff00000001
product  = deadbeeeec50fbcf24d137760f42ff01bbbbbbbc55667788
operands = 128/64 bits, product = 192 bits
algorithm= auto (3.2µs)
verified = yes (against schoolbook)

# 3. Open the paper's toyboxes and print its numbers.
$ cargo run --release -- motif 6            # Section 3 cancellation network, bit case
$ cargo run --release -- motif 5 --complex  # the same network over Gaussian dyadics
$ cargo run --release -- info               # kappa, constants, speedup tables
$ cargo run --release --example paper_math  # the README tables, computed by code

# 4. Watch the scaling.
$ cargo bench                               # schoolbook vs Karatsuba vs NTT
```

As a library (this is the crate's doctest, so it is compiled and tested on
every build):

```rust
use below_nlogn::{BigUint, mul};

// Operands can come from hex, decimal, or raw bits.
let a = BigUint::from_hex("cafebabe1234567890abcdef").unwrap();
let b = BigUint::from_dec("123456789012345678901234567890").unwrap();
let product = mul::mul(&a, &b, mul::Algorithm::Auto);

// The Section 9 corollary: exact division and integer square root.
let (q, r) = below_nlogn::division::fast_div_rem(&product, &b);
assert_eq!(q, a);
assert!(r.is_zero());
let root = below_nlogn::division::isqrt(&product);
assert!(below_nlogn::mul::auto(&root, &root) <= product);
```

## Benchmarks: honest numbers only

`cargo bench` (Windows desktop, single core, release build), random operands:

| bits | schoolbook | Karatsuba | NTT | NTT ns/bit |
|---:|---:|---:|---:|---:|
| 4 096 | 0.013 ms | 0.016 ms | 0.076 ms | 18.6 |
| 16 384 | 0.208 ms | 0.134 ms | 0.326 ms | 19.9 |
| 65 536 | 3.246 ms | 1.153 ms | 1.538 ms | 23.5 |
| 262 144 | 52.775 ms | 10.019 ms | 6.887 ms | 26.3 |
| 1 048 576 | — | 96.537 ms | 30.907 ms | 29.5 |

The NTT's `ns/bit` drifts upward like `log n`, because it is an `n log n`-class
algorithm. That is the only crossover this repository will ever see.

## Why you will never observe the theorem

The ratio of the paper's bound to `n lg n` is `(lg n)^(−κ)`. For a factor-2
saving:

```text
(lg n)^κ = 2  ⟺  lg n = 2^(1/κ) = 2^(2^182)  ⟺  n = 2^(2^(2^182)) bits,
```

i.e. `lg n ≈ 10^(1.85×10^54)`: the bit-length *alone* has about
`1.85×10^54` decimal digits. Even at `n = 2^(2^60)` bits the ratio is
`1.000000000000000` in `f64` — as `cargo run --release --example paper_math`
demonstrates. The theorem's victory is qualitative: the ratio tends to zero,
eventually. The universe supplies the "eventually".

## Layout

```text
src/
├── limbs.rs        arbitrary-precision integers (zero dependencies)
├── ntt.rs          exact convolution mod 2^64 - 2^32 + 1 (Goldilocks)
├── mul.rs          schoolbook / Karatsuba / NTT dispatch
├── fft.rs          complex FFT + Bluestein chirp transform
├── division.rs     Newton reciprocal, exact division, isqrt
├── crt.rs          CRT (Good–Thomas) convolution-preserving layout
├── motif.rs        Section 3 scalar motif networks (small-h, executable)
├── xor_net.rs      XOR swaps, three-shear interchange, exact factorization
├── layers.rs       H0 butterfly layers, packed selected-bit XOR
├── transforms.rs   synthetic transforms over C[y]/(y^r+1)
├── resample.rs     Gaussian resampling matrices and identity checks
└── asymptotics.rs  κ = 2^-182, parameters, speedup/crossover math
examples/  motif_section3, gaussian_resampling, paper_math
benches/   bench.rs (custom harness, no dependencies)
tests/     integration.rs
papers/    the source paper (PDF + LaTeX) and attribution, from openai/math
docs/      paper-map.md, architecture.md, honest-faq.md
```

## Testing

```console
$ cargo test
test result: ok. 46 passed; 0 failed; 1 ignored     # unit tests
test result: ok.  6 passed; 0 failed                 # integration tests
test result: ok.  1 passed; 0 failed                 # doctest
```

Highlights: cross-algorithm multiplication agreement and Mersenne-square
identities; exact rational `A = E1·Π·E2` reconstruction; Lemma 3.2 exchange in
both `F₂` and `Z[i,½]`; `F⁺F⁻ = I/t`; the Harvey–van der Hoeven identity
`T P_s F_s = P_t F_t S` and the composed resampling identity; Newton division
against reference long division; exact `(s+1)² > a ≥ s²` square-root bounds;
and speedup percentages pinned to the tables above (so the README cannot
quietly drift).

## Contributing

Issues and PRs are welcome. See `CONTRIBUTING.md`. The bar for new code:
tests against a naive reference, no new dependencies without discussion, and
honesty in comments — if something is a numerical demonstration rather than a
proof, say so.

## Citation, attribution, license

* Cite the paper for the mathematics; cite this repository only for the
  executable gadgets and the honest tables.

```bibtex
@misc{OAI:Integer-multiplication-below-n-log-n-September-23-2026,
  author = {{OpenAI}},
  title = {{Integer multiplication below $n\log n$}},
  howpublished = {OpenAI Math Release preprint},
  year = {2026}
}
```

* The paper and its LaTeX sources in `papers/` are Apache-2.0, © OpenAI
  (see `papers/SOURCES.md`). The Goldilocks-prime NTT is standard folklore and
  not taken from the preprint. Everything else is Apache-2.0, © the
  `below-nlogn` contributors. See `LICENSE`.
* Not affiliated with OpenAI. Any transcription mistakes in this repository
  are ours, not the paper's.
