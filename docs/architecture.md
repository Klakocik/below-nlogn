# Architecture

`below-nlogn` is a single Rust crate with **zero dependencies**. Everything —
big integers, the NTT, complex arithmetic, exact rational matrices — is
implemented locally so the repository is self-contained and reproducible on
any stable Rust toolchain.

```text
below-nlogn
├── limbs        BigUint on Vec<u32> (little-endian, normalized)
├── ntt          Goldilocks-prime NTT + exact convolution
├── mul          schoolbook / Karatsuba / NTT dispatch
├── fft          complex FFT, naive DFT, Bluestein chirp transform
├── division     Newton reciprocal, fast exact division, isqrt
├── crt          Good–Thomas CRT layout, multidimensional cyclic convolution
├── motif        Section 3 scalar networks (bit and Gaussian-dyadic)
├── xor_net      XOR swaps, three-shear interchange, rational Bruhat factorization
├── layers       H0 butterfly layers, packed selected-bit XOR
├── transforms   synthetic transforms over C[y]/(y^r+1)
├── resample     Gaussian resampling matrices, HvH identity, contraction bound
└── asymptotics  constants of the paper + crossover analysis
```

## Exact multiplier

`BigUint` stores base-`2^32` limbs. `mul::mul` dispatches by operand size:

* **Schoolbook** (`O(n²)`) below 24 limbs on the shorter operand.
* **Karatsuba** (`O(n^1.585)`) up to 4096 limbs.
* **NTT** above that: digits of 16 bits (15 if needed), convolution modulo
  `p = 2^64 - 2^32 + 1` (Goldilocks; `p - 1 = 2^32 · (2^32 - 1)`, primitive
  root 7), then carry folding base `2^16`. The bound
  `min(len_a,len_b) · max_a · max_b < p` is checked before transforming, and
  the code falls back to Karatsuba if a single prime would not suffice.

Why an NTT instead of the paper's complex FFT? The paper needs approximate
complex arithmetic because a Turing machine cannot cheaply do modular
arithmetic on wide words; Rust on a 64-bit CPU can. The NTT is exact, hence
no rounding proof is needed for the shipped multiplier. The complex machinery
still exists in `fft`, `transforms`, `layers` and `resample` because those
modules study the paper's construction itself.

## Executed gadget models

* `motif` enumerates `C(h,3)` triples and runs the eight-row schedule for
  small `h`, in `F_2` and in `Z[i,1/2]`. It verifies the cancellation of
  arbitrary scratch values and the three-invocation exchange. The paper's
  `h = 100` is not enumerable (1.77e20 wires), which is exactly why the result
  is a bit-complexity theorem.
* `transforms` implements the paper's DIF recursion literally: doubled
  butterflies, monomial twiddles by signed permutation, bit-reversed storage.
  Tests compare it with the naive transform and verify `F⁺F⁻ = I/t`.
* `resample` builds the Gaussian matrices `S`, `T`, the selector `C`, the
  diagonal `D`, computes `N = C T D`, checks `‖N−I‖` against the paper's bound,
  verifies the Harvey–van der Hoeven identity `T P_s F_s = P_t F_t S`, and
  composes `A = S/2`, `B_0 = D' J' C` into `P_s F_s = 2^{2α²} B0 P_t F_t A`.
  It uses `f64` and is clearly labelled as a numerical demonstration.
* `xor_net` implements the exact rational algebra behind the interchange
  count: full Bruhat factorization with `i128` fractions, verified by
  reconstructing `A = E1·Π·E2`, plus a modular shear executed via the
  factorization.

## Correctness strategy

* Every exact algorithm has a naive reference implementation in the same
  crate (`schoolbook`, `BigUint::div_rem`, `fft::dft_naive`,
  `transforms::naive_forward`, `crt::cyclic_convolution_naive`).
* Randomized tests compare fast paths against references across thresholds:
  NTT vs schoolbook at 20 000 bits, Newton division vs long division,
  Bluestein vs naive DFT, DIF vs naive synthetic transform.
* Final correction loops make the rounded results *self-verifying*: the
  reciprocal loop stops only when `2^p - c·x < c`, `isqrt` only when
  `s² ≤ a < (s+1)²`, and `fast_div_rem` verifies `q·b ≤ a < (q+1)·b`.

## Benchmarking

`benches/bench.rs` is a dependency-free harness (`harness = false`) that
multiplies pseudorandom operands and reports the best of a few repetitions.
It deliberately shows schoolbook/Karatsuba crossovers and the NTT's
`n log n`-like growth; it does not attempt to reproduce tape-machine costs.

## Scope boundaries

Deliberately not implemented: the fixed-tape cost accounting (counter
bookkeeping, head returns, descriptor stacks), the explicit `h = 100` motif
network, the full assembly with its rational parameter table, and the
transposition corollary. These are documented in `docs/paper-map.md` with
their formal statements. The repository is a tool for understanding and
executing the finite pieces, not a verification of the paper's asymptotics.
