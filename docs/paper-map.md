# Paper map

A section-by-section map of *Integer multiplication below n log n*
(OpenAI, September 23, 2026) to this repository. "Implemented" means: runnable
code with tests. "Documented" means: the module docs explain the construction
and why a faithful reproduction is out of scope.

The headline result (Theorem 1.1):

```text
T_A(n) = O(n (lg n)^(1-κ)),   κ = 2^-182,
```

on one deterministic Turing machine with a fixed finite alphabet and a fixed
finite number of one-dimensional tapes, exact for every input length. The
ratio to `n lg n` is `(lg n)^(-κ)`.

## §1 Introduction

The reduction chain: radix-digit splitting → convolution → Gaussian-resampling
reduction to power-of-two multidimensional transforms → Nussbaumer–Quandalle
synthetic transforms over `C[y]/(y^r+1)` → address rearrangement + simultaneous
butterflies → final error `< 1/2` and carry propagation.

Executable pieces: everything below. The two headline tape procedures are:

1. **Address interchange** of two contiguous chunks of `u` bits in
   `O(V u^τ)`, `τ = 1 - 2^-50`, with intermediate XOR combinations and an
   exact final permutation (Section 4).
2. **Simultaneous normalized butterfly layers** `H0(u,v) = ((u+v)/2,(u-v)/2)`
   on selected coordinate bits in `O(V d^{λ'})`, `λ' = 1 - 2^-54`, with
   coefficient error `< √2 2^-p` (Section 5).

The recursive saving: a call at parameter `mf` makes `s < Wm` calls at
parameter `f` on volume `V/W`, so `s/W < m^τ` (or `m^σ`), and geometric sums
save a power of the parameter.

## §2 The machine model and elementary streams

Definitions: lexicographic arrays, implicit addresses, `Γ` self-delimiting
integer code, exact scan costs including counters and head returns,
controlled shifts and fixed rational scalings in linear time.

*Implemented analogue:* `xor_net::cyclic_row_classes` (Lemma 2.2 item 4),
`layers::packed_selected_bit_xor`, and the digit/stream plumbing throughout.
The full tape-scan cost bookkeeping (counters, head returns, work-tape reuse)
is documented, not reimplemented.

## §3 Two finite networks with a rank saving

Parameters: `h = 100`, `v = C(h,3) = 161700`, `N = v^3 = 4227952113000000`.
Wires: banks `X,Y`, side wires `A_{S,T}` for neighbor triples, centers `C_i`
(plus `C_*` in the complex case). The eight-row schedule subtracts the old
side/scatter contributions, copies/gathers, adds the new contributions, and
undoes the copy/gather, so arbitrary initial scratch values cancel exactly.
Net effect of a forward invocation: `y_S ← y_S + x_S` (XOR in the bit case).
Three invocations (forward, inverse-forward, forward) send `(X,Y) → (-Y, X)`;
in the bit case that is a genuine exchange.

Wire counts: `W_b = 177176569091445000000`,
`W_c = 1873807244643542670000`. Rank deficits:
`η_b = 339/22587335000000`, `η_c = 73/19906842167500`; both recurrences use
`τ = σ = 1 - 2^-50`.

*Implemented:* `motif.rs` instantiates the scalar network for small `h` and
verifies Lemma 3.2 by execution (`cargo run --release -- motif 6`). `h = 100`
is impossible to enumerate (1.77e20 wires), but the cancellation property is
independent of `h`.

## §4 Faster interchange of address chunks

* Lemma 4.3: every rational `m×m` matrix of rank `a` factors as
  `A = E1 · Π · E2`, with `E1, E2` invertible lower triangular and `Π` a
  partial permutation with exactly `a` ones.
* Lemma 4.4: the shear `(H,D) → (H + A·D, D)` uses exactly `a` field
  interchanges plus linear work; each one is the three-operation
  `D_j += H_i; (H_i,D_j) = (D_j,H_i); D_j = H_i - D_j`.
* Proposition 4.10: a circuit with `W` wires and edge-rank sum `s` interchanges
  two chunks of `e = m^k` digits in `O(V e^τ)` when `s/W < m^τ`.

*Implemented:* `xor_net.rs` — XOR swap, three-shear digit interchange, exact
rational Bruhat factorization (`Frac`/`Mat`, no floats), and
`shear_via_factorization` checked against the direct matrix–vector product.
*Documented:* the fixed-tape depth-first schedule and its `O(V)` per-node
cost.

## §5 Fast simultaneous synthetic butterfly layers

`H0` is applied coefficientwise on selected coordinate bits; packed
rectangles perform `y[j_i] ← y[j_i] XOR x[j_i]` on selected bits across
records. Cost `O(V d^{λ'})` with `λ' = 1 - 2^-54`; completed layers have
coefficient error `< √2 2^-p`.

*Implemented:* `layers.rs` — `h0_layer`, `h0_layer_inverse` (verifying
`H0² = I/2`), `packed_selected_bit_xor`, `xor_one_position`.
*Documented:* the exceptional-address repair and the tape batching.

## §6 Synthetic transforms and their tape layout

Ring `R_r = C[y]/(y^r+1)`; roots `ω_t = y^{2r/t}` act by signed coefficient
permutations; normalized `F_t^-`/`F_t^+`; decimation-in-frequency recursion
with `H0` and twiddle `ω_t^{-k}`; character cancellation
`Σ ω_t^{jk} = t·[t|j]`; stored frequency order is bit-reversed
(`B_t`); `F_t^+ F_t^- = t^{-1} I`.

*Implemented:* `transforms.rs` — exact ring monomial multiplication,
DIF forward/inverse, bit-reversal permutation, all tested against the naive
definition, including the `r=2, t=2` and `r=8, t=8` cases.

## §7 Resampling with the permutations left in the transform

Gaussian maps

```text
(Su)_k = α^{-1} Σ_j e^{-π α^{-2} (j - sk/t)^2} u_j
(Tu)_k =       Σ_j e^{-π α^2 (tj/s - k)^2} u_j
```

with periodic extension modulo `s`; permutations `(P_s u)_j = u_{tj}`,
`(P_t v)_k = v_{-sk}`; identity `T P_s F_s = P_t F_t S`. Selection
`C`, diagonal `D_j = e^{πα²β_j²}`, `N = C T D`, `||N - I|| < 2.01 e^{-πα²θ/2}`
with `θ = t/s - 1`; `A = S/2`, `B_0 = D' J' C`, `J' = N^{-1}/2`,
`D' = 2^{-(2α²-2)}D`, giving `P_s F_s = 2^{2α²} B_0 P_t F_t A`
(Lemma 7.2 for `α² < p` and `θ > p/α⁴`).

*Implemented:* `resample.rs` builds every matrix explicitly, verifies
Theorem 4.2 of Harvey–van der Hoeven (`harvey_van_der_hoeven_resampling_identity`),
the contraction estimate, and the composed identity numerically in `f64`.
*Documented:* the exact fixed-point contraction algorithms and their error
bounds (`16p`, `3p²/4`, `4`) from the cited manuscript, and their tape costs.

## §8 Exact integer multiplication and its cost

Fixed parameters: `b = ceil(log2 n)`, `p = 6b`, `d = floor(b^ε)`,
`K = floor(d^c)`, `T ∈ [4n/b, 8n/b)` a power of two split into `d` axes of
length `r` or `r/2` with `r` a power of two; `α = ceil((12d²b)^{1/4})`,
`γ = 2dα²`, `η = 1/(4d)`; constants `ε = 2^-75`, `c = 2^-56`,
`λ = 1 - 2^-52`, `λ' = 1 - 2^-54`, `C1 = 20`, `κ = 2^-182`.
Algorithm: digit streams → CRT layout `Φ` → two normalized synthetic
transforms → pointwise product/truncation → opposite transform → scale by
`S` and `2^{2b+4}S` → round → carry propagation → exact product.
Final error bound `< 1/2` after all truncations.

*Implemented:* `crt.rs` (Good–Thomas map and convolution preservation, tested
against flat cyclic convolution), `min`, `ntt`, `mul` (practical exact
multiplier), `ntt::convolve` (exact convolution), `fft.rs` (Bluestein for the
chirp transform).
*Documented:* the full assembly, its cost table and the `< 1/2` error proof.
The practical multiplier uses an exact number-theoretic transform instead of
the paper's approximate complex transforms, so there is no rounding step at
all.

## §9 Exact division and integer square root

Brent's Newton iterations with `B(p) = p (lg p)^{1-κ}`, doubling-precision
schedule with guard bits, exact remainder and square tests.

*Implemented:* `division.rs` — `reciprocal_floor` (integer Newton, exact
final test), `fast_div_rem` (64 guard bits + correction), `isqrt`
(`x ← (x + a/x)/2`, exact corrections), all cross-checked against reference
long division and `(s+1)² > a ≥ s²`.

## §10 A consequence for matrix transposition

Any multiplication machine with cost `M(m)` yields transposition at
`O(M(m) + m lg lg max(n1,n2))`, hence `O(n² (lg n)^{1-κ})` for an `n×n`
binary matrix. Pure corollary; no separate code.

## §11 Packed additions on grouped rectangles

The tape-level lemma behind simultaneous layers, with descriptors, grouped
rectangles, and one final repaired sort for exceptional addresses.

*Implemented:* the semantics (`packed_selected_bit_xor`); the descriptor
format, repair pass and cost formula `O(V(L^τ+1) + M A³ + ...)` are
documented only.
