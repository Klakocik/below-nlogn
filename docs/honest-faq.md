# Honest FAQ

**Is this library faster than GMP / FLINT?**
No. It is not faster than anything production-grade, and it was never going to
be. The paper's improvement is `(lg n)^(-κ)` with `κ = 2^-182`; a factor-2 win
needs `n ≈ 2^(2^(2^182))` bits. The bundled benchmarks show the expected
`n log n`-like behaviour of an exact NTT and nothing more.

**Then what is the point of the repository?**
Three things: (1) the paper is a landmark result and its finite gadgets are
beautiful and small enough to run — the Section 3 cancellation network, the
Section 6 synthetic transform, the Section 7 resampling identity, the
Section 4 exact shear factorization; (2) a self-contained zero-dependency
exact multiplier with division/sqrt is genuinely useful for teaching and
experimentation; (3) the asymptotics module makes the "why this can never be
practical" argument concrete and quantitative.

**Is the paper's algorithm implemented?**
No, and honestly: a faithful implementation would be a research program, not
a crate. The paper is a construction on a fixed multitape Turing machine with
astronomically large fixed networks (`W_b ≈ 1.77·10^20` wires) and constants.
This repository implements the runnable algebraic content: everything that
can be instantiated at small parameters and checked by execution. The rest is
mapped, quoted and explained in `docs/paper-map.md`.

**Why does `motif` only run small `h`?**
The network's ground set is all 3-subsets of `[100]`: `C(100,3) = 161700`
triples, three banks of `v^3 ≈ 4.2·10^15` values and `≈ 1.77·10^20` wires.
That is a *fixed* network, which is fine for a theorem but impossible to
enumerate. The cancellation property proved by Lemma 3.2 holds for every `h`;
the tests check it exactly for `h` up to 10.

**Why is there an NTT when the paper uses complex Fourier transforms?**
The paper's approximate transforms exist because its machine cannot do wide
modular arithmetic cheaply; a `u64` CPU can, and a Goldilocks-prime NTT gives
*exact* convolution with a small, checkable coefficient bound. The paper's
complex machinery is still implemented where it matters for the theory
(`fft`, `transforms`, `layers`, `resample`).

**Do the "verified" claims in the CLI multiply the input again?**
`below-nlogn mul` cross-checks the product with a second algorithm where
practical (not for schoolbook on huge inputs). Tests additionally compare all
three algorithms on random inputs.

**Are the floating-point resampling checks rigorous?**
No — they are `f64` numerical demonstrations of identities and bounds whose
rigorous versions are stated over fixed-point dyadic arithmetic in the paper.
The Rust code and the docs say so; do not cite this repository for the error
analysis.

**Can I use the code in production?**
It is Apache-2.0, so legally yes; technically, prefer GMP/FLINT/`num-bigint`
for production. This crate optimises for clarity, testability and the paper
correspondence, not throughput. It has not been audited, and `mul::mul` may
panic if a single-prime NTT cannot cover an input (it falls back only for
sizes, not for adversarial constants; the fallback path is exercised in
tests).

**What about division and square root?**
They are exact, self-correcting Newton implementations built on the shipped
multiplier — the paper's corollary in spirit and in interface
(`O(M(n))` in the multiplier's cost), verified against reference long
division.

**Do I need anything besides Rust?**
No. Zero dependencies; `cargo test`, `cargo run`, `cargo bench` are enough.
The paper sources in `papers/` are text/PDF only.

**How should I cite this?**
Cite the paper (see `README.md`/`CITATION.cff`). If you use the code, link the
repository and note that it is an experimental, partial implementation.
