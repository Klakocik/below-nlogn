# Contributing

Thanks for looking! This is an experimental, educational project about a
specific preprint, so the contribution rules are a little unusual.

## What belongs here

* Executable, tested models of the paper's finite constructions — the smaller
  and more self-contained, the better.
* Reference implementations and tests. Every fast path should be checked
  against a naive one.
* Corrections to `docs/paper-map.md` or the README, especially anything that
  makes the "what is actually implemented" boundary clearer or more honest.
* Benchmarks, as long as they do not overstate what they measure.

## What does not

* Claims that this crate is production-ready, faster than GMP/FLINT, or a
  faithful implementation of the full Turing-machine construction.
* New dependencies. The crate is deliberately zero-dependency; open an issue
  first if you believe something genuinely cannot be written in std.
* Reproducing the paper's fixed-tape cost accounting at scale — that needs a
  research program, not a PR. Documenting it is welcome.

## Ground rules

1. `cargo test` must pass and `cargo clippy --all-targets` should be clean.
2. Keep the paper correspondence explicit: cite section/lemma numbers in doc
   comments, and say whether code is exact or a numerical demonstration.
3. Do not copy text from the cited papers into the repository. Reference
   numbers and paraphrase.
4. Format with `cargo fmt`.

## Reporting issues

Include your platform, `rustc --version`, the command you ran, and what you
expected. For correctness bugs, a minimal failing input (hex literals are
fine) is gold.
