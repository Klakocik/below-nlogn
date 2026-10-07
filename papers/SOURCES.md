# Sources and attribution

This directory contains reference material for the repository. It is **not**
part of the `below-nlogn` library crate.

## OpenAI preprint (primary source)

* `paper.pdf`, `sources/` and the original `README.md` are taken verbatim from
  [`openai/math`](https://github.com/openai/math) —
  `preprints/Integer-multiplication-below-n-log-n-September-23-2026` —
  which is licensed under the Apache License 2.0, © OpenAI.
* Title: *Integer multiplication below n log n*. Author: OpenAI.
  Date: September 23, 2026.

## Secondary references consulted

* D. Harvey and J. van der Hoeven, *Integer Multiplication in Time
  O(n log n)*, Annals of Mathematics 193(2):563–617, 2021; result numbering
  in the OpenAI preprint follows the author-hosted 45-page manuscript at
  <https://www.texmacs.org/joris/nlogn/nlogn.pdf>. Only section and lemma
  numbers are referenced by the code; no text is copied here.
* A. Schönhage and V. Strassen, *Schnelle Multiplikation großer Zahlen*,
  Computing 7:281–292, 1971 (the `n log n` optimality conjecture that the
  preprint refutes in the fixed finite-tape model).
* H. J. Nussbaumer and P. Quandalle, *Fast computation of discrete Fourier
  transforms using polynomial transforms*, IEEE Trans. ASSP, 1979; and
  H. J. Nussbaumer, *Fast Fourier Transform and Convolution Algorithms*,
  Springer, 1981 (synthetic polynomial transforms used in Section 6).
* P. B. Brent, *Multiple-precision zero-finding methods and the complexity of
  elementary function evaluation*, 1976 (Newton reductions for division and
  square root, Section 9).

## Notes

* The practical Goldilocks-prime NTT in `src/ntt.rs` is standard
  engineering; the modulus `2^64 - 2^32 + 1` and generator `7` are folklore
  in the SNARK/DFT community and not taken from the preprint.
* All other code in this repository is original work of the
  `below-nlogn` contributors, licensed Apache-2.0 (see `../LICENSE`).
