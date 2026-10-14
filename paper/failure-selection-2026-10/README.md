# Decoding Failure Biases IBLT Count Measurements: A Gaussian-Selection Law at the Peeling Threshold

Min Wu, October 2026.

English version of the working paper on failure-conditioned count measurements in the standard IBLT:
the joint Gaussian limit of the count energy and the precritical score in the critical window, the
resulting Gaussian-selection law after decoding failure, the sign reversal of its bias between k = 4
and k = 5, the upper confidence bound for the unknown difference size, and its extension along the
load axis.

- [main.pdf](main.pdf): full paper with appendices; source `main.tex`, `sections/`, `appendix/`, `references.bib`.
- [appendix-macro.pdf](appendix-macro.pdf): supplementary appendix with the key lemmas for the
  macroscopic conditional law; source `appendix/appendix-macro.tex`.
- `figures/`: figures used by the paper.

Build: `latexmk -pdf main.tex` in this directory; `latexmk -pdf appendix-macro.tex` in `appendix/`.

This version is posted to make the current results citable. The proofs are drafts and have not yet
been independently refereed. The paper and figures are not covered by the repository's MIT software
license; copyright remains with the author.
