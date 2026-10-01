# Peeling Lens

**An Interactive Lab for Peeling Thresholds and Failure-Conditioned Measurements**

Peeling Lens is a companion to the [working paper](paper/main.pdf) on the critical window of
degree-one peeling (the decoder of invertible Bloom lookup tables, IBLTs). It
runs peeling on random k-uniform hypergraphs in the browser, so you can
watch how a run near the threshold succeeds or gets stuck, and compare it with
the closed-form limits.

**Version 0.1** · [Paper (PDF)](paper/main.pdf) · [LaTeX source](paper/main.tex)

Open `index.html` in a browser. It needs no build step and no server, so it
can also be served directly with GitHub Pages.

## What the lab shows

| Panel | Object | What to look for |
|---|---|---|
| A | Degree-one count Q_n / M over the whole run | Paths follow the fluid limit ℓ(t); at r = 0 it touches zero at t*, which defines the threshold. |
| B | Bottleneck, time × M^{1/4}, height ÷ √M | Near t* paths look like a parabola shifted as a whole; a downward shift reaches zero and the run fails. The shift is the score R, and failure ≈ {R ≤ 0}. |
| C | Approximately compensated V_s − V_0 − A^V_s of the weighted degree profile | Backward weights cancel the leading drift; the first-order drift correction A^V stays near zero. |
| D | Standard deviation of the corrected fluctuation and of R | Initial randomness plus accumulated peeling noise approach τ_k and σ_k, with σ_k² = initial variance + τ_k². |

Controls: mapping model (plain: k distinct cells per key; configuration:
iid endpoints), k = 3…7, M = 10³…10⁶, the load offset r in
d/M = α_k + r/√M, and the number of runs. All black curves come from
analytic formulas; nothing is fitted to simulation output.

## Repository layout

```
index.html                    interactive lab (single file, runs offline)
sim/main.go                   trajectory simulator: Q_n, degree histogram, first-order drift corrections
analysis/peeling_constants.py closed-form constants x_c, α_k, τ_k², σ_k², η_k (standard library only)
analysis/trajio.py            loader and fluid-limit formulas
analysis/plot_trajectories.py Q_n against the fluid limit, bottleneck zoom, drift-corrected Q
analysis/doob_weighted.py     Doob decomposition of the weighted functional V and the score R
figures/                      reference figures and their printed summaries
paper/                        English working paper, LaTeX source, PDF, and figure assets
```

## Reproduce the reference figures

Requires Go ≥ 1.21 and Python 3 with NumPy and Matplotlib. The trajectory file
is about 90 MB and takes a few minutes.

```sh
cd sim && go build -o ../peeltraj . && cd ..
./peeltraj -M 1000000 -runs 400 -stride 500 -seed 11 -model plain > traj.txt
python3 analysis/doob_weighted.py traj.txt figures/doob_weighted.png
python3 analysis/plot_trajectories.py traj.txt figures/trajectories.png
python3 analysis/peeling_constants.py
```

## Scope

The limits shown hold for fixed k ≥ 3 as M → ∞ in the window
d/M = α_k + r/√M. Finite-size failure rates can differ from the limit. A threshold-center
correction of order M^{-1/6} in the window coordinate is known in the
configuration setting; the same-order transfer to Plain is not established here. k = 2, convergence rates, and
the rare-failure regime far below the threshold are outside the scope.

## Working paper

**Min Wu. Failure-Conditioned Count Measurements at the Peeling Threshold.**
Working paper, version 0.1, October 2026. This is an early public manuscript,
not a peer-reviewed publication. It contains the main statements, proof
arguments, and numerical results, with explicit asymptotic and finite-size scope.

```sh
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Requires a standard TeX Live installation with `latexmk`, `lmodern`, and
`microtype`. All figure assets needed to compile the paper are included.
The paper's calibration figures summarize separate larger experiments; the
trajectory simulator here reproduces the lab's reference figures, not the
entire calibration study.

## Interpretation of the fluctuation panels

The simulator accumulates a first-order uniform-partner drift approximation.
It is not the exact finite-graph Doob compensator, particularly under Plain
mapping constraints. The corrected process illustrates the martingale
approximation used in the asymptotic analysis. Reference image labels mentioning
a Doob martingale should be read with this qualification.

Panel D stops once fewer than 99% of runs remain active. This limits, but does
not eliminate, survivor selection; it is not an unbiased estimator of the
terminal variance. Numerical agreement is illustrative, not a proof or a
finite-size guarantee.

## Status

The interactive lab supports Plain and configuration mappings. A sign-imbalance
and failure-conditioned measurement panel is planned. The English manuscript
and its figure assets are included in this release.

## License

Code is licensed under the MIT License. The manuscript and figure assets in
`paper/` and `figures/` are excluded from that license and remain copyrighted
by their respective authors. See [LICENSE](LICENSE).
