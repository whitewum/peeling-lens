"""Doob decomposition of the weighted functional used in the proof.

usage: python3 analysis/doob_weighted.py TRAJ.txt [OUT.png]

V_s = sum_j f_j(a_c / a(s/M)) U_j(s), f_j(b) = j b (1 - (1-b)^(j-1)). The
backward weights cancel the leading drift, so V_s - V_0 is close to a
martingale; the simulator also records a first-order drift correction A^V,
so the plotted approximate martingale is V_s - V_0 - A^V_s. The score at time s is
R(s) = (M m(t*) - V_s) / (c* sqrt M), c* = x_c, and failure ~ {R <= 0}.
Standard deviations are reported only while >= 99% of runs are still alive,
because failed runs leave the sample before t*.
"""
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trajio  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
J = np.arange(41.0)


def f(b):
    return np.where(J >= 2, J * b * (1 - (1 - b) ** (J - 1)), 0.0)


def theory(k):
    out = subprocess.run([sys.executable, os.path.join(HERE, "peeling_constants.py")],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    head = out[0].split(",")
    for line in out[1:]:
        row = dict(zip(head, line.split(",")))
        if int(row["k"]) == k:
            return {key: float(v) for key, v in row.items()}
    raise ValueError(k)


def main(path, out=None):
    meta, runs, failed = trajio.load(path)
    M, k = meta["M"], meta["k"]
    fl = trajio.fluid(meta)
    ts, ac, cst, sq = fl["tstar"], fl["ac"], fl["x"], np.sqrt(M)
    th = theory(k)

    def V(row):
        return f(ac / fl["a"](row[0] / M)) @ row[6:]

    grid = np.linspace(0, ts * 0.995, 140)
    mart = np.full((len(runs), len(grid)), np.nan)
    comp = np.full_like(mart, np.nan)
    score = np.full_like(mart, np.nan)
    for i, a in enumerate(runs):
        t = a[:, 0] / M
        keep = t < ts
        v = np.array([V(row) for row in a[keep]])
        alive_to = t[-1]
        g = grid <= alive_to
        mart[i, g] = np.interp(grid[g], t[keep], (v - v[0] - a[keep, 5]) / (cst * sq))
        comp[i, g] = np.interp(grid[g], t[keep], a[keep, 5] / (cst * sq))
        score[i, g] = np.interp(grid[g], t[keep], (M * fl["mstar"] - v) / (cst * sq))
    alive = np.isfinite(mart).mean(0)
    ok = alive >= 0.99
    print(f"{meta['model']} model, M={M}, k={k}, runs={len(runs)}, P(fail)={failed.mean():.3f}")
    print(f"theory: sd(initial)={th['initial'] ** .5:.4f}  tau_k={th['tau2'] ** .5:.4f}  sigma_k={th['sigma2'] ** .5:.4f}")
    print("    t   alive  mean(mart)  mean(A^V)  sd(mart)  sd(R)  corr(R,fail)")
    for gi in np.linspace(0, ok.sum() - 1, 7).astype(int):
        m, s = mart[:, gi], score[:, gi]
        fin = np.isfinite(s)
        c = np.corrcoef(s[fin], failed[fin])[0, 1] if failed[fin].std() > 0 else np.nan
        print(f"  {grid[gi]:.3f}  {alive[gi]:.2f}  {np.nanmean(m):+.4f}     {np.nanmean(comp[:, gi]):+.4f}    "
              f"{np.nanstd(m):.4f}   {np.nanstd(s):.4f}  {c:+.3f}")
    if not out:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for i in range(min(60, len(runs))):
        ax[0].plot(grid, mart[i], "tab:red" if failed[i] else "tab:blue", lw=.4, alpha=.5)
    ax[0].plot(grid[ok], np.nanmean(mart, 0)[ok], "k", lw=2, label="mean")
    ax[0].plot(grid[ok], np.nanmean(comp, 0)[ok], "k--", lw=1.2, label="mean drift correction A^V")
    ax[0].axhline(0, color="gray", lw=.6)
    ax[0].set(xlabel="t", title="(V_s - V_0 - A^V_s) / (c* sqrt M)")
    ax[0].legend(loc="upper left")
    for i in range(min(60, len(runs))):
        ax[1].plot(grid, score[i], "tab:red" if failed[i] else "tab:blue", lw=.4, alpha=.5)
    ax[1].axhline(0, color="gray", lw=.6)
    ax[1].set(xlabel="t", title="score R(s); failure ~ {R <= 0} near t*")
    ax[2].plot(grid[ok], np.nanstd(mart, 0)[ok], "tab:green", lw=2, label="sd of corrected fluctuation")
    ax[2].axhline(th["tau2"] ** .5, color="tab:green", ls="--", label=f"tau_k = {th['tau2'] ** .5:.3f} (at t*)")
    ax[2].plot(grid[ok], np.nanstd(score, 0)[ok], "k", lw=2, label="sd of R")
    ax[2].axhline(th["sigma2"] ** .5, color="k", ls="--", label=f"sigma_k = {th['sigma2'] ** .5:.3f} (at t*)")
    ax[2].axvline(ts, color="gray", ls=":")
    ax[2].set(xlabel="t", ylim=(0, None), title="accumulated noise (>= 99% of runs alive)")
    ax[2].legend(fontsize=8, loc="lower right")
    fig.suptitle(f"{meta['model']} model, k={k}, M={M}, {len(runs)} runs")
    fig.tight_layout()
    fig.savefig(out, dpi=105)


if __name__ == "__main__":
    main(*sys.argv[1:3])
