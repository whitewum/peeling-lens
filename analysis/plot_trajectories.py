"""Degree-one trajectories against the fluid limit.

usage: python3 analysis/plot_trajectories.py TRAJ.txt OUT.png

Panels: (a) Q_n/M over the whole run; (b) the bottleneck window rescaled by
M^{1/4} in time and sqrt(M) in height; (c) an approximation to the Doob martingale of Q,
Q_n - Q_0 - A^Q_n; (d) the deviation Q_n - M ell(n/M). Red runs failed.
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trajio  # noqa: E402


def main(path, out):
    meta, runs, failed = trajio.load(path)
    M, k, d = meta["M"], meta["k"], meta["d"]
    fl = trajio.fluid(meta)
    ts, ell, sq, q4 = fl["tstar"], fl["ell"], np.sqrt(M), M ** 0.25
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    for a, fail in zip(runs, failed):
        n, Q, AQ = a[:, 0], a[:, 1], a[:, 4]
        t = n / M
        c = "tab:red" if fail else "tab:blue"
        ax[0, 0].plot(t, Q / M, c, lw=.5, alpha=.5)
        ax[0, 1].plot((t - ts) * q4, Q / sq, c, lw=.6, alpha=.6)
        ax[1, 0].plot(t, (Q - Q[0] - AQ) / sq, c, lw=.5, alpha=.6)
        ax[1, 1].plot(t, (Q - M * ell(t)) / sq, c, lw=.5, alpha=.6)
    tt = np.linspace(0, d / M, 2000)
    ax[0, 0].plot(tt, ell(tt), "k", lw=2, label="fluid limit ell(t)")
    ax[0, 0].set(xlabel="t = n/M", ylabel="Q_n / M", title="(a) degree-one count, whole run")
    ax[0, 0].legend()
    w = np.linspace(ts - 4 / q4, min(ts + 4 / q4, d / M), 400)
    ax[0, 1].plot((w - ts) * q4, M * ell(w) / sq, "k", lw=2, label="M ell / sqrt(M)")
    ax[0, 1].axhline(0, color="gray", lw=.8)
    ax[0, 1].set(xlim=(-4, 4), ylim=(-0.5, 12), xlabel="(t - t*) M^{1/4}", ylabel="Q_n / sqrt(M)",
                 title="(b) bottleneck window (red = failed)")
    ax[0, 1].legend()
    ax[1, 0].set(xlabel="t", ylabel="(Q_n - Q_0 - A^Q_n) / sqrt(M)", title="(c) approximately drift-corrected Q")
    ax[1, 1].set(xlabel="t", ylabel="(Q_n - M ell(n/M)) / sqrt(M)", title="(d) deviation from the fluid limit")
    for a_ in (ax[0, 0], ax[1, 0], ax[1, 1]):
        a_.axvline(ts, color="gray", ls=":")
    fig.suptitle(f"{meta['model']} model, k={k}, M={M}, d={d} (d/M={d / M:.4f})")
    fig.tight_layout()
    fig.savefig(out, dpi=110)

    mins = [a[np.abs(a[:, 0] / M - ts) < 4 / q4, 1].min() for a, f in zip(runs, failed) if not f]
    print(f"runs={len(runs)} P(fail)={failed.mean():.3f} t*={ts:.4f} Q0/M={runs[0][0, 1] / M:.4f}")
    if mins:
        print("successful runs, min Q in |t-t*|<4 M^-1/4: 10/50/90% =", np.percentile(mins, [10, 50, 90]))


if __name__ == "__main__":
    main(*sys.argv[1:3])
