"""Shared loader for the output of sim/main.go."""
import numpy as np


def load(path):
    """Return (meta, runs, failed): runs[i] is an array with columns
    n, Q, U2, S, AQ, AV, U_0..U_40; failed[i] is True when peeling stopped early."""
    meta, rows, ends = {}, {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                meta = dict(kv.split("=") for kv in line[1:].split())
                continue
            p = line.split()
            if p[0] == "END":
                ends[int(p[1])] = int(p[2]) < int(p[3])
                continue
            rows.setdefault(int(p[0]), []).append([float(x) for x in p[1:]])
    for key in ("M", "k", "d"):
        meta[key] = int(meta[key])
    ids = sorted(rows)
    return meta, [np.array(rows[i]) for i in ids], np.array([ends[i] for i in ids])


def fluid(meta):
    """Fluid-limit quantities for the realised (M, k, d)."""
    from peeling_constants import critical
    M, k, d = meta["M"], meta["k"], meta["d"]
    x, q, ac, _ = critical(k)
    mu = k * d / M
    tstar = mu / k * (1 - q ** k)

    def a(t):
        return np.clip((mu - k * t) / mu, 1e-12, None) ** ((k - 1) / k)

    def ell(t):
        m = mu - k * t
        y = mu * np.clip(m / mu, 0, None) ** ((k - 1) / k)
        return m - y * (1 - np.exp(-y))

    return dict(x=x, q=q, ac=ac, mu=mu, tstar=tstar, a=a, ell=ell, mstar=mu - k * tstar)
