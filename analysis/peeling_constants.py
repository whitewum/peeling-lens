"""Analytic full-profile peeling constants; Python standard library only.

No simulation inputs or fitted parameters. Three variance computations:
closed form, elementary one-dimensional quadrature, and Poisson degree sums.
Run: python3 theory/peeling_constants.py
"""

import argparse
import math


def critical(k):
    lo, hi = 0.001, 10.0
    for _ in range(100):
        x = (lo + hi) / 2
        if math.expm1(x) > (k - 1) * x:
            hi = x
        else:
            lo = x
    x = (lo + hi) / 2
    q = -math.expm1(-x)
    a = q ** (k - 1)
    return x, q, a, x / a


def simpson(f, lo, hi, n):
    assert n > 0 and n % 2 == 0
    step = (hi - lo) / n
    total = f(lo) + f(hi)
    total += math.fsum((4 if i % 2 else 2) * f(lo + i * step)
                       for i in range(1, n))
    return total * step / 3


def constants(k, steps=2000):
    x, q, ac, mu = critical(k)
    E = math.exp(-x)
    s = 1 - ac
    initial = ((1 / mu + s * s) * math.exp(mu * (s * s - 1))
               - E * E * (1 + (1 - x) ** 2 / mu))
    total = E * q * (1 - k * ac / (k - 1))
    dynamic = total - initial
    C = 1 + (x - 1) * E

    def elementary(b):
        B = (1 + 2 * (x - 1) * E
             + ((1 - x + x * b) ** 2 + x * b) * math.exp(-2 * x + x * b))
        return (B - C * C * (b / ac) ** (1 / (k - 1))) / x

    # Independent degree-sum expression in the original time-change variable.
    def degree_sum(a):
        b = ac / a
        y = mu * a
        m = mu * a ** (k / (k - 1))
        pj = math.exp(-y)
        avg, avg2 = 0.0, 0.0
        for j in range(1, 100):
            pj *= y / j
            if j == 1:
                continue
            g = b * (1 - (1 - b) ** (j - 1)
                     + (j - 1) * b * (1 - b) ** (j - 2))
            prob = j * pj / m
            avg += prob * g
            avg2 += prob * g * g
        return mu * a ** (1 / (k - 1)) * (avg2 - avg * avg) / (mu * ac) ** 2

    integral = simpson(elementary, ac, 1, steps)
    discrete = simpson(degree_sum, ac, 1, steps)
    refined = simpson(elementary, ac, 1, steps * 2)
    err = max(abs(integral - dynamic), abs(discrete - dynamic), abs(refined - integral))
    assert initial > 0 and dynamic > 0
    assert err < 1e-9, (k, err)
    eta = E * ac * (x - 2) / math.sqrt(2 * total)
    return dict(k=k, x=x, mu=mu, alpha=mu / k, t_over_d=1 - q ** k,
                initial=initial, tau2=dynamic, sigma2=total,
                lam=initial / total, eta=eta, error=err)


def forward_dynamic_variance(k, steps=1000, max_degree=30):
    """Independent finite histogram covariance ODE, integrated by RK4.

V'=A V + V A^T + B, V(0)=0, on degrees 2,...,max_degree.
This is a numerical check of the coefficient, not a proof of the FCLT.
"""
    x, q, ac, mu = critical(k)
    end = mu / k * (1 - q ** k)
    degrees = list(range(2, max_degree + 1))
    dim = len(degrees)
    dt = end / steps

    def rhs(t, V):
        m = mu - k * t
        a = (m / mu) ** ((k - 1) / k)
        y = mu * a
        c = (k - 1) / m
        p = [0.0] * (max_degree + 2)
        z = math.exp(-y)
        for j in range(1, max_degree + 2):
            z *= y / j
            p[j] = j * z / m
        v = [p[j + 1] - p[j] for j in degrees]
        result = []
        for ii, i in enumerate(degrees):
            row = []
            for jj, j in enumerate(degrees):
                noise = -v[ii] * v[jj]
                if i == j:
                    noise += p[i] + p[i + 1]
                if abs(i - j) == 1:
                    noise -= p[max(i, j)]
                value = -c * (i + j) * V[ii][jj] + (k - 1) * noise
                if ii + 1 < dim:
                    value += c * (i + 1) * V[ii + 1][jj]
                if jj + 1 < dim:
                    value += c * (j + 1) * V[ii][jj + 1]
                row.append(value)
            result.append(row)
        return result

    def add(A, B, scale):
        return [[a + scale * b for a, b in zip(ar, br)] for ar, br in zip(A, B)]

    V = [[0.0] * dim for _ in degrees]
    for index in range(steps):
        t = index * dt
        A = rhs(t, V)
        B = rhs(t + dt / 2, add(V, A, dt / 2))
        C = rhs(t + dt / 2, add(V, B, dt / 2))
        D = rhs(t + dt, add(V, C, dt))
        V = [[V[i][j] + dt / 6 * (A[i][j] + 2 * B[i][j] + 2 * C[i][j] + D[i][j])
              for j in range(dim)] for i in range(dim)]
    return math.fsum(i * j * V[ii][jj] for ii, i in enumerate(degrees)
                     for jj, j in enumerate(degrees)) / (mu * ac) ** 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps', type=int, default=2000)
    parser.add_argument('--check-forward', action='store_true')
    args = parser.parse_args()
    fields = ['k', 'x', 'mu', 'alpha', 't_over_d', 'initial', 'tau2', 'sigma2', 'lam', 'eta', 'error']
    print(','.join(fields))
    for k in range(3, 8):
        r = constants(k, args.steps)
        print(','.join(str(r[f]) if f == 'k' else f'{r[f]:.12g}' for f in fields))
    if args.check_forward:
        for k in (3, 4):
            target = constants(k, args.steps)['tau2']
            for n in (500, 1000):
                value = forward_dynamic_variance(k, n)
                assert abs(value - target) < 1e-9
                print(f'# forward k={k} steps={n} tau2={value:.14g} error={value-target:.3g}')


if __name__ == '__main__':
    main()
