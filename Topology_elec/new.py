import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import matplotlib.pyplot as plt

PI = np.pi

# ============================================================
# integrand
# ============================================================

def f(k2, k1, m, n, alpha):
    num = 1.0 - np.cos(m * k1) * np.cos(n * k2)
    den = alpha * np.sin(k1 / 2.0) ** 2 - np.sin(k2 / 2.0) ** 2
    return num / den


# ============================================================
# singular locations
# ============================================================

def singular_points(k1, alpha):

    rhs = np.sqrt(alpha) * abs(np.sin(k1 / 2.0))

    if rhs >= 1.0:
        return []

    ks = 2.0 * np.arcsin(rhs)

    pts = []

    if 0.0 < ks < PI:
        pts.extend([ks, -ks])

    return sorted(pts)


# ============================================================
# PV inner integral
# ============================================================

def inner_pv(k1, m, n, alpha, delta=1e-6):

    pts = singular_points(k1, alpha)

    if len(pts) == 0:
        val, _ = quad(
            lambda k2: f(k2, k1, m, n, alpha),
            -PI,
            PI,
            limit=400
        )
        return val

    intervals = []

    a = -PI

    for s in pts:

        if s - delta > a:
            intervals.append((a, s - delta))

        a = s + delta

    if a < PI:
        intervals.append((a, PI))

    total = 0.0

    for left, right in intervals:

        val, _ = quad(
            lambda k2: f(k2, k1, m, n, alpha),
            left,
            right,
            limit=400
        )

        total += val

    return total


# ============================================================
# full 2D PV integral
# ============================================================

def I(alpha, m, n):

    val, err = quad(
        lambda k1: inner_pv(k1, m, n, alpha),
        -PI,
        PI,
        limit=300
    )

    return val


# ============================================================
# scan alpha
# ============================================================

def scan_alpha(m, n,
               a_min=0.05,
               a_max=2.0,
               N=80):

    alphas = np.linspace(a_min, a_max, N)

    vals = []

    for a in alphas:

        try:
            v = I(a, m, n)

        except Exception:

            v = np.nan

        vals.append(v)

        print(
            f"alpha={a:.6f}  I={v:.8e}"
        )

    return alphas, np.array(vals)


# ============================================================
# root search
# ============================================================

def find_roots(m, n,
               a_min=0.05,
               a_max=2.0,
               N=80):

    alphas, vals = scan_alpha(
        m,
        n,
        a_min,
        a_max,
        N
    )

    roots = []

    for i in range(N - 1):

        f1 = vals[i]
        f2 = vals[i + 1]

        if np.isnan(f1) or np.isnan(f2):
            continue

        if f1 * f2 > 0:
            continue

        left = alphas[i]
        right = alphas[i + 1]

        try:

            root = brentq(
                lambda a: I(a, m, n),
                left,
                right,
                xtol=1e-6,
                rtol=1e-6,
                maxiter=100
            )

            roots.append(root)

            print(
                f"root found: {root:.10f}"
            )

        except Exception:
            pass

    return roots, alphas, vals


# ============================================================
# main
# ============================================================

if __name__ == "__main__":

    m = 2
    n = 1

    roots, alphas, vals = find_roots(
        m,
        n,
        a_min=0.05,
        a_max=2.0,
        N=100
    )

    print()
    print("================================")
    print(f"(m,n)=({m},{n})")
    print("roots:")
    print(roots)
    print("================================")

    plt.figure(figsize=(8,5))
    plt.plot(alphas, vals)
    plt.axhline(0.0)
    plt.xlabel("alpha")
    plt.ylabel("I(alpha)")
    plt.title(f"I(alpha), (m,n)=({m},{n})")
    plt.grid(True)
    plt.show()