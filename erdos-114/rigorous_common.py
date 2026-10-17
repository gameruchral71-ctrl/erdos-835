"""
Rigorous helpers (ball arithmetic via python-flint / Arb) for the degree-3 case
of Erdos problem #114.

Model cubic: q(zeta) = zeta^3 + zeta.  Critical points c = +-i/sqrt(3),
critical values +-i*sigma with sigma = 2/(3*sqrt(3)).

h(eta) = sum over the three roots zeta_k of q(zeta) = eta of 1/|q'(zeta_k)|.

Everything here returns *upper bounds* for h (the only direction the proof needs).
"""
import cmath, math
import numpy as np
from flint import arb, acb, ctx

ctx.prec = 80

SQRT3 = arb(3).sqrt()
SIGMA = 2 / (3 * SQRT3)            # critical value modulus
C_ABS = 1 / SQRT3                   # |c|
PI = arb.pi()


def _box(z0, r):
    return acb(arb(z0.real, r), arb(z0.imag, r))


def _inside(K, D):
    return D.real.contains_interior(K.real) and D.imag.contains_interior(K.imag)


def krawczyk(z0, eta, r):
    """Return the box D = z0 +- r (each coordinate) if Krawczyk proves that for every
    eta in the ball `eta`, q(zeta) = eta has exactly one root in D. Otherwise None."""
    z0a = acb(z0.real, z0.imag)
    D = _box(z0, arb(r))
    qp0 = 3 * z0a * z0a + 1
    Y = 1 / qp0
    K = z0a - Y * (z0a ** 3 + z0a - eta) + (1 - Y * (3 * D * D + 1)) * (D - z0a)
    if _inside(K, D):
        return D
    return None


def h_upper(eta):
    """Upper bound (arb, use .upper()) of h over the complex ball `eta`, or None if
    the three roots could not be isolated (caller should subdivide)."""
    em = complex(float(eta.real.mid()), float(eta.imag.mid()))
    erad = float(abs(eta.real.rad()) + abs(eta.imag.rad()))
    roots = np.roots([1, 0, 1, -em])
    boxes = []
    for z0 in roots:
        z0 = complex(z0)
        # one Newton step in floating point to polish the centre
        for _ in range(2):
            z0 = z0 - (z0 ** 3 + z0 - em) / (3 * z0 * z0 + 1)
        d = abs(3 * z0 * z0 + 1)
        base = 2.0 * erad / d + 1e-14 * (1 + abs(z0))
        D = None
        for mult in (1.0, 2.0, 4.0, 16.0):
            D = krawczyk(z0, eta, base * mult)
            if D is not None:
                break
        if D is None:
            return None
        boxes.append(D)
    for i in range(3):
        for j in range(i + 1, 3):
            if boxes[i].overlaps(boxes[j]):
                return None
    total = arb(0)
    for D in boxes:
        m = (3 * D * D + 1).abs_lower()
        if not m > 0:
            return None
        total += 1 / m
    return total


def h_float(e):
    rs = np.roots([1, 0, 1, -e])
    return float(sum(1 / abs(3 * r * r + 1) for r in rs))


def near_critical_constants(dmax):
    """For |eta - i*sigma| <= dmax (an upper bound, arb), return (K1, K3) such that
       h(eta) <= 2*K1/sqrt|eta - i sigma| + K3,
    via the factorisation q(zeta) - i sigma = (zeta - c)^2 (zeta + 2c), c = i/sqrt3,
    and Rouche on the disks |zeta - c| <= rho (two roots) and |zeta + 2c| <= rho2 (one root).
    Returns None if dmax is too large."""
    three_c = SQRT3
    dm = float(dmax.upper())
    rho = arb(min(0.25, 1.2 * math.sqrt(dm / 1.5) + 1e-12))
    rho2 = arb(min(0.25, 1.2 * dm / 2.0 + 1e-12))
    # Rouche: on |z-c| = rho, |(z-c)^2 (z+2c)| >= rho^2 (3|c| - rho) > |d|  => exactly 2 roots inside
    if not (rho * rho * (three_c - rho) > dmax):
        return None
    # on |z+2c| = rho2, |(z-c)^2 (z+2c)| >= (3|c| - rho2)^2 rho2 > |d|  => exactly 1 root inside
    if not ((three_c - rho2) ** 2 * rho2 > dmax):
        return None
    if not (rho + rho2 < three_c):
        return None
    # two roots near c: 1/|q'| = 1/(3|z-c||z+c|),  |z-c| = sqrt(|d|/|z+2c|) >= sqrt(|d|/(3|c|+rho)),
    #                   |z+c| >= 2|c| - rho
    K1 = (three_c + rho).sqrt() / (3 * (2 * C_ABS - rho))
    # third root near -2c:  |z-c| >= 3|c|-rho2,  |z+c| >= |c|-rho2
    K3 = 1 / (3 * (three_c - rho2) * (C_ABS - rho2))
    return K1, K3
