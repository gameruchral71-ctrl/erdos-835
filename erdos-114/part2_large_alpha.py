"""
Part 2 of the degree-3 verification: alpha in [ALPHA1, 3*sqrt(3)/2].

Family (after the Eremenko-Hayman reduction):  p(z) = z^3 + a z + b,  a = alpha^(2/3),
b = sqrt(1 - 4 alpha^2/27); critical values exp(+-i t_c), sin t_c = alpha*sigma.
Rescaling z = sqrt(a) zeta, w = b + alpha*eta gives
    L(p) = (2/a) * int_0^pi h(eta(t)) dt,   eta(t) = (e^{it} - b)/alpha.
We prove   L_upper(alpha-box) < L0 = 2^(1/3) B(1/6, 1/2)   on every box.

To keep the critical parameter fixed under the alpha-box, t is written as
    t = t_c * v              (v in [0, 1 - d1])          piece A
    t = t_c + (pi - t_c) v   (v in [d2, 1])               piece B
and the two windows of t-length l1 = d1*t_c, l2 = d2*(pi - t_c) next to t_c are
bounded analytically (near_critical_constants).
"""
import sys, math, time
from flint import arb, acb
from rigorous_common import h_upper, h_float, near_critical_constants, SIGMA, PI

L0 = arb(2) ** (arb(1) / 3) * (arb(1) / 6).gamma() * arb(0.5).gamma() / (arb(1) / 6 + arb(0.5)).gamma()
ALPHA_MAX = 3 * arb(3).sqrt() / 2


def expm1i(x):
    """e^{ix} - 1 = 2i e^{ix/2} sin(x/2), tight for small x (x an arb)."""
    h = x / 2
    return 2 * acb(0, 1) * acb(h.cos(), h.sin()) * h.sin()


def integrate_piece(eta_of_v, jac_upper, v_lo, v_hi, tol, fl):
    """Upper bound of jac * int_{v_lo}^{v_hi} sup_alpha h(eta(t(v))) dv.
    fl(v) = float value of h at the box centre (only used to decide refinement)."""
    total = arb(0)
    stack = [(v_lo, v_hi, 0)]
    ncell = 0
    L = v_hi - v_lo
    while stack:
        lo, hi, depth = stack.pop()
        mid = 0.5 * (lo + hi)
        v = arb(mid, 0.5 * (hi - lo))
        eta = eta_of_v(v)
        U = h_upper(eta)
        if U is not None:
            slack = float(U.upper()) - fl(mid)
            if slack <= tol / L or depth >= 14:
                total += U.upper() * (arb(hi) - arb(lo))
                ncell += 1
                continue
        if depth > 40:
            raise RuntimeError('cannot enclose near v=%r' % mid)
        stack.append((lo, mid, depth + 1))
        stack.append((mid, hi, depth + 1))
    return jac_upper * total, ncell


def check_box(tc_lo, tc_hi, tol=None, ell=1e-2, verbose=False):
    """Box in the critical angle t_c in [tc_lo, tc_hi] (0 < t_c <= pi/2).
    alpha = sin(t_c)/sigma, b = cos(t_c): smooth in t_c up to t_c = pi/2."""
    tc = arb(0.5 * (tc_lo + tc_hi), 0.5 * (tc_hi - tc_lo))
    assert tc.lower() <= arb(tc_lo) and tc.upper() >= arb(tc_hi)
    tc = tc.min(PI / 2)
    al = tc.sin() / SIGMA
    b = tc.cos()
    al_lo, al_hi = float(al.lower()), float(al.upper())
    if tol is None:
        tol = 0.002 * al_lo ** (2 / 3)
    a_lo = al.lower() ** (arb(2) / 3)
    tc_mid = min(0.5 * (tc_lo + tc_hi), math.pi / 2)
    al_mid = math.sin(tc_mid) / float(SIGMA.mid())
    b_mid = math.cos(tc_mid)

    def hf(t):
        return h_float((complex(math.cos(t), math.sin(t)) - b_mid) / al_mid)

    # window lengths (in t) on each side of t_c
    l1 = ell * al_lo
    assert tc.lower() > 0
    d1 = l1 / float(tc.lower())
    d2 = l1 / (math.pi - float(tc.lower()))
    l1u = arb(d1) * tc.upper()
    l2u = arb(d2) * (PI - tc.lower())
    dmax = (l1u.max(l2u)) / al.lower()
    consts = near_critical_constants(dmax)
    if consts is None:
        return None
    K1, K3 = consts
    # |eta - i sigma| = 2 sin(|dt|/2)/alpha >= |dt| sinc(lmax/2)/alpha  for |dt| <= lmax
    lmax = l1u.max(l2u)
    root = (al.upper() / (lmax / 2).sinc()).sqrt()
    win = sum((2 * K1 * root * 2 * l.sqrt() + K3 * l) for l in (l1u, l2u))
    # eta(t) - i sigma = (e^{it} - e^{i t_c})/alpha = e^{i t_c} (e^{i (t - t_c)} - 1)/alpha
    ec = acb(tc.cos(), tc.sin())
    isg = acb(0, SIGMA)
    IA, nA = integrate_piece(lambda v: isg + ec * expm1i(tc * (v - 1)) / al, tc.upper(), 0.0, 1 - d1, tol,
                             lambda v: hf(tc_mid * v))
    IB, nB = integrate_piece(lambda v: isg + ec * expm1i((PI - tc) * v) / al, (PI - tc).upper(), d2, 1.0, tol,
                             lambda v: hf(tc_mid + (math.pi - tc_mid) * v))
    Ltot = (2 / a_lo) * (IA + IB + win)
    ok = Ltot.upper() < L0.lower()
    if verbose:
        print('t_c [%.6g,%.6g] alpha [%.6g,%.6g] L<=%.6f (L0=%.6f) window=%.4g cells=%d ok=%s' % (
            tc_lo, tc_hi, al_lo, al_hi, float(Ltot.upper()), float(L0.mid()), float((2 / a_lo * win).upper()), nA + nB, ok), flush=True)
    return ok, float(Ltot.upper())


if __name__ == '__main__':
    a1 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
    ratio = float(sys.argv[2]) if len(sys.argv) > 2 else 1.01
    # boxes in t_c from asin(a1*sigma) up to (a float >=) pi/2
    # optional explicit t_c range (for running chunks in parallel)
    t_start = float(sys.argv[3]) if len(sys.argv) > 3 else float((arb(a1) * SIGMA).asin().lower())
    t_end = float(sys.argv[4]) if len(sys.argv) > 4 else math.nextafter(math.pi / 2, 4.0)
    t0 = time.time()
    todo = []
    x = t_start
    while x < t_end:
        y = min(x * ratio, t_end)
        todo.append((x, y))
        x = y
    worst, nbox = -1.0, 0
    while todo:
        lo, hi = todo.pop(0)
        r = check_box(lo, hi, verbose=True)
        if r is None or not r[0]:
            print('  split', lo, hi, flush=True)
            m = 0.5 * (lo + hi)
            todo[0:0] = [(lo, m), (m, hi)]
            continue
        nbox += 1
        worst = max(worst, r[1])
    print('PART 2 DONE: t_c in [%r, %r] (alpha >= %g): %d boxes, every box has L <= upper bound < L0; '
          'max upper bound %.6f; L0 = %s; %.1fs' % (t_start, t_end, a1, nbox, worst, L0, time.time() - t0))
