"""
Part 1 of the degree-3 verification: alpha in [0, ALPHA1] (the cusp at z^3 - 1).

With a = alpha^(2/3), b = sqrt(1 - 4 alpha^2/27), eps = (1-b)/alpha = (4/27) alpha/(1+b),
eta0(tau) = i tau/(1 - i alpha tau/2) (Cayley parametrisation of the circle, t = 2 arctan(alpha tau/2)),
w(tau) = 1/(1 + alpha^2 tau^2/4) = |d eta0/d tau|, one has the exact identity
    G(alpha) := (L0 - L(p_alpha))/sqrt(a)
              = 2 * int_0^oo [ |eta0|^(-2/3) - h(eta0 + eps) ] w dtau ,
where h(eta) = sum_{zeta^3+zeta=eta} 1/|3 zeta^2 + 1|. (G(0) = 2D = 0.93830708...)
We prove G(alpha) > 0 on boxes covering [0, ALPHA1], i.e. L(p_alpha) < L0.

Split  int_0^oo = int_0^{T(alpha)} + tail,  T(alpha) = tau_c(alpha) + S,
tau_c = 2 sigma/(1 + sqrt(1 - alpha^2 sigma^2)) the parameter of the critical value i*sigma.
"""
import sys, math, time
from flint import arb, acb
from rigorous_common import h_upper, h_float, near_critical_constants, SIGMA

I = acb(0, 1)


def tail_constant(R):
    """C with |h(eta) - |eta|^(-2/3)| <= C |eta|^(-2) for all |eta| >= R (Lemma T)."""
    U = arb(R) ** (arb(1) / 3)
    s0 = arb(1) / 60
    # Rouche condition at U (RHS*U^3 is decreasing in U, so it holds for all larger U)
    lhs = 3 * s0
    rhs = s0 / U ** 2 + 3 * s0 ** 2 / U ** 6 + (arb(1) / 3 + s0 / U ** 4) ** 3
    assert lhs > rhs, 'Rouche fails'
    assert arb(3).sqrt() * U - 2 / (3 * U) > 2 * s0 / U ** 5     # disks disjoint
    s = s0 / U ** 5
    Y = 2 * s / U + (1 / (3 * U) + s) ** 2 / U ** 2
    v = 1 / (3 * U ** 2) + Y
    g = 2 * v + v * v
    assert g < 1
    E = v * v / 2 + arb(3) / 8 * g * g / (1 - g) ** (arb(5) / 2)
    return (Y + E) * U ** 4


def integrate(eta_w_of_v, v_lo, v_hi, tol, fl, maxdepth=14):
    """Upper bound of int_{v_lo}^{v_hi} sup h(eta(v)) * w(v) dv; eta_w_of_v(v) -> (eta, w_upper, jac_upper)."""
    total = arb(0)
    stack = [(v_lo, v_hi, 0)]
    L = v_hi - v_lo
    ncell = 0
    while stack:
        lo, hi, depth = stack.pop()
        mid = 0.5 * (lo + hi)
        v = arb(mid, 0.5 * (hi - lo))
        eta, wu, ju = eta_w_of_v(v)
        U = h_upper(eta)
        if U is not None:
            val = U.upper() * wu * ju
            if float(val.upper()) - fl(mid) <= tol / L or depth >= maxdepth:
                total += val * (arb(hi) - arb(lo))
                ncell += 1
                continue
        if depth > 40:
            raise RuntimeError('cannot enclose at v=%r' % mid)
        stack.append((lo, mid, depth + 1))
        stack.append((mid, hi, depth + 1))
    return total, ncell


def check_box(al_lo, al_hi, S=16.0, ell=0.01, tol=0.01, verbose=False):
    al = arb(0.5 * (al_lo + al_hi), 0.5 * (al_hi - al_lo))
    assert al.lower() <= arb(al_lo) and al.upper() >= arb(al_hi)
    if al_lo == 0.0:
        al = al.nonnegative_part()
    b = (1 - 4 * al * al / 27).sqrt()
    eps = (arb(4) / 27) * al / (1 + b)
    tauc = 2 * SIGMA / (1 + (1 - al * al * SIGMA * SIGMA).sqrt())
    almid = 0.5 * (al_lo + al_hi)
    bmid = math.sqrt(1 - 4 * almid ** 2 / 27)
    epsmid = (4 / 27) * almid / (1 + bmid)
    sg = float(SIGMA.mid())
    tcmid = 2 * sg / (1 + math.sqrt(1 - almid ** 2 * sg ** 2))
    D = 1 - I * al * tauc / 2

    def eta_at(off, tau):
        # eta = i sigma + i (tau - tau_c) / ((1 - i al tau/2)(1 - i al tau_c/2))
        return acb(0, SIGMA) + I * off / ((1 - I * al * tau / 2) * D)

    def w_up(tau):
        x = (al * tau).nonnegative_part().lower()   # 0 <= x <= alpha*tau
        return 1 / (1 + x ** 2 / 4)

    def flo(tau):
        e0 = 1j * tau / (1 - 1j * almid * tau / 2)
        return h_float(e0 + epsmid) / (1 + almid ** 2 * tau ** 2 / 4)

    # piece A: tau = tau_c * v, v in [0, 1 - ell/tau_c_lo]
    vA = 1 - ell / float(tauc.lower())

    def fA(v):
        tau = tauc * v
        return eta_at(tauc * (v - 1), tau), w_up(tau).upper(), tauc.upper()
    IA, nA = integrate(fA, 0.0, vA, tol, lambda v: flo(tcmid * v) * tcmid)
    # window lengths actually excluded: left = tauc*(1-vA) in [ell, ell*tauc_hi/tauc_lo], right = ell
    lL = tauc.upper() * (1 - arb(vA))
    # piece B: tau = tau_c + s, s in [ell, S]
    def fB(s):
        tau = tauc + s
        return eta_at(s, tau), w_up(tau).upper(), arb(1)
    IB, nB = integrate(fB, ell, S, tol, lambda s: flo(tcmid + s))
    # windows: |eta - i sigma| >= |tau - tau_c| / |(1 - i al tau/2)(1 - i al tau_c/2)| >= |tau-tau_c|/(1+al^2 tmax^2/4)
    tmax = tauc.upper() + ell
    fac = (1 + al.upper() ** 2 * tmax ** 2 / 4)
    dmax = lL.max(arb(ell))
    consts = near_critical_constants(dmax)
    assert consts is not None
    K1, K3 = consts
    win = sum(2 * K1 * fac.sqrt() * 2 * l.sqrt() + K3 * l for l in (lL, arb(ell)))
    Bup = IA + IB + win
    # term A: int_0^{T_lo} tau^(-2/3) (1+al^2 tau^2/4)^(-2/3) dtau, lower bound, T_lo = tauc_lo + S
    Tlo = tauc.lower() + S
    Alow = arb(0)
    edges = [arb(0)] + [Tlo * arb(2) ** (-k) for k in range(40, -1, -1)]
    for x0, x1 in zip(edges[:-1], edges[1:]):
        Alow += 3 * (x1 ** (arb(1) / 3) - x0 ** (arb(1) / 3)) * (1 + al.upper() ** 2 * x1 ** 2 / 4) ** (-arb(2) / 3)
    # tail over tau >= T(alpha) >= Tlo
    R0 = Tlo / (1 + al.upper() ** 2 * Tlo ** 2 / 4).sqrt()       # |eta0| >= R0 there
    epsu = eps.upper()
    R = R0 - epsu
    C = tail_constant(float(R.lower()) * (1 - 1e-9))   # valid for |eta| >= a value below R
    k = R0 / R
    tail = C * k ** 2 / Tlo + epsu * k ** (arb(5) / 3) * Tlo ** (-arb(2) / 3)
    Glow = 2 * (Alow.lower() - Bup.upper() - tail.upper())
    ok = Glow > 0
    if verbose:
        print('alpha [%.6g,%.6g]: G >= %.5f  (A>=%.4f, h-part<=%.4f [window %.4f], tail<=%.2e, C=%.3f) cells=%d ok=%s' % (
            al_lo, al_hi, float(Glow.lower()), float(Alow.lower()), float(Bup.upper()), float(win.upper()),
            float(tail.upper()), float(C.upper()), nA + nB, ok), flush=True)
    return ok, float(Glow.lower())


if __name__ == '__main__':
    a1 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
    nb = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    t0 = time.time()
    edges = [a1 * i / nb for i in range(nb + 1)]
    todo = list(zip(edges[:-1], edges[1:]))
    worst, n = 1e9, 0
    while todo:
        lo, hi = todo.pop(0)
        r = check_box(lo, hi, verbose=True)
        if not r[0]:
            m = 0.5 * (lo + hi)
            todo[0:0] = [(lo, m), (m, hi)]
            print('  split', lo, hi, flush=True)
            continue
        n += 1
        worst = min(worst, r[1])
    print('PART 1 DONE: alpha in [0, %g]: %d boxes, G(alpha) >= %.5f > 0 on all; %.1fs' % (a1, n, worst, time.time() - t0))
