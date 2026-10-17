#!/usr/bin/env python3
"""
Maximum likelihood, CFN model, 4 taxa, k = 4 sites.

Part A (exact): a rule M* that differs from ML only on three kinds of data sets satisfies
P(M* = T) >= P(ML = T) for every tree T and all branch lengths, strictly somewhere,
PROVIDED that on those three kinds of data ML ties all three trees (ML = 1/3 each).
Checked with exact rational Bernstein coefficients.

Part B (exact lower bound + numerical upper bound): for each of the three data kinds and
each tree, the likelihood reaches the same value V at an explicit boundary point (exact),
and 400 random-restart local optimisations never exceed V. Part B is NOT a proof that the
supremum equals V; that is the remaining gap for the ML conjecture.
"""
import itertools, sys
from fractions import Fraction as Fr
from math import comb
import sympy as sp
import numpy as np
np.seterr(all='ignore')
import warnings; warnings.filterwarnings('ignore')
from scipy.optimize import minimize

TH = sp.symbols('th0:5')            # th0..th3 pendant edges of leaves 0..3, th4 internal
TREES = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
def tree_with_cherry(p):
    for i, (ab, cd) in enumerate(TREES):
        if set(ab) == set(p) or set(cd) == set(p): return i

def pprob(ti, x, th):
    (a, b), (c, d) = TREES[ti]
    tr = lambda s1, s2, t: (1 + t) / 2 if s1 == s2 else (1 - t) / 2
    s = 0
    for u, v in itertools.product((0, 1), repeat=2):
        s += sp.Rational(1, 2) * tr(u, v, th[4]) * tr(u, x[a], th[a]) * tr(u, x[b], th[b]) * tr(v, x[c], th[c]) * tr(v, x[d], th[d])
    return s
CLASSES = [x for x in itertools.product((0, 1), repeat=4) if x[0] == 0]
def side(x):
    s = frozenset(i for i in range(4) if x[i]); o = frozenset(range(4)) - s
    return min(s, o, key=lambda z: (len(z), sorted(z)))
def cprob(ti, x, th):
    return pprob(ti, x, th) + pprob(ti, tuple(1 - v for v in x), th)

def partner(split, u):
    """the leaf on u's side of a 2|2 split (split given by one of its sides)"""
    S = split if u in split else frozenset(range(4)) - split
    return next(iter(S - {u}))

def new_delta(data):
    """data: tuple of 4 classes. Returns {tree: change relative to ML's uniform 1/3} or None."""
    S = [side(x) for x in data]
    const = [s for s in S if len(s) == 0]; sing = [next(iter(s)) for s in S if len(s) == 1]
    pair = [canon2(s) for s in S if len(s) == 2]
    # (i) one constant site, one singleton site (any leaf), two different 2|2 splits:
    #     +1/14 on each of the two trees displaying an observed split, -1/7 on the third tree
    if len(const) == 1 and len(sing) == 1 and len(pair) == 2 and pair[0] != pair[1]:
        t1, t2 = tree_with_cherry(pair[0]), tree_with_cherry(pair[1]); t3 = 3 - t1 - t2
        return {t1: Fr(1, 14), t2: Fr(1, 14), t3: Fr(-1, 7)}
    # (ii) singleton at u twice, singleton at v once, one 2|2 split separating u from v
    #      (u's partner w): -1/7 on cherry {u,v}, +1/7 on cherry {v,w}
    if len(sing) == 3 and len(pair) == 1:
        from collections import Counter
        c = Counter(sing)
        if sorted(c.values()) == [1, 2]:
            u = [i for i in c if c[i] == 2][0]; v = [i for i in c if c[i] == 1][0]
            w = partner(pair[0], u)
            if w != v:
                return {tree_with_cherry({u, v}): Fr(-1, 7), tree_with_cherry({v, w}): Fr(1, 7)}
    # (iii) singleton at u, split A twice, a different split B once; v = u's partner in A,
    #       w = u's partner in B: -1/7 on cherry {u,w}, +1/7 on cherry {v,w}
    if len(sing) == 1 and len(pair) == 3:
        from collections import Counter
        u = sing[0]; c = Counter(pair)
        if sorted(c.values()) == [1, 2]:
            A = [p for p in c if c[p] == 2][0]; B = [p for p in c if c[p] == 1][0]
            v = partner(A, u); w = partner(B, u)
            return {tree_with_cherry({u, w}): Fr(-1, 7), tree_with_cherry({v, w}): Fr(1, 7)}
    return None

def canon2(s):
    o = frozenset(range(4)) - s
    return min(s, o, key=sorted)

def bern(G):
    Pl = sp.Poly(G, *TH); degs = [max(Pl.degree(s), 1) for s in TH]
    cur = {m: Fr(int(sp.numer(c)), int(sp.denom(c))) for m, c in Pl.terms()}
    for vi, dg in enumerate(degs):
        new = {}
        for m, c in cur.items():
            i = m[vi]
            for j in range(i, dg + 1):
                mm = list(m); mm[vi] = j; mm = tuple(mm); new[mm] = new.get(mm, 0) + c * Fr(comb(j, i), comb(dg, i))
        cur = new
    full = 1
    for dg in degs: full *= dg + 1
    vals = list(cur.values()) + ([Fr(0)] if len(cur) < full else [])
    return min(vals), max(vals)

def main():
    from math import factorial
    from collections import Counter
    datasets = list(itertools.combinations_with_replacement(range(len(CLASSES)), 4))
    kinds = Counter()
    changed = []
    for d in datasets:
        dd = tuple(CLASSES[i] for i in d)
        r = new_delta(dd)
        if r is None: continue
        assert sum(r.values()) == 0
        mult = factorial(4)
        for v in Counter(d).values(): mult //= factorial(v)
        changed.append((d, mult, r))
    print("Part A: data sets where M* differs from ML:", len(changed))
    allok = True
    for t in range(3):
        cp = [sp.expand(cprob(t, x, TH)) for x in CLASSES]
        assert sp.expand(sum(cp) - 1) == 0
        G = 0
        for d, mult, r in changed:
            dl = r.get(t, 0)
            if dl:
                G += sp.Rational(dl.numerator, dl.denominator) * mult * cp[d[0]] * cp[d[1]] * cp[d[2]] * cp[d[3]]
        lo, hi = bern(sp.expand(G))
        print(f"  true tree {TREES[t]}: min Bernstein coefficient {lo}, max {float(hi):.5f}")
        allok &= lo >= 0 and hi > 0
    print("  Part A", "VERIFIED (conditional on the ties below)" if allok else "FAILED")

    print("Part B: likelihood values on the three data kinds (leaves 0..3; u=1, v=2, w=3)")
    reps = {'(i)  const, s1, {1,2}, {1,3}': [(0,0,0,0), (0,1,0,0), (0,1,1,0), (0,1,0,1)],
            '(ii) s1, s1, s2, {1,3}':        [(0,1,0,0), (0,1,0,0), (0,0,1,0), (0,1,0,1)],
            '(iii) s1, {1,2}, {1,2}, {1,3}': [(0,1,0,0), (0,1,1,0), (0,1,1,0), (0,1,0,1)]}
    rng = np.random.default_rng(0)
    for name, data in reps.items():
        assert new_delta(tuple(data)) is not None
        vals = []
        for t in range(3):
            Lsym = sp.Mul(*[cprob(t, x, TH) for x in data])
            f = sp.lambdify([TH], -sp.log(Lsym), 'numpy')
            best = None
            for _ in range(400):
                r = minimize(lambda z: f(z), rng.random(5), method='L-BFGS-B', bounds=[(1e-12, 1)] * 5)
                if best is None or r.fun < best.fun: best = r
            # exact value at a rounded maximiser (rational point in the closed box)
            pt = [Fr(v).limit_denominator(8) for v in best.x]
            exact = sp.nsimplify(Lsym.subs({TH[i]: sp.Rational(pt[i].numerator, pt[i].denominator) for i in range(5)}))
            vals.append(exact)
            print(f"  {name} tree {TREES[t]}: numerical sup {np.exp(-best.fun):.12f}; exact value {exact} at point {[str(p) for p in pt]}")
        print("   -> exact values equal for all three trees:", len(set(vals)) == 1)

if __name__ == "__main__":
    main()
