#!/usr/bin/env python3
"""
Exact computer verification that Maximum Parsimony (ties broken uniformly at random)
is an INADMISSIBLE tree-reconstruction method for 5 taxa and k = 3 sites under the
symmetric 2-state (CFN / Neyman) model, answering one case of Mike Steel's
"phylogenetic inadmissibility" conjecture.

Everything is computed independently of the search code: trees are enumerated by
their two cherries, site-pattern probabilities are summed over internal states,
parsimony scores come from brute force over internal states, and nonnegativity of
the gain polynomial is certified by exact (rational) Bernstein coefficients on [0,1]^7.

Run:  python3 verify.py        (needs sympy; takes about a minute)
"""
import itertools, sys
from fractions import Fraction as Fr
from math import comb, factorial
from collections import Counter
import sympy as sp

LEAVES = range(5)
L = frozenset(LEAVES)

# ---- the 15 unrooted binary trees on 5 leaves: two disjoint cherries + middle leaf
TREES = []
for e in LEAVES:
    rest = sorted(L - {e})
    for ab in itertools.combinations(rest, 2):
        cd = tuple(sorted(set(rest) - set(ab)))
        t = frozenset([frozenset(ab), frozenset(cd)])
        if t not in TREES:
            TREES.append(t)
assert len(TREES) == 15
def middle(t): return next(iter(L - set().union(*t)))
def tree(ch1, ch2):
    return TREES.index(frozenset([frozenset(ch1), frozenset(ch2)]))

# ---- CFN model. theta = 1 - 2p (p = change probability on an edge), in (0,1).
TH = sp.symbols('th0:7')   # th0..th4 pendant edges of leaves 0..4; th5, th6 internal
def leaf_pattern_prob(t, x):
    """P(leaf states = x) on tree t, edges: th[i] pendant of leaf i,
    th5 internal edge of the first cherry (sorted), th6 of the second."""
    t = TREES[t]
    (ab, cd) = sorted(t, key=sorted)
    a, b = sorted(ab); c, d = sorted(cd); m = middle(t)
    def tr(s1, s2, th):
        return (1 + th) / 2 if s1 == s2 else (1 - th) / 2
    tot = 0
    for w, u, v in itertools.product((0, 1), repeat=3):
        tot += sp.Rational(1, 2) * tr(w, u, TH[5]) * tr(w, v, TH[6]) * tr(w, x[m], TH[m]) \
             * tr(u, x[a], TH[a]) * tr(u, x[b], TH[b]) * tr(v, x[c], TH[c]) * tr(v, x[d], TH[d])
    return sp.expand(tot)

# site classes = binary patterns up to complementing all states (CFN is symmetric)
CLASSES = []
for x in itertools.product((0, 1), repeat=5):
    if x[0] == 0:
        CLASSES.append(x)
def ones_small(x):
    """the smaller side of the bipartition induced by a pattern"""
    s = frozenset(i for i in LEAVES if x[i] == 1)
    return min(s, L - s, key=lambda z: (len(z), sorted(z)))

def class_prob(t, x):
    xb = tuple(1 - v for v in x)
    return leaf_pattern_prob(t, x) + leaf_pattern_prob(t, xb)

def pscore(t, x):
    t = TREES[t]
    (ab, cd) = sorted(t, key=sorted)
    a, b = sorted(ab); c, d = sorted(cd); m = middle(t)
    best = 99
    for w, u, v in itertools.product((0, 1), repeat=3):
        s = (w != u) + (w != v) + (w != x[m]) + (u != x[a]) + (u != x[b]) + (v != x[c]) + (v != x[d])
        best = min(best, s)
    return best

def mp_rule(data):
    sc = [sum(pscore(t, x) for x in data) for t in range(15)]
    m = min(sc); W = [t for t in range(15) if sc[t] == m]
    return {t: Fr(1, len(W)) for t in W}

# ---- the dominating rule: equal to MP except on four kinds of 3-site data sets
def new_rule(data):
    S = [ones_small(x) for x in data]
    const = [s for s in S if len(s) == 0]
    sing = [s for s in S if len(s) == 1]
    pair = [s for s in S if len(s) == 2]
    # A: one constant site + two pair-splits {a,b},{a,c} sharing one leaf.
    #    MP ties 6 trees; choose ab|de or ac|de (de = the two untouched leaves), 1/2 each.
    if len(const) == 1 and len(pair) == 2 and len(pair[0] & pair[1]) == 1:
        P, Q = pair; de = L - (P | Q)
        return {tree(P, de): Fr(1, 2), tree(Q, de): Fr(1, 2)}
    # B: singleton sites at two different leaves u, v and one pair-split P with u in P, v not in P.
    #    MP ties the 3 trees with cherry P; choose the two whose other cherry contains v.
    if len(sing) == 2 and len(pair) == 1 and sing[0] != sing[1]:
        P = pair[0]; U = sing[0] | sing[1]
        if len(U & P) == 1:
            v = next(iter(U - P)); y, z = sorted(L - P - {v})
            return {tree(P, {v, y}): Fr(1, 2), tree(P, {v, z}): Fr(1, 2)}
    # C: one singleton site at leaf u + two pair-splits {a,b},{a,c} not containing u; e = 5th leaf.
    #    MP ties 6 trees; put 1/3 on ab|ue and ac|ue, 1/6 on ab|ce and ac|be, 0 on ab|uc, ac|ub.
    if len(sing) == 1 and len(pair) == 2 and len(pair[0] & pair[1]) == 1:
        u = next(iter(sing[0]))
        if u not in pair[0] | pair[1]:
            P, Q = pair; e = next(iter(L - (P | Q) - {u}))
            out = {}
            for X, Y in ((P, Q), (Q, P)):
                c = next(iter(Y - X))
                out[tree(X, {u, e})] = Fr(1, 3)
                out[tree(X, {e, c})] = Fr(1, 6)
            return out
    # D: three distinct pair-splits forming a path q1-q2-r2-r1 ... i.e. Q, R disjoint and P
    #    meeting each in one leaf; e = the leaf in none. MP's unique tree is Q|R.
    #    Keep Q|R with prob 2/3; put 1/6 on Q|{e, P∩R} and 1/6 on R|{e, P∩Q}.
    if len(pair) == 3 and len(set(pair)) == 3:
        for P in pair:
            Q, R = [p for p in pair if p != P]
            if not (Q & R) and len(P & Q) == 1 and len(P & R) == 1:
                e = next(iter(L - (P | Q | R)))
                return {tree(Q, R): Fr(2, 3),
                        tree(Q, {e} | (P & R)): Fr(1, 6),
                        tree(R, {e} | (P & Q)): Fr(1, 6)}
    return None

def bernstein_min(G):
    Pl = sp.Poly(G, *TH)
    degs = [max(Pl.degree(s), 1) for s in TH]
    cur = {m: Fr(int(sp.numer(c)), int(sp.denom(c))) for m, c in Pl.terms()}
    for vi, dg in enumerate(degs):
        new = {}
        for m, c in cur.items():
            i = m[vi]
            for j in range(i, dg + 1):
                mm = list(m); mm[vi] = j
                mm = tuple(mm)
                new[mm] = new.get(mm, 0) + c * Fr(comb(j, i), comb(dg, i))
        cur = new
    # include zero coefficients for missing multi-indices
    full = 1
    for dg in degs: full *= dg + 1
    vals = list(cur.values()) + ([Fr(0)] if len(cur) < full else [])
    return min(vals), max(vals)

def main():
    datasets = list(itertools.combinations_with_replacement(range(len(CLASSES)), 3))
    changed = []
    for d in datasets:
        data = [CLASSES[i] for i in d]
        r = new_rule(data)
        if r is None: continue
        assert sum(r.values()) == 1 and all(v >= 0 for v in r.values())
        mult = factorial(3)
        for v in Counter(d).values(): mult //= factorial(v)
        changed.append((d, mult, r, mp_rule(data)))
    print("data sets (unordered, up to complement) where the rule differs from MP:", len(changed))
    probs = {}
    allok = True
    pt = {s: sp.Rational(1, 3) + sp.Rational(i, 17) for i, s in enumerate(TH)}  # interior test point
    for t in range(15):
        cp = [class_prob(t, x) for x in CLASSES]
        assert sp.expand(sum(cp) - 1) == 0
        G = 0
        for d, mult, r, mp in changed:
            dl = r.get(t, 0) - mp.get(t, 0)
            if dl:
                G += sp.Rational(dl.numerator, dl.denominator) * mult * cp[d[0]] * cp[d[1]] * cp[d[2]]
        G = sp.expand(G)
        lo, hi = bernstein_min(G)
        val = G.subs(pt)
        ok = lo >= 0 and val > 0
        allok &= ok
        print(f"true tree {sorted(map(sorted, TREES[t]))} middle {middle(TREES[t])}: "
              f"min Bernstein coef = {lo}, max = {float(hi):.5f}, gain at interior point = {float(val):.3e}  {'OK' if ok else 'FAIL'}")
    print("VERIFIED: new rule >= MP for every tree and every branch length, strictly better somewhere"
          if allok else "VERIFICATION FAILED")
    return 0 if allok else 1

if __name__ == "__main__":
    sys.exit(main())
