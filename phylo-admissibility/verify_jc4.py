#!/usr/bin/env python3
"""
Exact verification: under the Jukes-Cantor model, Maximum Parsimony with ties broken
uniformly at random is INADMISSIBLE for 4 taxa and k = 2 sites.

A rule M* (equal to MP except on four kinds of 2-site data, all of which are MP ties among
all three trees) satisfies  P(M* = T) - P(MP = T) >= 0  for every tree T and every choice of
the five branch lengths, with strict inequality on an open set. The difference is a polynomial
in th_e = exp(-4 t_e / 3) in [0,1]; nonnegativity is certified by exact rational Bernstein
coefficients on [0,1]^5. Leaves are numbered 0..3 here (1..4 in the write-up).

Run:  python3 verify_jc4.py      (needs sympy; a few seconds)
"""
import itertools, sympy as sp, sys
from fractions import Fraction as Fr
from math import comb, factorial
from collections import Counter
TH=sp.symbols('th0:5')
TREES=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
def tr(a,b,t): return (sp.Rational(1,4)+sp.Rational(3,4)*t) if a==b else (sp.Rational(1,4)-sp.Rational(1,4)*t)
def pprob(ti,x):
    (a,b),(c,d)=TREES[ti]; s=0
    for u,v in itertools.product(range(4),repeat=2):
        s+=sp.Rational(1,4)*tr(u,v,TH[4])*tr(u,x[a],TH[a])*tr(u,x[b],TH[b])*tr(v,x[c],TH[c])*tr(v,x[d],TH[d])
    return sp.expand(s)
def ps(ti,x):
    (a,b),(c,d)=TREES[ti]
    return min((u!=v)+(u!=x[a])+(u!=x[b])+(v!=x[c])+(v!=x[d]) for u,v in itertools.product(range(4),repeat=2))
def canon(x):
    mp={}; return tuple(mp.setdefault(v,len(mp)) for v in x)
CL=sorted(set(canon(x) for x in itertools.product(range(4),repeat=4)))
def cprob(ti,c):
    b=len(set(c)); return factorial(4)//factorial(4-b)*pprob(ti,c)
def blocks(c):
    B={}
    for i,v in enumerate(c): B.setdefault(v,set()).add(i)
    return sorted(map(frozenset,B.values()),key=lambda s:(len(s),sorted(s)))
def tree_with_cherry(pair):
    for i,(ab,cd) in enumerate(TREES):
        if set(ab)==set(pair) or set(cd)==set(pair): return i
def kind(c):
    b=blocks(c); sizes=sorted(len(s) for s in b)
    if sizes==[4]: return ('const',)
    if sizes==[1,3]: return ('single',next(iter(b[0])))
    if sizes==[1,1,2]: return ('paironly',frozenset(b[2]))   # e.g. 0012: pair {0,1} identical, others distinct
    if sizes==[2,2]: return ('split',frozenset(b[0]))
    return ('other',)
def new_rule(d, w):
    """d: tuple of 2 class patterns. returns dict tree->prob or None (= MP)."""
    k1,k2=kind(d[0]),kind(d[1])
    ks=sorted([k1,k2],key=str)
    W1,W2,W3,W4=w
    # type 1: constant site + 'paironly' site with identical pair P
    if {k1[0],k2[0]}=={'const','paironly'}:
        P=(k1 if k1[0]=='paironly' else k2)[1]; t=tree_with_cherry(P)
        return {i:(Fr(1,3)+W1 if i==t else Fr(1,3)-W1/2) for i in range(3)}
    # type 2: singleton sites at two different leaves i,j
    if k1[0]=='single' and k2[0]=='single' and k1[1]!=k2[1]:
        t=tree_with_cherry({k1[1],k2[1]})
        return {i:(Fr(1,3)-W2 if i==t else Fr(1,3)+W2/2) for i in range(3)}
    # types 3/4: two 'paironly' sites with identical pairs P,Q
    if k1[0]=='paironly' and k2[0]=='paironly' and k1[1]!=k2[1]:
        P,Q=k1[1],k2[1]
        if P&Q:   # type 3: pairs overlap (e.g. 0012 & 0102): favour the tree supported by neither
            tp,tq=tree_with_cherry(P),tree_with_cherry(Q); t=3-tp-tq
            return {i:(Fr(1,3)+W3 if i==t else Fr(1,3)-W3/2) for i in range(3)}
        else:     # type 4: complementary pairs (0012 & 0122): both support the same tree
            t=tree_with_cherry(P)
            return {i:(Fr(1,3)+W4 if i==t else Fr(1,3)-W4/2) for i in range(3)}
    return None
def mp(d):
    sc=[sum(ps(t,x) for x in d) for t in range(3)]; m=min(sc); W=[t for t in range(3) if sc[t]==m]
    return {t:Fr(1,len(W)) for t in W}
def bern(G):
    Pl=sp.Poly(G,*TH); degs=[max(Pl.degree(s),1) for s in TH]
    cur={mm:Fr(int(sp.numer(c)),int(sp.denom(c))) for mm,c in Pl.terms()}
    for vi,dg in enumerate(degs):
        new={}
        for mm,c in cur.items():
            i=mm[vi]
            for j in range(i,dg+1):
                m2=list(mm); m2[vi]=j; m2=tuple(m2); new[m2]=new.get(m2,0)+c*Fr(comb(j,i),comb(dg,i))
        cur=new
    full=1
    for dg in degs: full*=dg+1
    vals=list(cur.values())+([Fr(0)] if len(cur)<full else [])
    return min(vals),max(vals)
def check(w, verbose=True):
    data=list(itertools.combinations_with_replacement(range(len(CL)),2))
    ok=True
    for t in range(3):
        cp=[cprob(t,c) for c in CL]
        assert sp.expand(sum(cp)-1)==0
        G=0
        for d in data:
            dd=(CL[d[0]],CL[d[1]])
            r=new_rule(dd,w)
            if r is None: continue
            assert sum(r.values())==1 and min(r.values())>=0
            m=mp(dd)
            assert len(m)==3   # every changed data set is an MP tie among all three trees
            dl=r.get(t,0)-m.get(t,0)
            mult=2 if d[0]!=d[1] else 1
            if dl: G+=sp.Rational(dl.numerator,dl.denominator)*mult*cp[d[0]]*cp[d[1]]
        G=sp.expand(G)
        lo,hi=bern(G)
        if verbose: print('tree',TREES[t],'min bern',lo,'max',float(hi))
        ok&= lo>=0 and hi>0
    return ok
if __name__=='__main__':
    args=sys.argv[1:5] or ['1/3','2/9','1/18','1/18']
    w=tuple(Fr(x) for x in args)
    ok=check(w)
    print('weights',[str(x) for x in w])
    print('VERIFIED: M* >= MP for every tree and all branch lengths, strictly better somewhere' if ok else 'FAILED')
    sys.exit(0 if ok else 1)
