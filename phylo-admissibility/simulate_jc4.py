# Monte Carlo check: simulate JC on 4-taxon trees, compare MP(uniform ties) with the new rule.
import numpy as np, itertools
from fractions import Fraction as Fr
from verify_jc4 import TREES, ps, mp, new_rule, canon
W=(Fr(1,3),Fr(2,9),Fr(1,18),Fr(1,18))
rng=np.random.default_rng(42)
def sim_site(ti,th):
    (a,b),(c,d)=TREES[ti]
    def ev(s,t):
        return s if rng.random()<0.25+0.75*t else rng.choice([x for x in range(4) if x!=s])
    u=rng.integers(4); v=ev(u,th[4]); x=[0]*4
    x[a]=ev(u,th[a]); x[b]=ev(u,th[b]); x[c]=ev(v,th[c]); x[d]=ev(v,th[d])
    return canon(tuple(x))
cache={}
def probs(d):
    key=tuple(sorted(d))
    if key not in cache:
        r=new_rule(key,W); m=mp(key)
        cache[key]=(m, r if r is not None else m)
    return cache[key]
for th in [np.array([0.9,0.9,0.3,0.3,0.5]),np.array([0.3,0.3,0.3,0.3,0.6]),np.array([0.5,0.2,0.5,0.2,0.4])]:
    for ti in range(3):
        N=200000; smp=0; snew=0
        for _ in range(N):
            d=(sim_site(ti,th),sim_site(ti,th))
            m,r=probs(d)
            smp+=float(m.get(ti,0)); snew+=float(r.get(ti,0))
        print('theta',th,'true',TREES[ti],'P(correct) MP=%.4f new=%.4f diff=%+.4f (se ~%.4f)'%(smp/N,snew/N,(snew-smp)/N,0.0015))
