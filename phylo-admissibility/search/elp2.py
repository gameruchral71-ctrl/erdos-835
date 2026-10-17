import sys, itertools, numpy as np, time
from scipy.optimize import linprog
from elp import setup, mp_rule, data_probs
from gen import *
def solve2(M,D,orb,norb,rule,th,mode,objw=None):
    nT=len(M.T)
    P=data_probs(M,D,th,0)
    A=np.zeros((len(th),norb))
    for xi in range(len(D)): A[:,orb[(xi,0)]]+=P[:,xi]
    A=A/np.abs(A).max(1,keepdims=True)
    r=np.zeros(norb)
    for (xi,t),o in orb.items(): r[o]=rule[xi,t]
    Aeq=np.zeros((len(D),norb))
    for xi in range(len(D)):
        for t in range(nT): Aeq[xi,orb[(xi,t)]]+=1
    Aeq=np.unique(Aeq,axis=0)
    obj=-A.sum(0) if objw is None else objw
    if mode=='unique':  # change only data with a unique MP tree
        uniq=np.zeros(norb,bool)
        for (xi,t),o in orb.items():
            if rule[xi].max()==1: uniq[o]=True
        bounds=[(-r[o],1-r[o]) if uniq[o] else (0,0) for o in range(norb)]
    elif mode=='ties':   # only re-break ties: keep zero where MP is zero
        bounds=[(-r[o],1-r[o]) if r[o]>0 else (0,0) for o in range(norb)]
    else:
        bounds=[(-r[o],1-r[o]) for o in range(norb)]
    res=linprog(obj,A_ub=-A,b_ub=np.zeros(len(A)),A_eq=Aeq,b_eq=np.zeros(len(Aeq)),bounds=bounds,method='highs')
    return res,A
def sample(m,G,seed=1):
    rng=np.random.default_rng(seed)
    grid=np.array(list(itertools.product([0.02,0.5,0.98],repeat=m)))
    return np.vstack([grid, rng.random((G,m)), rng.random((G,m))**0.2, 1-rng.random((G,m))**0.2])
if __name__=='__main__':
    n,k,mode=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
    M,D,orb,norb=setup(n,k); rule=mp_rule(M,D)
    th=sample(2*n-3,3000)
    res,A=solve2(M,D,orb,norb,rule,th,mode)
    print(n,k,mode,'status',res.status,'gain',-res.fun)
    np.save(f'sol_{n}_{k}.npy',res.x)
