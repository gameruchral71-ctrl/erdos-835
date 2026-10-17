import sys, itertools, numpy as np, time
from scipy.optimize import linprog
from gen import *

def setup(n,k):
    M=Model(n); M.build_ps()
    D=M.datasets(k); didx={d:i for i,d in enumerate(D)}
    tidx={t:i for i,t in enumerate(M.T)}
    perms=list(itertools.permutations(range(n)))
    def mapc(s,c): return M.cidx[canon_split({s[i] for i in c},n)]
    cmap=[[mapc(s,c) for c in M.C] for s in perms]
    def mapt(s,t): return tidx[frozenset(canon_split({s[i] for i in S},n) for S in M.T[t])]
    tmap=[[mapt(s,t) for t in range(len(M.T))] for s in perms]
    orb={}; norb=0
    for xi,d in enumerate(D):
        for t in range(len(M.T)):
            if (xi,t) in orb: continue
            for si in range(len(perms)):
                d2=tuple(sorted(cmap[si][c] for c in d)); t2=tmap[si][t]
                orb[(didx[d2],t2)]=norb
            norb+=1
    return M,D,orb,norb

def mp_rule(M,D):
    R=np.zeros((len(D),len(M.T)))
    for j,d in enumerate(D):
        sc=[sum(M.ps(M.C[i],t) for i in d) for t in range(len(M.T))]
        m=min(sc); w=[t for t in range(len(M.T)) if sc[t]==m]
        R[j,w]=1/len(w)
    return R

def data_probs(M,D,th,t=0):
    cp=M.class_probs(t,th)
    return np.stack([multinom(d)*np.prod(cp[:,list(d)],axis=1) for d in D],1)

def solve(M,D,orb,norb,rule,th):
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
    obj=-A.sum(0)
    bounds=[(-r[o],1-r[o]) for o in range(norb)]
    res=linprog(obj,A_ub=-A,b_ub=np.zeros(len(A)),A_eq=Aeq,b_eq=np.zeros(len(Aeq)),bounds=bounds,method='highs')
    if res.x is None: print(res.message); return res,None
    return res,-res.fun

if __name__=='__main__':
    n,k,G=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
    t0=time.time()
    M,D,orb,norb=setup(n,k)
    print('n',n,'k',k,'data',len(D),'orbits',norb,'setup',round(time.time()-t0,1)); sys.stdout.flush()
    rng=np.random.default_rng(1)
    m=n+ (n-3)
    vals=[0.02,0.5,0.98]
    grid=np.array(list(itertools.product(vals,repeat=m)))
    th=np.vstack([grid, rng.random((G,m)), rng.random((G,m))**0.2, 1-rng.random((G,m))**0.2])
    rule=mp_rule(M,D)
    res,gain=solve(M,D,orb,norb,rule,th)
    print('status',res.status,'gain',gain,'time',round(time.time()-t0,1))
    np.save(f'sol_{n}_{k}.npy',res.x)
