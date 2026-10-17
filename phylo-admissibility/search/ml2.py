import sys, numpy as np
from scipy.optimize import minimize
from multiprocessing import Pool
from gen import Model
from elp import setup
n,k=int(sys.argv[1]),int(sys.argv[2])
M,D,orb,norb=setup(n,k); nT=len(M.T); m=2*n-3
eps=1e-10
def negll(th,t,cnt):
    cp=M.class_probs(t,th[None,:])[0]
    return -np.sum(cnt*np.log(np.maximum(cp,1e-300)))
rep={}
for (xi,t),o in orb.items():
    if o not in rep: rep[o]=(xi,t)
def fit(o):
    j,t=rep[o]
    cnt=np.bincount(D[j],minlength=len(M.C)).astype(float)
    rng=np.random.default_rng(o)
    best=np.inf
    starts=[np.full(m,0.5),np.full(m,0.99),np.full(m,0.01)]+[rng.random(m) for _ in range(9)]
    for s in starts:
        r=minimize(negll,s,args=(t,cnt),method='L-BFGS-B',bounds=[(eps,1-eps)]*m)
        best=min(best,r.fun)
    return best
if __name__=='__main__':
    with Pool(4) as p: res=p.map(fit,range(norb))
    L=np.zeros((len(D),nT))
    for (xi,t),o in orb.items(): L[xi,t]=res[o]
    np.save(f'mlnegll_{n}_{k}.npy',L); print('done',norb)
