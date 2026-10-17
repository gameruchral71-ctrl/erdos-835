import sys,numpy as np
from elp import setup, mp_rule
from elp2 import solve2, sample
n,k=int(sys.argv[1]),int(sys.argv[2]); tol=float(sys.argv[3]) if len(sys.argv)>3 else 1e-6
M,D,orb,norb=setup(n,k)
L=np.load(f'mlnegll_{n}_{k}.npy')
R=np.zeros_like(L)
for j in range(len(D)):
    W=np.nonzero(L[j]<=L[j].min()+tol)[0]; R[j,W]=1/len(W)
# check equivariance of ML rule (consistency of numerics)
vals={}
for (xi,t),o in orb.items(): vals.setdefault(o,set()).add(round(R[xi,t],6))
print('non-constant orbits',sum(len(v)>1 for v in vals.values()))
MP=mp_rule(M,D)
print('data where ML != MP:',int((np.abs(R-MP).max(1)>1e-9).sum()),'of',len(D))
th=sample(2*n-3,3000)
res,A=solve2(M,D,orb,norb,R,th,'all')
print(n,k,'ML gain',-res.fun if res.x is not None else res.message)
np.save(f'mlsol_{n}_{k}.npy',res.x)
