import sys,numpy as np
from elp import setup, mp_rule
from elp2 import solve2, sample
n,k=int(sys.argv[1]),int(sys.argv[2]); trials=int(sys.argv[3])
M,D,orb,norb=setup(n,k); R0=mp_rule(M,D)
th=sample(2*n-3,3000)
rng=np.random.default_rng(0)
for trial in range(trials):
    vertex = trial%2==1
    w=rng.random(norb)
    if vertex: w=(w>0.5)*1.0+1e-12*0
    R=np.zeros_like(R0)
    for j in range(len(D)):
        W=np.nonzero(R0[j])[0]
        ws=np.array([w[orb[(j,t)]] for t in W])
        if ws.sum()==0: ws=np.ones(len(W))
        R[j,W]=ws/ws.sum()
    res,A=solve2(M,D,orb,norb,R,th,'all')
    print(trial,'vertex' if vertex else 'interior','gain',None if res.x is None else round(-res.fun,6))
