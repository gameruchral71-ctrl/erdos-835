import sys,pickle,numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from elp import setup, mp_rule
from bernlp import gain_matrix
n,k,base=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
M,D,orb,norb=setup(n,k)
th,degs,idx,Bx,polys=pickle.load(open(f'bern_{n}_{k}.pkl','rb'))
G=gain_matrix(D,orb,norb,idx,Bx)
if base=='mp': R=mp_rule(M,D); allowed=np.ones(len(D),bool)
else:
    R=np.load(f'mlrule_{n}_{k}.npy'); allowed=np.array([R[j].min()>0.3 for j in range(len(D))]) if n==4 else np.ones(len(D),bool)
nT=len(M.T)
r=np.zeros(norb); allow=np.zeros(norb,bool)
for (xi,t),o in orb.items(): r[o]=R[xi,t]; allow[o]|=allowed[xi]
Aeq=np.zeros((len(D),norb))
for xi in range(len(D)):
    for t in range(nT): Aeq[xi,orb[(xi,t)]]+=1
Aeq=np.unique(Aeq,axis=0)
scale=np.abs(G).max(1,keepdims=True); keep=scale[:,0]>0
Gs=G[keep]/scale[keep]
bounds=[(-r[o],1-r[o]) if allow[o] else (0,0) for o in range(norb)]
res=linprog(-Gs.sum(0),A_ub=-Gs,b_ub=np.zeros(len(Gs)),A_eq=Aeq,b_eq=np.zeros(len(Aeq)),bounds=bounds,method='highs')
print('LP status',res.status,'objective',-res.fun)
x=res.x
np.save(f'bsol_{base}_{n}_{k}.npy',x)
# sparse: minimize L1 subject to Bernstein >=0 and total >= 1/4 of max
N=norb
tgt=-res.fun*0.25
Aub=np.vstack([np.hstack([-Gs,Gs]), np.hstack([-Gs.sum(0),Gs.sum(0)])[None,:]])
bub=np.concatenate([np.zeros(len(Gs)),[-tgt]])
bnds=[(0,max(0,b[1])) for b in bounds]+[(0,max(0,-b[0])) for b in bounds]
res2=linprog(np.ones(2*N),A_ub=Aub,b_ub=bub,A_eq=np.hstack([Aeq,-Aeq]),b_eq=np.zeros(len(Aeq)),bounds=bnds,method='highs')
xs=res2.x[:N]-res2.x[N:]
print('sparse status',res2.status,'support',int((np.abs(xs)>1e-9).sum()))
np.save(f'bsparse_{base}_{n}_{k}.npy',xs)
# rationalize & exact check
from itertools import product
def exact_check(xr):
    tot={}
    for xi in range(len(D)):
        o=orb[(xi,0)]
        if xr[o]==0: continue
        for mi,c in Bx[xi].items(): tot[mi]=tot.get(mi,0)+c*xr[o]
    vals=list(tot.values())
    return min(vals) if vals else None, max(vals) if vals else None
for den in (6,12,24,60,120,840,5040):
    xr=[Fr(v).limit_denominator(den) for v in xs]
    # fix row sums exactly: check
    ok=all(sum(xr[orb[(xi,t)]] for t in range(nT))==0 for xi in range(len(D)))
    okb=all(-r[o]-1e-12<=float(xr[o])<=1-r[o]+1e-12 for o in range(N))
    lo,hi=exact_check(xr)
    print(den,'rowsum ok',ok,'bounds ok',okb,'min bern',lo,'max',float(hi) if hi is not None else None)
    if ok and okb and lo is not None and lo>=0 and hi>0:
        pickle.dump(xr,open(f'brat_{base}_{n}_{k}.pkl','wb')); print('FOUND exact certificate with denominator',den); break
