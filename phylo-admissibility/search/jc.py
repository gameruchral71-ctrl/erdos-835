# Jukes-Cantor (4-state) version for 4 or 5 taxa, equivariant LP for dominating MP.
import sys, itertools, numpy as np
from math import factorial
from collections import Counter
from scipy.optimize import linprog
n=int(sys.argv[1]); k=int(sys.argv[2]); G=int(sys.argv[3]) if len(sys.argv)>3 else 3000
L=list(range(n))
# trees: n=4: (cherry ab, cherry cd) ; n=5: cherries ab, cd, middle e
def trees():
    T=[]
    if n==4:
        for ab in [(0,1),(0,2),(0,3)]:
            cd=tuple(sorted(set(L)-set(ab))); T.append((ab,cd,None))
    else:
        seen=set()
        for e in L:
            rest=[i for i in L if i!=e]
            for ab in itertools.combinations(rest,2):
                cd=tuple(sorted(set(rest)-set(ab)))
                key=frozenset([ab,cd])
                if key in seen: continue
                seen.add(key); T.append((ab,cd,e))
    return T
TREES=trees(); nT=len(TREES); m=2*n-3
# edges: pendant i -> th[i]; internal edges -> th[n], th[n+1] (n=5) or th[4] (n=4)
def pattern_prob(t, x, th):
    ab,cd,e=TREES[t]
    def tr(s1,s2,te): return np.where(s1==s2, 0.25+0.75*te, 0.25-0.25*te) if False else ((0.25+0.75*te) if s1==s2 else (0.25-0.25*te))
    tot=0
    if n==4:
        for u,v in itertools.product(range(4),repeat=2):
            tot=tot+0.25*tr(u,v,th[:,4])*tr(u,x[ab[0]],th[:,ab[0]])*tr(u,x[ab[1]],th[:,ab[1]])*tr(v,x[cd[0]],th[:,cd[0]])*tr(v,x[cd[1]],th[:,cd[1]])
    else:
        for w,u,v in itertools.product(range(4),repeat=3):
            tot=tot+0.25*tr(w,u,th[:,5])*tr(w,v,th[:,6])*tr(w,x[e],th[:,e])*tr(u,x[ab[0]],th[:,ab[0]])*tr(u,x[ab[1]],th[:,ab[1]])*tr(v,x[cd[0]],th[:,cd[0]])*tr(v,x[cd[1]],th[:,cd[1]])
    return tot
def pscore(t,x):
    ab,cd,e=TREES[t]; best=99
    if n==4:
        for u,v in itertools.product(range(4),repeat=2):
            best=min(best,(u!=v)+(u!=x[ab[0]])+(u!=x[ab[1]])+(v!=x[cd[0]])+(v!=x[cd[1]]))
    else:
        for w,u,v in itertools.product(range(4),repeat=3):
            best=min(best,(w!=u)+(w!=v)+(w!=x[e])+(u!=x[ab[0]])+(u!=x[ab[1]])+(v!=x[cd[0]])+(v!=x[cd[1]]))
    return best
# classes: set partitions in restricted growth form with <=4 blocks
CL=[]
for x in itertools.product(range(4),repeat=n):
    # restricted growth: x[0]=0 and each new label = max+1
    ok=True; mx=-1
    for v in x:
        if v>mx+1: ok=False;break
        mx=max(mx,v)
    if ok: CL.append(x)
nb=[len(set(x)) for x in CL]
labelings=[factorial(4)//factorial(4-b) for b in nb]
def canon(x):
    mp={}; out=[]
    for v in x:
        if v not in mp: mp[v]=len(mp)
        out.append(mp[v])
    return tuple(out)
cidx={c:i for i,c in enumerate(CL)}
def class_probs(t,th):
    return np.stack([labelings[i]*pattern_prob(t,CL[i],th) for i in range(len(CL))],1)
D=list(itertools.combinations_with_replacement(range(len(CL)),k))
didx={d:i for i,d in enumerate(D)}
def multinom(d):
    r=factorial(len(d))
    for v in Counter(d).values(): r//=factorial(v)
    return r
perms=list(itertools.permutations(L))
def tkey(t):
    ab,cd,e=TREES[t]; return frozenset([frozenset(ab),frozenset(cd)])
tid={tkey(t):t for t in range(nT)}
cmap=[[cidx[canon(tuple(CL[c][s.index(i)] for i in L))] for c in range(len(CL))] for s in perms]
# careful: permutation s maps leaf i -> s[i]; new pattern y[s[i]]=x[i]
cmap=[]
for s in perms:
    row=[]
    for c in range(len(CL)):
        y=[0]*n
        for i in L: y[s[i]]=CL[c][i]
        row.append(cidx[canon(tuple(y))])
    cmap.append(row)
tmap=[[tid[frozenset([frozenset(s[i] for i in TREES[t][0]),frozenset(s[i] for i in TREES[t][1])])] for t in range(nT)] for s in perms]
orb={}; norb=0
for xi,d in enumerate(D):
    for t in range(nT):
        if (xi,t) in orb: continue
        for si in range(len(perms)):
            orb[(didx[tuple(sorted(cmap[si][c] for c in d))],tmap[si][t])]=norb
        norb+=1
print('classes',len(CL),'data',len(D),'orbits',norb); sys.stdout.flush()
PS=np.array([[pscore(t,CL[c]) for t in range(nT)] for c in range(len(CL))])
R=np.zeros((len(D),nT))
for j,d in enumerate(D):
    sc=PS[list(d)].sum(0); w=np.nonzero(sc==sc.min())[0]; R[j,w]=1/len(w)
import os
if os.environ.get('ML') and os.path.exists(f'jcml_{n}_{k}.npy'):
    Lm=np.load(f'jcml_{n}_{k}.npy')
    for j in range(len(D)):
        Wt=np.nonzero(Lm[j]<=Lm[j].min()+1e-6)[0]; R[j]=0; R[j,Wt]=1/len(Wt)
elif os.environ.get('ML'):
    from scipy.optimize import minimize
    rep={}
    for (xi,t),o in orb.items():
        if o not in rep: rep[o]=(xi,t)
    val=np.zeros(norb)
    rngm=np.random.default_rng(5)
    for o,(xi,t) in rep.items():
        cnt=np.bincount(D[xi],minlength=len(CL)).astype(float)
        f=lambda z: -np.sum(cnt*np.log(np.maximum(class_probs(t,z[None,:])[0],1e-300)))
        best=np.inf
        for s0 in [np.full(m,0.5),np.full(m,0.99),np.full(m,0.01)]+[rngm.random(m) for _ in range(12)]:
            rr=minimize(f,s0,method='L-BFGS-B',bounds=[(0,1)]*m); best=min(best,rr.fun)
        val[o]=best
    Lm=np.zeros((len(D),nT))
    for (xi,t),o in orb.items(): Lm[xi,t]=val[o]
    np.save(f'jcml_{n}_{k}.npy',Lm)
    for j in range(len(D)):
        Wt=np.nonzero(Lm[j]<=Lm[j].min()+1e-6)[0]; R[j]=0; R[j,Wt]=1/len(Wt)
    print('ML rule built; data where ML != MP-uniform:', 'n/a')
if os.environ.get('TIE'):
    trng=np.random.default_rng(int(os.environ['TIE']))
    w=trng.random(norb)
    if int(os.environ['TIE'])%2: w=(w>0.5)*1.0
    for j in range(len(D)):
        Wt=np.nonzero(R[j])[0]; ws=np.array([w[orb[(j,t)]] for t in Wt])
        if ws.sum()==0: ws=np.ones(len(Wt))
        R[j]=0; R[j,Wt]=ws/ws.sum()
rng=np.random.default_rng(1)
grid=np.array(list(itertools.product([0.02,0.5,0.98],repeat=m)))
th=np.vstack([grid,rng.random((G,m)),rng.random((G,m))**0.2,1-rng.random((G,m))**0.2])
cp=class_probs(0,th)
print('prob sum check',np.abs(cp.sum(1)-1).max())
P=np.stack([multinom(d)*np.prod(cp[:,list(d)],axis=1) for d in D],1)
A=np.zeros((len(th),norb))
for xi in range(len(D)): A[:,orb[(xi,0)]]+=P[:,xi]
A=A/np.abs(A).max(1,keepdims=True)
r=np.zeros(norb)
for (xi,t),o in orb.items(): r[o]=R[xi,t]
Aeq=np.zeros((len(D),norb))
for xi in range(len(D)):
    for t in range(nT): Aeq[xi,orb[(xi,t)]]+=1
Aeq=np.unique(Aeq,axis=0)
import os
allowset=os.environ.get('ALLOW')
if allowset:
    want=[tuple(sorted(cidx[canon(tuple(int(ch) for ch in p))] for p in w.split(','))) for w in allowset.split()]
    perm_can=lambda d: min(tuple(sorted(cmap[si][c] for c in d)) for si in range(len(perms)))
    wantc=set(perm_can(w) for w in want)
    allowed=np.array([perm_can(d) in wantc for d in D])
else: allowed=np.ones(len(D),bool)
if os.environ.get('UNIQUE'): allowed=allowed & (R.max(1)==1)
if os.environ.get('GOOD'):
    good=np.zeros(len(D),bool)
    for j,d in enumerate(D):
        Wt=list(np.nonzero(R[j]>0)[0]); gap=min([Lm[j,t]-Lm[j].min() for t in range(nT) if t not in Wt],default=9)
        stab=[si for si in range(len(perms)) if tuple(sorted(cmap[si][c] for c in d))==tuple(d)]
        orbs=set(frozenset(tmap[si][t] for si in stab) for t in Wt)
        sym=len(orbs)==1 and all(o<=set(Wt) for o in orbs)
        good[j]=(len(Wt)==1 or sym) and gap>float(os.environ['GOOD'])
    print('good data',good.sum(),'of',len(D))
    np.save(f'jcgood_{n}_{k}.npy',good)
    allowed=allowed&good
al=np.zeros(norb,bool)
for (xi,t),o in orb.items(): al[o]|=allowed[xi]
bnd=[(-r[o],1-r[o]) if al[o] else (0,0) for o in range(norb)]
res=linprog(-A.sum(0),A_ub=-A,b_ub=np.zeros(len(A)),A_eq=Aeq,b_eq=np.zeros(len(Aeq)),bounds=bnd,method='highs')
print('JC n',n,'k',k,'status',res.status,'gain',-res.fun if res.x is not None else res.message)
np.save(f'jcsol_{n}_{k}.npy',res.x if res.x is not None else np.zeros(1))
# sparse solution
N=norb
tgt=-res.fun*0.25
Aub=np.vstack([np.hstack([-A,A]), np.hstack([-A.sum(0),A.sum(0)])[None,:]])
bub=np.concatenate([np.zeros(len(A)),[-tgt]])
bnds=[(0,b[1]) for b in bnd]+[(0,-b[0]) for b in bnd]
res2=linprog(np.ones(2*N),A_ub=Aub,b_ub=bub,A_eq=np.hstack([Aeq,-Aeq]),b_eq=np.zeros(len(Aeq)),bounds=bnds,method='highs')
xs=res2.x[:N]-res2.x[N:]
np.save(f'jcsparse_{n}_{k}.npy',xs)
seen=set()
def pname(c): return ''.join(map(str,CL[c]))
def tname(t): ab,cd,e=TREES[t]; return f"{ab[0]}{ab[1]}|{cd[0]}{cd[1]}"+(f" mid {e}" if e is not None else '')
for (xi,t),o in sorted(orb.items(),key=lambda z:z[1]):
    if o in seen: continue
    seen.add(o)
    if abs(xs[o])>1e-6*np.abs(xs).max():
        d=D[xi]
        Lg=Lm[xi]-Lm[xi].min() if os.environ.get('ML') else np.zeros(nT)
        print('orbit',o,'ties',int((R[xi]>0).sum()),'gaps',np.round(Lg,4),[pname(c) for c in d],'tree',tname(t),'PS',PS[list(d)].sum(0)[t],'minPS',PS[list(d)].sum(0).min(),'MP',round(R[xi,t],4),'delta',round(xs[o]/np.abs(xs).max(),4))
