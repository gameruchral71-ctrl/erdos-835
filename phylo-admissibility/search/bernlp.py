# LP over Bernstein coefficients: find perturbation Delta of a base rule whose gain
# polynomial (true tree T0) has all Bernstein coefficients >= 0 on [0,1]^m.
import sys, itertools, numpy as np, sympy as sp, pickle
from math import comb
from fractions import Fraction as Fr
from scipy.optimize import linprog
from elp import setup
from gen import multinom

def sym_class_probs(M,t,th):
    E=M.edges(t); n=M.n; out=[]
    for c in M.C:
        s=0
        for A in M.even:
            term=1
            for e,S in enumerate(E):
                if len(A&S)%2==1: term*=th[e]
            s+=(-1)**len(A&c)*term
        out.append(sp.Poly(sp.Rational(2,2**n)*s,*th))
    return out

def to_bern_array(poly, th, degs):
    """exact Bernstein coefficients as dict multiindex->Fraction"""
    cur={m:Fr(int(sp.numer(c)),int(sp.denom(c))) for m,c in poly.terms()}
    for v,d in enumerate(degs):
        new={}
        for m,c in cur.items():
            i=m[v]
            for j in range(i,d+1):
                mm=list(m); mm[v]=j; mm=tuple(mm)
                new[mm]=new.get(mm,0)+c*Fr(comb(j,i),comb(d,i))
        cur=new
    return cur

def build(n,k,cache=None):
    M,D,orb,norb=setup(n,k); m=2*n-3
    th=sp.symbols(f't0:{m}')
    cp=sym_class_probs(M,0,th)
    # data prob polys for T0
    polys=[]
    for d in D:
        p=sp.Poly(multinom(d),*th)
        for c in d: p=p*cp[c]
        polys.append(p)
    degs=[k]*m
    idx={mi:i for i,mi in enumerate(itertools.product(*[range(d+1) for d in degs]))}
    Bx=[]
    for p in polys:
        b=to_bern_array(p,th,degs)
        Bx.append(b)
    return M,D,orb,norb,th,degs,idx,Bx,polys

def gain_matrix(D,orb,norb,idx,Bx):
    G=np.zeros((len(idx),norb))
    for xi in range(len(D)):
        o=orb[(xi,0)]
        for mi,c in Bx[xi].items(): G[idx[mi],o]+=float(c)
    return G
if __name__=='__main__':
    n,k=int(sys.argv[1]),int(sys.argv[2])
    out=build(n,k)
    pickle.dump(out[4:],open(f'bern_{n}_{k}.pkl','wb'))
    print('built')
