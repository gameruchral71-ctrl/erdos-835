# General-n CFN, equivariant LP for dominating a given rule (MP / ML).
import itertools, numpy as np
from math import factorial
from collections import Counter

def canon_split(S, n):
    S=frozenset(S)
    return S if 0 not in S else frozenset(range(n))-S

def trees(n):
    """unrooted binary trees as frozensets of nontrivial splits (canonical side without leaf 0)."""
    # build by stepwise addition; represent tree as edge list on nodes
    res=set()
    def splits_of(edges, n):
        adj={}
        for a,b in edges: adj.setdefault(a,[]).append(b); adj.setdefault(b,[]).append(a)
        out=[]
        for a,b in edges:
            # leaves reachable from b without crossing a
            seen={a,b}; st=[b]; L=set()
            while st:
                u=st.pop()
                if isinstance(u,int) and u<n: L.add(u)
                for w in adj[u]:
                    if w not in seen: seen.add(w); st.append(w)
            out.append(canon_split(L,n))
        return out
    def rec(edges, nxt, cnt):
        if nxt==n:
            sp=frozenset(s for s in splits_of(edges,n) if 1<len(s)<n-1)
            res.add(sp); return
        for i,(a,b) in enumerate(edges):
            v=('v',cnt)
            ne=edges[:i]+edges[i+1:]+[(a,v),(v,b),(v,nxt)]
            rec(ne,nxt+1,cnt+1)
    rec([(0,('v',0)),(1,('v',0)),(2,('v',0))],3,1)
    return sorted(res, key=lambda t: sorted(sorted(s) for s in t))

class Model:
    def __init__(self,n):
        self.n=n
        self.T=trees(n)
        # site classes: subsets not containing 0 (pattern = indicator); 2^(n-1)
        self.C=[frozenset(c) for r in range(n) for c in itertools.combinations(range(1,n),r)]
        self.cidx={c:i for i,c in enumerate(self.C)}
        self.even=[frozenset(A) for r in range(0,n+1,2) for A in itertools.combinations(range(n),r)]
    def edges(self,t):
        # pendant edges as singletons, then internal splits
        return [canon_split({i},self.n) for i in range(self.n)]+sorted(self.T[t],key=sorted)
    def class_probs(self,t,th):
        E=self.edges(t); n=self.n
        # path set: edge with split S included iff |A∩S| odd
        terms=[]
        for A in self.even:
            mask=[len(A&S)%2==1 for S in E]
            term=np.ones(th.shape[0])
            for e,m in enumerate(mask):
                if m: term=term*th[:,e]
            terms.append((A,term))
        out=np.zeros((th.shape[0],len(self.C)))
        for i,c in enumerate(self.C):
            s=0
            for A,term in terms:
                s=s+(-1)**len(A&c)*term
            out[:,i]=2*s/2**n
        return out
    def ps(self,c,t):
        # parsimony score of binary character c (subset) on tree t
        # Fitch on tree: use splits: score = min number of edges to cut so that each part is monochromatic
        # brute: for small n, compute via Fitch over an explicit tree build
        return self._ps_cache[(c,t)]
    def build_ps(self):
        self._ps_cache={}
        for t in range(len(self.T)):
            E=self.edges(t)
            m=len(E)
            for c in self.C:
                best=None
                # choose edge subset F; a character is consistent with cutting F if every component monochromatic
                # equivalently c's induced split must be generated: c = XOR-combination of chosen splits' sides
                for r in range(m+1):
                    for F in itertools.combinations(range(m),r):
                        x=frozenset()
                        for e in F: x=x^E[e]
                        if canon_split(x,self.n)==canon_split(c,self.n) or (len(c)==0 and len(x)==0):
                            best=r; break
                    if best is not None: break
                self._ps_cache[(c,t)]=best
    def datasets(self,k):
        return list(itertools.combinations_with_replacement(range(len(self.C)),k))

def multinom(d):
    m=factorial(len(d))
    for v in Counter(d).values(): m//=factorial(v)
    return m
