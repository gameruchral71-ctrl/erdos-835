"""
Computer checks for the note on Erdos Problem 835 (chromatic number of J(2k, k)).

1. Clique lemma: for every (k-1)-set L of [2k], the k-sets containing L form a clique
   of size k+1 in J(2k, k). (Checked for k = 2..6.)
2. k = 2: J(4,2) has a proper 3-colouring and every colour class is an S(1,2,4).
3. k = 4 end to end: every proper 5-colouring of J(8,4) would be a partition into five
   S(3,4,8). We list all 30 S(3,4,8) on 8 labelled points and show at most 2 are
   pairwise block-disjoint, so chi(J(8,4)) > 5 (agrees with the known value 6).
4. Derived designs: the derived design of S(3,4,8) at any point is an S(2,3,7).
5. Bookkeeping: repeated derivation of S(k-1, k, 2k) reaches S(4,5,15) for k = 10
   and S(4,5,17) for k = 12.
"""
from itertools import combinations, permutations
from math import comb

def is_steiner(blocks, t, k, v):
    cover = {}
    for B in blocks:
        assert len(B) == k
        for T in combinations(sorted(B), t):
            cover[T] = cover.get(T, 0) + 1
    return len(cover) == comb(v, t) and all(c == 1 for c in cover.values())

# 1. clique lemma
for k in range(2, 7):
    n = 2 * k
    for L in combinations(range(n), k - 1):
        C = [frozenset(L) | {x} for x in range(n) if x not in L]
        assert len(C) == k + 1
        assert all(len(a & b) == k - 1 for a, b in combinations(C, 2))
print('1. clique lemma holds for k = 2..6')

# 2. k = 2
V = [frozenset(s) for s in combinations(range(4), 2)]
adj = lambda a, b: len(a & b) == 1
found = None
for col in __import__('itertools').product(range(3), repeat=len(V)):
    if all(col[i] != col[j] for i, j in combinations(range(len(V)), 2) if adj(V[i], V[j])):
        found = col; break
classes = [[V[i] for i in range(len(V)) if found[i] == c] for c in range(3)]
assert all(is_steiner(c, 1, 2, 4) for c in classes)
print('2. J(4,2) is 3-colourable and each class is an S(1,2,4):', [sorted(map(sorted, c)) for c in classes])

# 3. k = 4: all S(3,4,8), maximal pairwise-disjoint family
# extended Hamming code: weight-4 codewords of the [8,4,4] code
gens = [0b11110000, 0b11001100, 0b10101010, 0b11111111]
words = set()
for m in range(16):
    w = 0
    for i in range(4):
        if m >> i & 1: w ^= gens[i]
    words.add(w)
base = [frozenset(i for i in range(8) if w >> i & 1) for w in words if bin(w).count('1') == 4]
assert len(base) == 14 and is_steiner(base, 3, 4, 8)
systems = set()
for perm in permutations(range(8)):
    systems.add(frozenset(frozenset(perm[x] for x in B) for B in base))
systems = list(systems)
assert len(systems) == 30 and all(is_steiner(list(s), 3, 4, 8) for s in systems)
disjoint = lambda a, b: not (a & b)
best = 1
for r in range(2, 6):
    if any(all(disjoint(a, b) for a, b in combinations(fam, 2)) for fam in combinations(systems, r)):
        best = r
    else:
        break
print('3. there are', len(systems), 'S(3,4,8) on 8 points; at most', best,
      'are pairwise disjoint, so the 70 4-sets cannot split into 5 of them: chi(J(8,4)) > 5')
assert best < 5

# 4. derived design
for x in range(8):
    der = [B - {x} for B in base if x in B]
    assert is_steiner(der, 2, 3, 7)
print('4. every derived design of S(3,4,8) is an S(2,3,7)')

# 5. bookkeeping for k = 10, 12
for k, target in [(10, (4, 5, 15)), (12, (4, 5, 17))]:
    t, kk, v = k - 1, k, 2 * k
    chain = [(t, kk, v)]
    while t > 4:
        t, kk, v = t - 1, kk - 1, v - 1
        chain.append((t, kk, v))
    assert chain[-1] == target
    print('5. k = %d: ' % k + ' -> '.join('S(%d,%d,%d)' % c for c in chain))
print('ALL CHECKS PASSED')
