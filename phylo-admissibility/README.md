# Maximum parsimony is inadmissible: one of Mike Steel's US$100 phylogenetics conjectures

**Credit.** The idea of tackling an open problem in biology came from **Tugsuu** (Nyamka's sister).
The literature search, the computations and the proofs below were done by **Claude** (an AI model by
Anthropic) in a project run by Nyamka. Nothing here has been reviewed by an expert yet.

## The problem

Mike Steel (University of Canterbury) lists "Phylogenetic inadmissibility conjectures" among his
$100 challenges (prize: US$100, or NZ$100 plus a bottle of NZ wine, for the first correct solution):
<https://www.math.canterbury.ac.nz/~m.steel/challenges.shtml>, statement at
<https://www.math.canterbury.ac.nz/~m.steel/Non_UC/files/research/conjecture2.pdf>, and item SC3 of his
"Phylogenetics: Challenges and Conjectures" (updated July 2018), where it is still listed as open.

Setting: k independent sites evolve under a finite-state Markov model (for example Jukes–Cantor) on
a fully resolved tree T with branch lengths λ. A tree-reconstruction method M is **inadmissible** if
some other method M′ has P(M′ = T | T, λ) ≥ P(M = T | T, λ) for every tree T and every positive,
finite λ, with strict inequality for at least one (T, λ).

* Conjecture 1: Maximum Parsimony (MP) is inadmissible.
* Conjecture 2: Maximum Likelihood (ML) is inadmissible.

A search of journals, arXiv/bioRxiv listings reachable from here, forums and news (October 2026)
found no proof or disproof of either since they were posed (2006).

## Result 1 (proved, exact): MP is inadmissible under Jukes–Cantor with 4 taxa and 2 sites

Take MP with ties broken uniformly at random. Let M* agree with MP except on the following 2-site
data sets (in every one of them MP ties all three trees, so MP picks each with probability 1/3).
A site of type "AABC" means two leaves share a state and the other two leaves each have their own,
different state.

| data (2 sites) | M* puts probability |
|---|---|
| one constant site + one AABC site whose identical pair is {a,b} | 2/3 on ab\|cd, 1/6 on each other tree |
| two singleton sites (one leaf differs from the rest) at different leaves i, j | 1/9 on ij\|kl, 4/9 on each other tree |
| two AABC sites with overlapping identical pairs {a,b}, {a,c} | 7/18 on ad\|bc (the tree neither supports), 11/36 on each other |
| two AABC sites with complementary identical pairs {a,b}, {c,d} | 7/18 on ab\|cd, 11/36 on each other |

**Theorem.** For each of the three trees T and every choice of the five branch lengths,
P(M* = T) − P(MP = T) ≥ 0, and it is strictly positive on an open set of branch lengths.

So MP is inadmissible. `verify_jc4.py` proves this exactly: the difference is a polynomial in
θ_e = exp(−4t_e/3) ∈ (0,1), and all its Bernstein coefficients on [0,1]^5, computed in exact rational
arithmetic, are ≥ 0 and not all 0. The largest gain is about 2.1 percentage points of
probability of recovering the true tree (average about 0.3 points over random branch lengths).
`simulate_jc4.py` independently confirms the gain by Monte Carlo simulation.

M* is itself a version of MP with a different tie-breaking rule: it uses information that
parsimony ignores (sites that are parsimony-uninformative but not likelihood-uninformative).

## Result 2 (proved, exact): the same holds under the 2-state (CFN) model with 5 taxa and 3 sites

`verify_cfn5.py` checks an explicit rule that differs from MP on 210 of the 3-site data sets and
beats it for all 15 trees and all 7 branch lengths (again by exact Bernstein coefficients). Here the
dominating rule also sometimes leaves the unique MP tree, so it is not just a different tie-break.
For the 2-state model with 4 taxa, a linear-programming search found **no** rule dominating MP for
k = 2, 3, 4 sites, so MP appears admissible there: whether MP is inadmissible depends on the size
of the problem and on the model.

## What is not proved

* **Other tie-breaking rules for MP.** The theorems are about MP with ties broken uniformly at random
  (the usual convention; if ties count as failures, MP is trivially inadmissible). Every one of 8 random
  alternative tie-breaking rules tested under Jukes–Cantor (4 taxa, 2 sites) and 12 under CFN
  (5 taxa, 4 sites) was also dominated in the linear-programming search, but a proof for every possible
  tie-breaking rule is open.
* **Maximum likelihood (Conjecture 2): strong evidence, not a proof.** Under CFN with 4 taxa and 4
  sites, `verify_ml_cfn4.py` proves exactly that a rule M* beats ML for all trees and branch lengths
  **provided** that on three kinds of data ML ties all three trees. For each of those, every tree
  reaches the same likelihood at an explicit point (exact arithmetic), and 400 random restarts per tree
  never exceed it, but the upper bound (that no tree does better) is only numerical. The maximisers lie on
  curves on the boundary of parameter space, which is why a simple certificate fails. With that gap
  closed, ML would also be inadmissible. Under Jukes–Cantor, the search finds a rule beating ML at 4 taxa
  and 3 sites, also resting on ML ties.

## Files

* `verify_jc4.py`: exact proof of Result 1 (seconds).
* `verify_cfn5.py`: exact proof of Result 2 (about 2 minutes).
* `verify_ml_cfn4.py`: the conditional ML result and the numerical tie check.
* `simulate_jc4.py`: Monte Carlo check of Result 1.
* `search/`: linear-programming search code that found the rules (equivariant LP over site-pattern data).

Requirements: Python 3 with `sympy` (verification) and `numpy`, `scipy` (search, simulation).
