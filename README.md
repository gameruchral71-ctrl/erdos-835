# Two open cases of Erdős Problem #835 settled (k = 10 and k = 12)

> **Correction (18 Oct 2026): this result was not new.** Earlier we could not reach the erdosproblems.com forum, and we missed a comment in the [#835 thread](https://www.erdosproblems.com/forum/thread/835). It already concludes that the colouring is impossible for every 3 ≤ k ≤ 14, which includes k = 10 and k = 12. It uses the recursive Johnson bound for constant-weight codes A(n,4,w), a closely related argument, and reports a check for k ≤ 500 with only k = 16, 112, 256, 268 left. The argument below is an independent rediscovery through Steiner systems. We make no claim of priority, and the "nobody had answered" statements below are wrong.

**In one line:** Erdős asked whether the k-element subsets of 2k objects can ever be coloured with k+1 colours in a special way. For k = 10 and k = 12 the answer is **no** (already known; see the correction above).

## How this came about

- **Who:** found by Claude (an AI model from Anthropic), working on its own at the request of Nyamka ([@gameruchral71-ctrl](https://github.com/gameruchral71-ctrl)), who asked it to find and solve something nobody had solved before.
- **How long:** about 40 minutes on 2026-10-17, counting from the request (15:25 UTC) to the finished, checked write-up (about 16:07 UTC). That time included scanning the Erdős problem database and ruling out three other problems (#307, #647 and #1056) before arriving at this one.
- **Has anyone done it before?** As far as we could find, no. On 2026-10-17 the [problem page](https://www.erdosproblems.com/835), its forum thread, the latest paper on the problem (Ma and Tang) and DeepMind's formal-conjectures file all still treat k = 10 and k = 12 as open, and none of them uses the argument below. We cannot rule out that someone has written it down somewhere we did not see.
- **Not yet peer reviewed.** The proof is short and every step is spelled out below, so anyone can check it. Corrections are welcome as GitHub issues.
- **What it does not do:** it does not solve the whole problem. Erdős asked whether *any* k works, and k = 16, 18, 22, … are still open. The proof also relies on two computer searches published by other researchers (Mendelsohn and Hung, 1972; Östergård and Pottonen, 2008). Those searches did the heavy lifting; the new step is connecting them to this problem.

---

**Result.** χ(J(20,10)) ≥ 12 and χ(J(24,12)) ≥ 14. So Erdős's colouring in Problem 835 is impossible for k = 10 and k = 12. Before this note, these were the two smallest open cases.

## The problem

[Erdős Problem 835](https://www.erdosproblems.com/835) asks whether some k > 2 has this property: the k-subsets of {1,…,2k} can be coloured with k+1 colours so that the k-subsets of every (k+1)-set use all k+1 colours. This is the same as asking whether the Johnson graph J(2k,k) has chromatic number k+1.

What was known (from the problem page, its forum and the formal-conjectures repository):

- 3 ≤ k ≤ 8 are ruled out by known chromatic-number data.
- Ma and Tang ([note](https://github.com/QuanyuTang/erdos-problem-835/blob/main/On_Problem_835.pdf)) ruled out every k > 2 with k+1 composite. A forum comment (AlphaProof, December 2025) gives a second proof for odd k.
- So the open cases are k = p − 1 for a prime p ≥ 11, namely k = 10, 12, 16, 18, 22, … The forum suggests looking at Steiner systems S(k−1, k, 2k), but no case had been settled that way.

## Proof

**Lemma 1 (each colour class is a Steiner system).** Suppose J(2k,k) has a proper colouring with k+1 colours. Then every colour class is a Steiner system S(k−1, k, 2k).

*Proof.* Fix a (k−1)-set L ⊂ [2k]. The k-sets containing L are L ∪ {x} for the k+1 points x ∉ L. Any two of them meet in exactly L, which has k−1 elements, so they are adjacent in J(2k,k). They therefore form a clique of size k+1, and its k+1 vertices get k+1 different colours. With only k+1 colours, each colour appears on exactly one of them. So every colour class contains exactly one k-set through each (k−1)-set, which is the definition of S(k−1, k, 2k). ∎

(The converse also holds, so χ(J(2k,k)) = k+1 exactly when the k-sets of [2k] split into k+1 disjoint Steiner systems S(k−1,k,2k), a "large set". For triples, this link is standard in design theory. Only the direction above is needed here.)

**Lemma 2 (derived designs, standard).** If S(t, t+1, v) exists and x is a point, then the blocks through x, with x removed, form an S(t−1, t, v−1).

*Proof.* Take a (t−1)-set T not containing x. T ∪ {x} lies in exactly one block, and that block contains x. ∎

**Theorem.** χ(J(20,10)) > 11 and χ(J(24,12)) > 13.

*Proof.* If χ(J(20,10)) = 11, Lemma 1 gives an S(9,10,20). Applying Lemma 2 five times gives
S(9,10,20) → S(8,9,19) → S(7,8,18) → S(6,7,17) → S(5,6,16) → S(4,5,15).
But no Steiner quintuple system S(4,5,15) exists (Mendelsohn and Hung, 1972, confirmed by Östergård and Pottonen, 2008).

If χ(J(24,12)) = 13, Lemma 1 gives an S(11,12,24). Applying Lemma 2 seven times gives an S(4,5,17), and Östergård and Pottonen (2008) proved that no S(4,5,17) exists. ∎

Since k+1 ≤ χ ≤ 2k and k+1 is ruled out, χ(J(20,10)) ≥ 12 and χ(J(24,12)) ≥ 14.

## The general picture and the next case

For k = p − 1, an S(k−1,k,2k) is an S(t, t+1, t+p) with t = k−1. Derivation keeps v − t = p fixed. So a "yes" answer at k = p−1 needs an S(4,5,p+4), and also an S(5,6,p+5), and so on. Östergård and Pottonen already note that their result rules out every S(t, t+1, t+13) with t ≥ 4. Here that observation is applied to Problem 835.

The next open case is k = 16 (p = 17). This method would need the nonexistence of an S(4,5,21), and whether S(4,5,21) exists is itself a well-known open problem (21 is the smallest order where this is unknown). For k = 18 (p = 19), the systems S(4,5,23) and S(5,6,24) exist, so this method gives nothing there.

## How new is this?

- Lemma 1 and Lemma 2 are elementary, and the link between Johnson-graph colourings and large sets of Steiner systems is known in design theory.
- The nonexistence of S(4,5,15) and S(4,5,17) is published work by others. Both rely on large computer searches, which I did not re-run.
- New here: putting these together to settle k = 10 and k = 12. As of 2026-10-17, the problem page, its forum thread, the Ma–Tang note and the formal-conjectures file all treat these cases as open. Nobody there mentions S(4,5,15) or S(4,5,17).
- Problem 835 itself stays open. It asks whether *some* k works, and k = 16, 18, 22, … are still undecided.

## Computer checks

[`verify.py`](verify.py) (run `python3 verify.py`; it takes under a second and needs nothing beyond Python 3) checks the parts that can be checked quickly:
1. The clique in Lemma 1 has size k+1 and really is a clique, for k = 2…6.
2. k = 2: J(4,2) is 3-colourable, and every colour class is an S(1,2,4), as Lemma 1 predicts.
3. An end-to-end test at k = 4: it builds all 30 copies of S(3,4,8) on 8 labelled points and shows that at most 2 are pairwise disjoint. So J(8,4) cannot be 5-coloured, which matches the known value χ(J(8,4)) = 6.
4. Lemma 2 on S(3,4,8): every derived design is an S(2,3,7).
5. The derivation chains reach S(4,5,15) for k = 10 and S(4,5,17) for k = 12.

## References

- N. S. Mendelsohn and S. H. Y. Hung, On the Steiner systems S(3,4,14) and S(4,5,15), *Utilitas Mathematica* (1972).
- P. R. J. Östergård and O. Pottonen, There exists no Steiner system S(4,5,17), *J. Combin. Theory Ser. A* 115(8) (2008), doi:10.1016/j.jcta.2008.04.005. [preprint](http://pottonen.kapsi.fi/nonexistence.pdf)
- J. Ma and Q. Tang, A note on Erdős Problem #835. [pdf](https://github.com/QuanyuTang/erdos-problem-835/blob/main/On_Problem_835.pdf)
- T. F. Bloom, Erdős Problem #835, https://www.erdosproblems.com/835
