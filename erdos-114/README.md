# Erdős Problem #114 in degree 3: $z^3-1$ maximises the lemniscate length

**Claim.** For every monic cubic $p$, the length of the lemniscate $E(p)=\{z\in\mathbb C:|p(z)|=1\}$ is at most the length for $z^3-1$:
$$
\ell(E(p))\le \ell(E(z^3-1)) = 2^{1/3}B\!\left(\tfrac16,\tfrac12\right)=9.17972422234315\ldots
$$
So the Erdős–Herzog–Piranian conjecture ([Erdős Problem #114](https://www.erdosproblems.com/114), \$250 prize for all degrees) holds for $n=3$.

**Status.** This is a computer-assisted proof. It has **not been peer-reviewed**. It settles one degree only, so the problem (and its prize) stays open.

**Previously known.** $n=1$ is trivial. $n=2$ is due to Eremenko–Hayman (1999). Tao ([arXiv:2512.12455](https://arxiv.org/abs/2512.12455), Dec 2025) proved the conjecture for all *sufficiently large* $n$, with no explicit threshold. To my knowledge, $n=3$ was open as of 17 Oct 2026:
- the erdosproblems.com database and the Lean statement (google-deepmind/formal-conjectures, updated 2026-10-09) list only $n=1$, $n=2$ and large $n$;
- the [MendozaLab census for $3\le n\le 14$](https://github.com/MendozaLab/erdos-experiments) was publicly challenged ([issue #4](https://github.com/MendozaLab/erdos-experiments/issues/4)) and its author now describes it as "a computational census, not a proof".

The erdosproblems.com forum thread could not be read from the environment used here. It should be checked before this result is announced.

**Credits.** Project initiated, directed and published by Uchral Nyam-Osor. The mathematics and code were done by Claude (Anthropic's AI model), in Oct 2026.

---

## 1. Reduction to a one-parameter family

**Eremenko–Hayman** ([Michigan Math. J. 46 (1999)](https://arxiv.org/abs/0805.2295)):
- **Lemma 4:** a maximiser of $\ell(E(p))$ over monic degree-$n$ polynomials exists.
- **Lemma 5:** some maximiser has all of its critical points on $E(p)$, i.e. all of its critical values lie on the unit circle $\mathbb T$.

(This is also Proposition 1.2 in Tao's paper.)

Length is unchanged by translating $z$ and by replacing $p(z)$ with $\omega^{-n}p(\omega z)$ for $|\omega|=1$. So for $n=3$ we may assume $p(z)=z^3+az+b$.

**The double-critical-point case.** If $p$ has a double critical point, then $a=0$ and $p=z^3+b$ with $|b|=1$. This is a rotated copy of $z^3-1$.

**The two-critical-point case.** Otherwise the critical points are $\pm c$ with $c^2=-a/3$, and $p(\pm c)=b\pm\tfrac{2a}{3}c$.
- Rotating, we may take the critical values to be $e^{\pm i\theta/2}$ with $0<\theta\le\pi$.
- Solving gives
  $$
  p_\theta(z)=z^3+az+b,\qquad b=\cos\tfrac\theta2,\qquad a^{3/2}=\alpha:=\tfrac{3\sqrt3}{2}\sin\tfrac\theta2\in(0,\tfrac{3\sqrt3}{2}].
  $$
- The three cube-root choices for $a$ give the same polynomial up to $z\mapsto e^{2\pi i/3}z$.

**Conclusion of the reduction.** It suffices to prove
$$
L(\theta):=\ell(E(p_\theta)) < L_0:=\ell(E(z^3-1))\qquad(0<\theta\le\pi). \tag{$*$}
$$
Then the maximiser from Lemma 5 is (a rotation and translation of) $z^3-1$, and every monic cubic has length at most $L_0$.

## 2. Length formula and rescaling

**Length formula.** Away from the finitely many critical points, $p$ maps $E(p)$ 3-to-1 onto $\mathbb T$ and stretches arc length by the factor $|p'|$. So
$$
\ell(E(p))=\int_{-\pi}^{\pi}\sum_{p(z)=e^{it}}\frac{dt}{|p'(z)|}.
$$
For $p_0=z^3+1$ this is $\int|e^{it}-1|^{-2/3}dt=L_0$.

**Rescaling.** Put $z=\sqrt a\,\zeta$ and $e^{it}=b+\alpha\eta$, so that $p_\theta(z)=\alpha(\zeta^3+\zeta)+b$. Then
$$
L(\theta)=\sqrt a\int_{C_\alpha}h(\eta)\,|d\eta|,\qquad h(\eta):=\sum_{\zeta^3+\zeta=\eta}\frac1{|3\zeta^2+1|},
$$
where $C_\alpha$ is the circle $|b+\alpha\eta|=1$.
- $C_\alpha$ passes through the two critical values $\pm i\sigma$ of $q(\zeta)=\zeta^3+\zeta$, where $\sigma=2/(3\sqrt3)$.
- The critical angle $t_c=\theta/2$ satisfies $\sin t_c=\alpha\sigma$ and $b=\cos t_c$.

## 3. The cusp at $\theta=0$ and the exact constant $2D$

**The cusp.** As $\theta\to0$, $L(\theta)=L_0-2D\sqrt a+o(\sqrt a)$ (numerically, the error looks like $O(a^3)$), with $2D=0.938307085876\ldots$,
where
$$
D=\int_{-\infty}^{\infty}\Big(\tfrac{2}{\sqrt3}-\sqrt{1+\tfrac{x^2}{3(1+x^2)}}\Big)dx .
$$

**Geometric meaning.** $2D$ is the length deficit of the hyperbola $3y^2-x^2=1$ against its asymptotes. Near the origin, $E(p_\theta)$ looks like $\sqrt a\cdot\{\operatorname{Re}(\zeta^3+\zeta)=0\}$, which is the imaginary axis plus that hyperbola. For $\theta=0$ it looks like three lines through $0$ instead.

Because of this cusp, $L$ has infinite slope at $\theta=0$, so naive Lipschitz-based checks cannot certify a neighbourhood of $z^3-1$. We avoid the problem with an exact identity that is regular at $\alpha=0$.

**The regular identity.** Parametrise $C_\alpha$ by
$$
\eta=\eta_0(\tau)+\varepsilon,\qquad \eta_0(\tau)=\frac{i\tau}{1-i\alpha\tau/2},\qquad \varepsilon=\frac{1-b}{\alpha}=\frac{4\alpha/27}{1+b},\qquad \tau\in\mathbb R,
$$
which is the Cayley substitution $t=2\arctan(\alpha\tau/2)$, with $|d\eta|=w(\tau)\,d\tau$ and $w=1/(1+\alpha^2\tau^2/4)$. Then
$$
G(\alpha):=\frac{L_0-L(\theta)}{\sqrt a}=2\int_0^\infty\Big[|\eta_0(\tau)|^{-2/3}-h(\eta_0(\tau)+\varepsilon)\Big]w(\tau)\,d\tau ,
$$
and $G(0)=2D$. Condition $(*)$ is equivalent to $G(\alpha)>0$ for $0<\alpha\le 3\sqrt3/2$.

## 4. Rigorous verification

All numerics use ball arithmetic (Arb, through python-flint 0.9, 80-bit precision). Only upper bounds on $h$ are ever needed.

**(a) Root enclosures.** For a complex ball of values $\eta$:
- each root of $\zeta^3+\zeta=\eta$ is enclosed in a box proven by the Krawczyk test to contain a root for every $\eta$ in the ball;
- the three boxes are checked to be disjoint, so they hold all three roots;
- then $h\le\sum 1/\min|3\zeta^2+1|$ over the boxes.

**(b) Near a critical value** ($|\eta-i\sigma|\le d$), using $q(\zeta)-i\sigma=(\zeta-c)^2(\zeta+2c)$ with $c=i/\sqrt3$ and Rouché on two small disks:
$$
h(\eta)\le \frac{2K_1}{\sqrt{|\eta-i\sigma|}}+K_3,
$$
with explicit $K_1\approx 3^{1/4}/(2\sqrt3)$ and $K_3\approx1/3$. This integrable bound covers a short window on each side of the critical parameter.

**(c) Tail lemma.** For $|\eta|\ge R$:
$$
|h(\eta)-|\eta|^{-2/3}|\le C(R)\,|\eta|^{-2},
$$
with $C(16)\le 0.444$. The proof uses Rouché around $u-\tfrac1{3u}$ ($u^3=\eta$), with radius $1/(60|u|^5)$, plus the cancellation $\sum_k \operatorname{Re}(u_k^{-2})=0$.

**(d) Part 1** (`part1_small_alpha.py`), for $\alpha\in[0,0.02]$ (40 boxes). Here $G$ is bounded below by:
- an exact lower bound for the $|\eta_0|^{-2/3}$ term;
- minus adaptive upper bounds (a)/(b) on $[0,\tau_c+16]$;
- minus the tail (c).

The result is $G(\alpha)\ge 0.72>0$ on every box. The true value is about $0.938$.

**(e) Part 2** (`part2_large_alpha.py`), for $t_c\in[\arcsin(0.02\sigma),\pi/2]$, i.e. $\alpha\in[0.02,3\sqrt3/2]$, in 536 geometric boxes in $t_c$ (ratio 1.01, run as 4 contiguous chunks split at the values in `cuts.txt`; no box needed splitting). Here
$$
L(\theta)=\frac2a\int_0^\pi h(\eta(t))\,dt
$$
is bounded above by (a) and (b). The variable is written relative to $t_c$, as $\eta-i\sigma=e^{it_c}(e^{i(t-t_c)}-1)/\alpha$, which avoids the dependency blow-up near the critical value. Every box gives $L\le$ (upper bound) $<L_0$.

**(f) Sanity checks** (`sanity_checks.py`). Random floating-point tests of lemmas (b) and (c) and of the enclosures (a).

**Results of the run (17 Oct 2026):** Part 1: $G(\alpha)\ge0.7234$ on all 40 boxes. Part 2: every one of the 536 boxes gives $L\le 9.1411<L_0=9.17972\ldots$ (the closest box is at $\alpha\approx0.02$). Logs are in `logs/` (the chunk logs print `pi/2` as the end of every chunk; each chunk actually ends at the next value in `cuts.txt`). To reproduce: `pip install python-flint==0.9.0 numpy mpmath`, then `python3 part1_small_alpha.py 0.02 40` and `python3 part2_large_alpha.py 0.02 1.01` (single process, about 80 min).

## 5. What this does and doesn't show

**Shows.** $z^3-1$ maximises the lemniscate length among monic cubics, so the conjecture holds in degree 3.

**Also shows.** Among cubics whose critical values lie on $\mathbb T$, $z^3-1$ is the unique maximiser up to rotation and translation.

**Does not show.**
- Uniqueness among *all* maximisers. Lemma 5 produces only one maximiser of that form.
- Anything for degree $\ge4$. There the same reduction gives a $(n-2)$-parameter family with several combinatorial types.

**Dependence.** The proof relies on Eremenko–Hayman Lemmas 4–5, which use quasiconformal deformations.

## Files

| file | role |
|---|---|
| `rigorous_common.py` | Krawczyk root enclosures, near-critical bound (b) |
| `part1_small_alpha.py` | region $\alpha\le0.02$, tail lemma (c) |
| `part2_large_alpha.py` | region $\alpha\ge0.02$ |
| `sanity_checks.py` | independent float checks of the lemmas |
| `explore3b.py`, `proto_G.py` | non-rigorous exploration ($L(\theta)$ table, $G(\alpha)$) |
| `logs/` | full output of the verification runs |
