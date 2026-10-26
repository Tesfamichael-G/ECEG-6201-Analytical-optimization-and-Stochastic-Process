# Chapter 1 — Optimization: Preliminaries

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 tutorial · **Home study:** 5 h
    
    **Coverage (syllabus):** characterization of a minimum · local vs. global minima ·
    Pareto optimality · gradient and Hessian · slope and curvature · properties of
    linear and quadratic functions · optimality conditions (unconstrained, linear
    equality, linear inequality, nonlinear) · line-search subproblem (bracketing and
    sectioning) · numerically and analytically derived gradients · adjoint-variable
    approach.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Characterize a *local*, *strict*, and *global* minimizer and recognize them on a plot.
2. Define Pareto optimality and identify the Pareto front of a small multi-objective problem.
3. Compute gradients and Hessians analytically; interpret them as slope and curvature.
4. State and use the first- and second-order optimality conditions (unconstrained and constrained, including KKT).
5. Solve the line-search subproblem via bracketing and golden-section sectioning.
6. Estimate gradients numerically and explain when to prefer analytical or adjoint derivations.
7. Implement basic gradient/Hessian checks and golden-section search in Python, Java, or C#.

---

## 1.1 Characterization of a minimum

Let \(f:\mathbb{R}^n \to \mathbb{R}\) be a (smooth) objective and let \(S \subseteq \mathbb{R}^n\) be the feasible set.

!!! note "Definitions"
    - **Local minimizer.** \(x^\star \in S\) is a *weak* local minimizer if there is \(\varepsilon>0\) with
      \(f(x^\star) \le f(x)\) for all \(x \in S\), \(\|x-x^\star\|<\varepsilon\).
      If the inequality is strict for \(x \ne x^\star\), \(x^\star\) is a **strict local minimizer**.
    - **Global minimizer.** \(f(x^\star) \le f(x)\) for all \(x \in S\).
    - **Isolated minimizer.** A neighborhood of \(x^\star\) contains no other local minimizer.

A point that satisfies the necessary condition \(\nabla f(x^\star)=0\) is a **stationary point** — it may be a minimizer, a maximizer, or a saddle.

### Worked example 1.1 — three critical points of a cubic

Minimize \(f(x) = x^3 - 3x + 1\) over \(x \in \mathbb{R}\).

1. **First-order condition.**
   \[
     f'(x) = 3x^2 - 3 = 3(x-1)(x+1) = 0 \;\Rightarrow\; x = \pm 1 .
   \]
2. **Second-order condition.**
   \[
     f''(x) = 6x,\qquad f''(1) = 6>0,\qquad f''(-1)=-6<0.
   \]
3. **Classification.**
   - \(x = 1\): strict local minimum with \(f(1) = 1 - 3 + 1 = -1\).
   - \(x = -1\): strict local maximum with \(f(-1) = -1 + 3 + 1 = 3\).
4. **Global behavior.** Because \(f(x) \to +\infty\) as \(x \to +\infty\) and \(f(x) \to -\infty\) as \(x \to -\infty\), neither a global minimum nor a global maximum exists on \(\mathbb{R}\).

!!! tip "Takeaway"
    A stationary point is *not automatically* a minimum. The Hessian (here the scalar \(f''\))
    decides between minimum, maximum, and saddle.

---

## 1.2 Local versus global minima

Local minima are **neighbourhood** statements — the comparison is only with nearby points. To decide *global* status you either

- prove **convexity** of \(f\) on a convex \(S\) (then every local minimizer is global), or
- enumerate all local minima and compare their values (only feasible for small or structured problems), or
- use a certified global solver (branch-and-bound, interval methods, etc.).

### Worked example 1.2 — a quartic with two global minima

Minimize \(f(x) = x^4 - 2x^2 + 1 = (x^2-1)^2\) over \(\mathbb{R}\).

1. \(f'(x) = 4x(x^2 - 1) = 0\) ⇒ critical points \(x = 0, \pm 1\).
2. \(f''(x) = 12x^2 - 4\).
   - \(x = 0\): \(f''(0) = -4 < 0\), so \(x=0\) is a **local maximum** (\(f(0)=1\)).
   - \(x = \pm 1\): \(f''(\pm 1) = 8 > 0\), so both are **strict local minima** with \(f(\pm 1) = 0\).
3. Since \(f(x) \ge 0\) everywhere and \(f(\pm 1)=0\), **both local minima are global**.
   The **global minimizer is not unique**; the set of minimizers is \(\{-1, +1\}\).

### Worked example 1.3 — a multimodal function where locals are *not* global

For \(f(x) = 0.5 x^2 + 3\sin(3x)\) on \([-3, 3]\) the derivative

\[
f'(x) = x + 9\cos(3x)
\]

has several zeros (see the Python scan below). Two local minima exist, but only one is global.

```python
import numpy as np

f  = lambda x: 0.5*x**2 + 3*np.sin(3*x)
df = lambda x: x + 9*np.cos(3*x)

# coarse scan for sign changes of f'
xs = np.linspace(-3, 3, 4001)
sign_change = np.where(np.diff(np.sign(df(xs))) != 0)[0]
roots = [ (xs[i]+xs[i+1])/2 for i in sign_change ]
print("stationary points:", roots)
print("values:", [round(f(r), 4) for r in roots])
