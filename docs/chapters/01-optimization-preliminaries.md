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
```

The smallest computed value is the **global** minimum; the others are **local** minima
(local but not global). This is the *typical* situation in engineering design problems.

---

## 1.3 Pareto-optimality in multi-objective optimization

When we minimize **several** objectives at once, there is rarely a single point that improves all of them simultaneously.

!!! note "Definitions"
    Let \(\mathbf f(x) = (f_1(x),\dots,f_m(x))\). A point \(x^\star\) **Pareto-dominates** \(y\) if
    \[
      f_i(x^\star) \le f_i(y) \ \forall i,\qquad\text{and } f_j(x^\star) < f_j(y) \ \text{for some }j.
    \]
    \(x^\star\) is **Pareto-optimal** if no feasible \(y\) dominates it.
    The set of objective vectors of all Pareto-optimal points is the **Pareto front**.

### Worked example 1.4 — the Pareto set of a two-objective problem

\[
\min_{x\in\mathbb{R}}\ \bigl(f_1(x), f_2(x)\bigr) = \bigl(x^2,\ (x-2)^2\bigr).
\]

- If \(x < 0\), moving to \(x=0\) reduces **both** objectives ⇒ every \(x<0\) is dominated.
- If \(x > 2\), moving to \(x=2\) reduces both ⇒ every \(x>2\) is dominated.
- If \(0 \le x \le 2\): moving right increases \(f_1\) but decreases \(f_2\); moving left does the opposite. **No point in \([0,2]\) is dominated.**

Hence the **Pareto set** is \([0,2]\) and the **Pareto front** is

\[
\mathcal F = \Bigl\{\,(x^2,\ (x-2)^2) : 0 \le x \le 2 \,\Bigr\}.
\]

On the front:

| \(x\) | \(f_1=x^2\) | \(f_2=(x-2)^2\) |
|:---:|:---:|:---:|
| 0.0 | 0.00 | 4.00 |
| 0.5 | 0.25 | 2.25 |
| 1.0 | 1.00 | 1.00 |
| 1.5 | 2.25 | 0.25 |
| 2.0 | 4.00 | 0.00 |

!!! tip "Engineering usage"
    When you face a trade-off between, say, **link throughput** and **energy per bit**, you
    do not look for a single number — you look for the Pareto front and let the design
    authority select a knee point on it.

---

## 1.4 Gradient and Hessian matrix

For \(f:\mathbb{R}^n \to \mathbb{R}\),

\[
\nabla f(x) = \begin{bmatrix} \partial f/\partial x_1 \\ \vdots \\ \partial f/\partial x_n\end{bmatrix},
\qquad
\nabla^2 f(x) = H(x) = \begin{bmatrix}
\dfrac{\partial^2 f}{\partial x_i \partial x_j}
\end{bmatrix}_{i,j=1}^{n}.
\]

\(H\) is symmetric when \(f\) has continuous second derivatives (Schwarz's theorem).

### Worked example 1.5 — gradient and Hessian of a polynomial

Let \(f(x_1, x_2) = x_1^2 x_2 + x_2^3 - 2 x_1\). Then

\[
\nabla f(x_1,x_2) = \begin{bmatrix} 2x_1 x_2 - 2 \\ x_1^2 + 3x_2^2 \end{bmatrix},
\qquad
H(x_1,x_2) = \begin{bmatrix} 2x_2 & 2x_1 \\ 2x_1 & 6x_2 \end{bmatrix}.
\]

At the point \((1,2)\): \(\nabla f = (2, 13)^\top\), \(H = \begin{bmatrix} 4 & 2 \\ 2 & 12\end{bmatrix}\).
The eigenvalues of \(H\) are \(3\pm\sqrt{5}>0\), so \(H\succ 0\) there — useful when we later use Newton's method (Chapter 2).

---

## 1.5 Slope and curvature

- The **directional derivative** of \(f\) at \(x\) along a unit vector \(d\) is \(\nabla f(x)^\top d\) — the *slope*.
- The **curvature** along \(d\) is \(d^\top H(x)\, d\):
  - \(d^\top H d > 0\): \(f\) curves *up* along \(d\).
  - \(d^\top H d < 0\): \(f\) curves *down* along \(d\).
  - \(d^\top H d = 0\): locally flat along \(d\) to second order.

The **level set** \(f(x) = \text{const}\) near a regular point is perpendicular to \(\nabla f(x)\), so moving along a *level set* is a zero-slope direction. This is the geometric basis of constrained optimality (Section 1.8).

---

## 1.6 Properties of linear functions

For \(f(x) = a^\top x + b\) with \(a \in \mathbb{R}^n\):

| Property | Value |
|---|---|
| Gradient | \(\nabla f(x) = a\) (constant) |
| Hessian | \(H(x) = 0\) |
| Level sets \(f=\text{const}\) | Hyperplanes orthogonal to \(a\) |
| Local minima | Only if \(a = 0\) (then every point is a minimum) |
| Global minima on a polytope | Attained on a face / at a vertex (LPs) |

A linear objective has **no interior unconstrained minimum** — you need constraints (Chapter 3) to get a well-posed problem.

---

## 1.7 Properties of quadratic functions

A quadratic function has the form

\[
f(x) = \tfrac12 x^\top Q x - c^\top x + d,
\qquad Q = Q^\top \in \mathbb{R}^{n\times n},\quad c\in\mathbb{R}^n,\ d\in\mathbb{R}.
\]

| Property | Value |
|---|---|
| Gradient | \(\nabla f(x) = Qx - c\) |
| Hessian | \(H(x) = Q\) (constant) |
| Stationary point | \(x^\star = Q^{-1}c\) (if \(Q\) invertible) |
| Unique global minimum | Exists ⇔ \(Q \succ 0\) |
| Level sets (with \(Q \succ 0\)) | Ellipsoids aligned with eigenvectors of \(Q\); semi-axis lengths \(\propto 1/\sqrt{\lambda_i}\) |
| Condition number \(\kappa(Q)\) | Determines how difficult the problem is for iterative methods |

### Worked example 1.6 — solve a quadratic via the eigen-structure

Minimize \(f(x) = 2x_1^2 + x_2^2 + x_1 x_2 - 3 x_1 + 4 x_2\).

Write it as \(\tfrac12 x^\top Q x - c^\top x\) with

\[
Q = \begin{bmatrix} 4 & 1 \\ 1 & 2\end{bmatrix},
\qquad
c = \begin{bmatrix} 3 \\ -4 \end{bmatrix}.
\]

**Eigenvalues of \(Q\):** \(\det(Q-\lambda I) = (4-\lambda)(2-\lambda) - 1 = \lambda^2 - 6\lambda + 7 = 0 \Rightarrow \lambda = 3 \pm \sqrt{2}\), both positive ⇒ \(Q \succ 0\), so \(f\) has a **unique global minimum**.

**Stationary point** solves \(Qx = c\):

\[
\begin{cases} 4x_1 + x_2 = 3 \\ x_1 + 2x_2 = -4 \end{cases}
\ \Rightarrow\ x^\star = \Bigl(\tfrac{10}{7}, -\tfrac{19}{7}\Bigr).
\]

**Optimal value** (using \(f^\star = -\tfrac12 c^\top Q^{-1} c\)):

\[
f(x^\star) = -\tfrac{53}{7}\approx -7.571.
\]

The condition number is \(\kappa(Q) = \frac{3+\sqrt2}{3-\sqrt2}\approx 1.97\) — well-conditioned; a gradient method would converge reasonably fast.

---

## 1.8 Optimality conditions

### 1.8.1 Unconstrained

!!! abstract "First-order necessary condition (FONC)"
    If \(x^\star\) is a local minimizer of \(f\) and \(f\) is differentiable in a neighbourhood, then
    
    \[
      \nabla f(x^\star) = 0 .
    \]
    
    If in addition \(f\in C^2\),
    
    \[
      \nabla^2 f(x^\star) \succeq 0
    \]
    
    (the Hessian is positive semi-definite) is **necessary**. The strict version
    \(\nabla^2 f(x^\star) \succ 0\) is **sufficient** for a strict local minimum.

### Worked example 1.7 — first-order and second-order conditions

For \(f(x_1,x_2)=x_1^2 x_2 + x_2^3 - 2x_1\) find all stationary points and classify them.

\[
\nabla f = 0 \;\Longleftrightarrow\;
\begin{cases} 2x_1 x_2 = 2 \\ x_1^2 + 3x_2^2 = 0 \end{cases}
\]

The second equation forces \(x_1 = x_2 = 0\), which then violates \(2x_1 x_2 = 2\). **No real stationary point** exists ⇒ \(f\) is unbounded below. This is an important diagnostic: an optimizer that diverges is telling you the problem is unbounded.

### 1.8.2 Linearly constrained — equality constraints

Minimize \(f(x)\) subject to \(Ax = b\), \(A\in\mathbb{R}^{m\times n}\) of rank \(m\).

Introduce a **Lagrange multiplier vector** \(\lambda\in\mathbb{R}^m\) and the Lagrangian

\[
\mathcal L(x,\lambda) = f(x) + \lambda^\top (Ax - b).
\]

!!! abstract "Equality-constrained FONC"
    If \(x^\star\) is a local minimizer of \(f\) on \(\{x: Ax=b\}\), there exists \(\lambda^\star\) with
    
    \[
      \nabla_x \mathcal L(x^\star,\lambda^\star) = \nabla f(x^\star) + A^\top \lambda^\star = 0,
      \qquad Ax^\star = b.
    \]

Geometric interpretation: at \(x^\star\), the negative gradient \(-\nabla f(x^\star)\) lies in the row space of \(A\) — equivalently, \(\nabla f(x^\star)\) is orthogonal to the *tangent space* of the constraint (which is \(\ker A\)).

### Worked example 1.8 — projection onto an affine subspace

Minimize \(f(x_1, x_2) = x_1^2 + x_2^2\) subject to \(x_1 + x_2 = 2\).

Lagrangian: \(\mathcal L = x_1^2 + x_2^2 + \lambda(x_1 + x_2 - 2)\).
First-order conditions:
\[
2x_1 + \lambda = 0,\qquad 2x_2 + \lambda = 0,\qquad x_1 + x_2 = 2.
\]
The first two imply \(x_1 = x_2\); with the constraint \(x_1 = x_2 = 1\). Then \(\lambda = -2\), and
\[
f^\star = 1 + 1 = 2.
\]
Geometric check: \(-\nabla f(1,1) = (-2,-2) = 2\cdot(1,1)\) — the negative gradient is (up to scale) the constraint normal, as required.

### 1.8.3 Linearly constrained — inequality constraints

Minimize \(f(x)\) subject to \(Ax \le b\).

The **Karush–Kuhn–Tucker (KKT) conditions** for this case read
\[
\nabla f(x^\star) + A^\top \mu^\star = 0,\qquad
\mu^\star \ge 0,\qquad
\mu_i^\star\,[Ax^\star - b]_i = 0,\qquad
Ax^\star \le b.
\]
The last column condition is **complementary slackness**: an active constraint carries a non-negative multiplier; an inactive one carries zero.

### Worked example 1.9 — a small inequality problem

Minimize \((x-2)^2 + 1\) subject to \(x \le 1\).

- **Unconstrained solution:** \(x = 2\), but \(2 \le 1\) is false ⇒ infeasible.
- **Active-constraint ansatz:** \(x^\star = 1\). Then \(f'(1) = 2(1-2) = -2\), and KKT with \(g(x)=x-1\le 0\), \(\nabla g = 1\):
  \[
  -2 + \mu \cdot 1 = 0 \Rightarrow \mu = 2 \ge 0. \;\checkmark
  \]
  \(x^\star = 1\) satisfies all KKT conditions. **Optimal value:** \((1-2)^2+1 = 2\).

### 1.8.4 Nonlinearly constrained — full KKT

Minimize \(f(x)\) subject to \(g_i(x) \le 0\ (i=1,\dots,p)\) and \(h_j(x) = 0\ (j=1,\dots,q)\), all \(C^1\).

!!! abstract "KKT conditions (necessary, under a constraint qualification)"
    There exist multipliers \(\mu^\star\in\mathbb{R}^p\) and \(\lambda^\star\in\mathbb{R}^q\) such that
    \[
      \nabla f(x^\star) + \sum_i \mu_i^\star \nabla g_i(x^\star) + \sum_j \lambda_j^\star \nabla h_j(x^\star) = 0,
    \]
    \[
      g_i(x^\star) \le 0,\ \ \mu_i^\star \ge 0,\ \ \mu_i^\star g_i(x^\star) = 0,
      \qquad
      h_j(x^\star) = 0.
    \]

Common constraint qualifications: LICQ (gradients of active constraints are linearly independent), Slater (for convex problems).

### Worked example 1.10 — a nonlinear inequality

Minimize \(f(x,y) = x^2+y^2\) subject to \(g(x,y) = x^2 + 2y^2 - 4 \le 0\).

The unconstrained minimizer \((0,0)\) satisfies \(g(0,0) = -4 < 0\) (strictly feasible), so \(x^\star=(0,0)\), \(\mu^\star = 0\) — inactive constraint. Optimal value \(0\). This illustrates complementary slackness: **when the constraint does not bite, the unconstrained solution is the answer**.

Now modify: minimize \(x^2+y^2\) subject to \(x^2+2y^2 - 1\ge 0\) (i.e. \(-x^2 - 2y^2 + 1 \le 0\)).
The unconstrained minimizer \((0,0)\) is infeasible. The active-constraint case forces \(x^2+2y^2=1\). By symmetry one can show \(y^\star = 0\), \(x^\star = \pm 1\) with \(\mu=1\), \(f^\star = 1\). (Exercise: verify by writing the KKT system.)

---

## 1.9 The line-search subproblem

Many algorithms (steepest descent, quasi-Newton, conjugate gradient) update

\[
x_{k+1} = x_k + \alpha_k d_k
\]

where \(d_k\) is a *descent direction* (\(\nabla f(x_k)^\top d_k < 0\)) and \(\alpha_k\) is a **step length** chosen by solving the scalar subproblem

\[
\min_{\alpha>0}\ \phi(\alpha) := f(x_k + \alpha d_k).
\]

### 1.9.1 Bracketing

First find three points \(\alpha^{(1)} < \alpha^{(2)} < \alpha^{(3)}\) with
\(\phi(\alpha^{(2)}) < \phi(\alpha^{(1)})\) and \(\phi(\alpha^{(2)}) < \phi(\alpha^{(3)})\). This bracket contains a local minimum.

Algorithm (bracket by expansion):

```
α₀ = 0, φ₀ = φ(α₀)
α₁ = s, φ₁ = φ(α₁)            # s = user-supplied initial step
if φ₁ > φ₀: swap and negate s  # make sure we descend
α₂ = α₁ + s, φ₂ = φ(α₂)
while φ₂ < φ₁:
    α₀, α₁, α₂ = α₁, α₂, α₂ + 2·s
    φ₀, φ₁, φ₂ = φ₁, φ₂, φ(α₂)
return (α₀, α₁, α₂)
```

### 1.9.2 Sectioning — golden-section search

Once a bracket \([a,b]\) is known, we shrink it. **Golden-section** places two interior points

\[
c = b - \tau(b-a),\qquad d = a + \tau(b-a),\qquad \tau = \frac{\sqrt 5 - 1}{2} \approx 0.618,
\]

evaluates \(\phi\) there, and keeps the subinterval that contains the smaller value. Each iteration reduces the interval by a factor \(\tau \approx 0.618\) **with one new evaluation** — the interval length \(L_k\) satisfies \(L_k \approx \tau^k L_0\), so the sub-interval reaching a tolerance \(\epsilon\) takes \(k \approx \log(\epsilon/L_0)/\log\tau\) iterations.

### Worked example 1.11 — golden-section on a smooth quadratic

Minimize \(\phi(\alpha) = (\alpha - 2)^2 + 1\) on \([0,5]\) using golden-section search (\(\tau=0.618\)).

| Iter | Interval | \(c\) | \(d\) | \(\phi(c)\) | \(\phi(d)\) | Keep |
|:---:|:---:|:---:|:---:|:---:|:---:|---|
| 1 | \([0, 5]\) | 1.910 | 3.090 | 1.008 | 2.188 | \(b \leftarrow d\), \([0, 3.09]\) |
| 2 | \([0, 3.09]\) | 1.180 | 1.910 | 1.672 | 1.008 | \(a \leftarrow c\), \([1.18, 3.09]\) |
| 3 | \([1.18, 3.09]\) | 1.910 | 2.360 | 1.008 | 1.130 | \(b \leftarrow d\), \([1.18, 2.36]\) |
| 4 | \([1.18, 2.36]\) | 1.630 | 1.910 | 1.136 | 1.008 | \(a \leftarrow c\), \([1.63, 2.36]\) |
| 5 | \([1.63, 2.36]\) | 1.910 | 2.081 | 1.008 | 1.007 | \(\approx [1.91, 2.08]\) |

After 5 iterations the bracket contains the exact minimum \(\alpha^\star = 2\) to within \(\sim 0.1\), having used just \(5+1=6\) evaluations of \(\phi\) (the reuse of the "carried" interior point is what makes golden-section efficient).

### 1.9.3 Convergence criteria for line search

The two classical conditions on \(\alpha_k\) are:

- **Armijo (sufficient decrease):** \(\phi(\alpha_k) \le \phi(0) + c_1 \alpha_k \phi'(0)\), with \(c_1 \approx 10^{-4}\).
- **Wolfe (curvature):** \(|\phi'(\alpha_k)| \le c_2 |\phi'(0)|\), with \(0<c_1<c_2<1\).

They together guarantee that the line-search subproblem terminates in finitely many iterations.

---

## 1.10 Derivatives of the objective function

### 1.10.1 Numerically derived gradient

Finite-difference formulas:

| Formula | Expression | Error |
|---|---|---|
| Forward | \(\dfrac{f(x+he_i) - f(x)}{h}\) | \(O(h)\) |
| Backward | \(\dfrac{f(x) - f(x-he_i)}{h}\) | \(O(h)\) |
| Central | \(\dfrac{f(x+he_i) - f(x-he_i)}{2h}\) | \(O(h^2)\) |

With floating-point arithmetic the truncation error \(O(h)\) and round-off error \(O(\epsilon_{\text{mach}}/h)\) balance; the optimal choice is roughly \(h \sim \sqrt{\epsilon_{\text{mach}}} \approx 10^{-8}\) for forward differences, and \(h \sim \epsilon_{\text{mach}}^{1/3} \approx 10^{-5}\) for central differences.

Cost: computing \(\nabla f\) by finite differences costs \(n+1\) (forward) or \(2n\) (central) objective evaluations — prohibitively expensive when \(n\) is large.

### Worked example 1.12 — compare numeric and analytic gradients

For \(f(x_1,x_2) = x_1^2 x_2 + x_2^3 - 2x_1\) at \(x = (1,2)\):

- Analytic: \(\nabla f = (2x_1 x_2 - 2,\ x_1^2 + 3x_2^2) = (2, 13)\).
- Central-difference with \(h = 10^{-6}\) gives, to seven significant figures, \((2.000000,\ 13.000000)\). Forward-difference gives about \((2.000001,\ 13.000001)\) — six correct digits.

**Rule of thumb.** Use central differences when \(n\) is small and accuracy matters; use forward differences inside large-scale optimizers where \(n\) is large and expensive.

### 1.10.2 Analytically derived gradient

Preferred whenever available:

- exact (no truncation error),
- typically an order of magnitude faster,
- required for *guaranteed* convergence of Newton-type methods.

Use **automatic differentiation** (AD) when the model is a big computational graph — tools such as JAX, PyTorch, TensorFlow, Stan, and CasADi give machine-precision derivatives at the cost of one forward and one reverse pass.

### 1.10.3 Adjoint variable approach

Adjoint methods compute the *gradient with respect to many parameters at the cost of one (or two) forward+adjoint solves*, regardless of how many parameters there are. This is the reason they dominate optimal-control and inverse problems.

Consider the discrete-time system

\[
x_{k+1} = F_k(x_k, u_k),\qquad x_0 \text{ given},\qquad
J(u) = \sum_{k=0}^{N-1} L_k(x_k, u_k) + \Phi(x_N).
\]

Introduce **adjoint variables** (costates) \(\lambda_k\) defined by the backward recursion

\[
\lambda_N = \nabla \Phi(x_N),\qquad
\lambda_k = \nabla_x L_k(x_k,u_k) + \Bigl(\frac{\partial F_k}{\partial x_k}\Bigr)^\top \lambda_{k+1}.
\]

Then the gradient w.r.t. each control is

\[
\boxed{\;\nabla_{u_k} J = \nabla_{u_k} L_k(x_k,u_k) + \Bigl(\frac{\partial F_k}{\partial u_k}\Bigr)^\top \lambda_{k+1}\;}
\]

The forward pass produces \(x_0,\dots,x_N\); the backward pass produces \(\lambda_N,\dots,\lambda_0\). Total cost: **one forward + one backward pass**, independent of \(\dim u\).

### Worked example 1.13 — adjoint derivation verified by hand

Discrete problem with \(x_0 = 2\), \(x_{k+1} = x_k + u_k\) (\(k=0,1\)), and
\[
J(u_0,u_1) = x_2^2 + u_0^2 + u_1^2.
\]

**Direct elimination** (to have a reference answer):

\[
x_1 = 2 + u_0,\qquad x_2 = 2 + u_0 + u_1,\qquad
J = (2+u_0+u_1)^2 + u_0^2 + u_1^2.
\]

\[
\frac{\partial J}{\partial u_0} = 2(2+u_0+u_1) + 2u_0,\qquad
\frac{\partial J}{\partial u_1} = 2(2+u_0+u_1) + 2u_1.
\]

Setting both to zero: \(4 + 4u_0 + 2u_1 = 0\), \(4 + 2u_0 + 4u_1 = 0\); subtracting gives \(u_0 = u_1\), hence \(u_0 = u_1 = -2/3\).

**Adjoint method.**

- Terminal: \(\lambda_2 = 2 x_2 = 2(2+u_0+u_1)\).
- Recurrence: \(\lambda_k = \partial L_k/\partial x_k + (\partial F_k/\partial x_k)\lambda_{k+1}\).
  Since \(L_0 = u_0^2\) and \(L_1 = u_1^2\) do not depend on \(x\), we get
  \(\lambda_1 = \lambda_2\) and \(\lambda_0 = \lambda_1 = \lambda_2\).
- Gradients:
  \[
  \frac{\partial J}{\partial u_0} = 2u_0 + \lambda_1,\qquad
  \frac{\partial J}{\partial u_1} = 2u_1 + \lambda_2.
  \]
- At \(u_0 = u_1 = -2/3\): \(x_2 = 2/3\), \(\lambda_1=\lambda_2 = 4/3\), so
  \(\partial J/\partial u_0 = -4/3 + 4/3 = 0\) and \(\partial J/\partial u_1 = 0\) — matching the direct computation.

!!! tip "Continuous-time version"
    For \(\dot x = F(x,u)\), \(J = \int_0^T L(x,u)\,dt + \Phi(x(T))\), the continuous adjoint is
    \[
      -\dot\lambda = \nabla_x H,\quad \lambda(T) = \nabla \Phi(x(T)),\quad H(x,u,\lambda) = L + \lambda^\top F,
    \]
    and \(\nabla_u J = \nabla_u H\). Use a forward ODE solve for \(x(t)\) and a backward ODE solve for \(\lambda(t)\).

---

## 1.11 Computational implementation

The snippets below implement *gradient checking* (comparing analytic and numeric derivatives) and *golden-section search* in three languages. You can drop them straight into a Jupyter notebook, a Java console app, or a C# program.

### Python (NumPy)

```python
import numpy as np

# ---- objective and analytic derivatives -------------------------------------
def f(x):                                # f(x,y) = x^2 y + y^3 - 2x
    return x[0]**2 * x[1] + x[1]**3 - 2 * x[0]

def grad_analytic(x):
    return np.array([2*x[0]*x[1] - 2, x[0]**2 + 3*x[1]**2])

def grad_forward(f, x, h=1e-6):
    g = np.zeros_like(x, dtype=float)
    for i in range(len(x)):
        e = np.zeros_like(x, dtype=float); e[i] = h
        g[i] = (f(x + e) - f(x)) / h
    return g

def grad_central(f, x, h=1e-6):
    g = np.zeros_like(x, dtype=float)
    for i in range(len(x)):
        e = np.zeros_like(x, dtype=float); e[i] = h
        g[i] = (f(x + e) - f(x - e)) / (2 * h)
    return g

def hessian_central(f, x, h=1e-4):
    n = len(x); H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ei = np.zeros(n); ei[i] = h
            ej = np.zeros(n); ej[j] = h
            H[i, j] = (f(x+ei+ej) - f(x+ei-ej) - f(x-ei+ej) + f(x-ei-ej)) / (4*h*h)
    return H

# ---- golden-section search --------------------------------------------------
def golden_section(phi, a, b, tol=1e-6, maxit=200):
    tau = (np.sqrt(5) - 1) / 2
    c, d = b - tau*(b-a), a + tau*(b-a)
    fc, fd = phi(c), phi(d)
    for _ in range(maxit):
        if b - a < tol:
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - tau*(b-a); fc = phi(c)
        else:
            a, c, fc = c, d, fd
            d = a + tau*(b-a); fd = phi(d)
    return 0.5 * (a + b)

# ---- demo -------------------------------------------------------------------
if __name__ == "__main__":
    x = np.array([1.0, 2.0])
    print("analytic :", grad_analytic(x))
    print("forward  :", grad_forward(f, x))
    print("central  :", grad_central(f, x))
    print("Hessian  :\n", hessian_central(f, x))
    print("golden   :", golden_section(lambda a: (a-2)**2 + 1, 0.0, 5.0))
```

**Typical output:**

```
analytic : [2. 13.]
forward  : [ 2.000001 13.000001]
central  : [ 2. 13.]
Hessian  :
 [[ 4.  2.]
  [ 2. 12.]]
golden   : 1.999999...
```

### Java

```java
public class GradientCheck {

    static double f(double[] x) {
        return x[0]*x[0]*x[1] + Math.pow(x[1], 3) - 2*x[0];
    }

    static double[] gradAnalytic(double[] x) {
        return new double[]{ 2*x[0]*x[1] - 2, x[0]*x[0] + 3*x[1]*x[1] };
    }

    static double[] gradCentral(double[] x, double h) {
        double[] g = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            double[] xp = x.clone(); xp[i] += h;
            double[] xm = x.clone(); xm[i] -= h;
            g[i] = (f(xp) - f(xm)) / (2*h);
        }
        return g;
    }

    static double goldenSection(java.util.function.DoubleUnaryOperator phi,
                                double a, double b, double tol) {
        double tau = (Math.sqrt(5) - 1) / 2;
        double c = b - tau*(b-a), d = a + tau*(b-a);
        double fc = phi.applyAsDouble(c), fd = phi.applyAsDouble(d);
        while (b - a > tol) {
            if (fc < fd) { b=d; d=c; fd=fc; c=b - tau*(b-a); fc=phi.applyAsDouble(c); }
            else        { a=c; c=d; fc=fd; d=a + tau*(b-a); fd=phi.applyAsDouble(d); }
        }
        return 0.5*(a+b);
    }

    public static void main(String[] args) {
        double[] x = { 1.0, 2.0 };
        System.out.println(java.util.Arrays.toString(gradAnalytic(x)));
        System.out.println(java.util.Arrays.toString(gradCentral(x, 1e-6)));
        System.out.println(goldenSection(a -> (a-2)*(a-2) + 1, 0.0, 5.0, 1e-6));
    }
}
```

### C#

```csharp
using System;

class GradientCheck {

    static double F(double[] x) => x[0]*x[0]*x[1] + Math.Pow(x[1], 3) - 2*x[0];

    static double[] GradAnalytic(double[] x) =>
        new[]{ 2*x[0]*x[1] - 2, x[0]*x[0] + 3*x[1]*x[1] };

    static double[] GradCentral(double[] x, double h) {
        var g = new double[x.Length];
        for (int i = 0; i < x.Length; i++) {
            var xp = (double[])x.Clone(); xp[i] += h;
            var xm = (double[])x.Clone(); xm[i] -= h;
            g[i] = (F(xp) - F(xm)) / (2*h);
        }
        return g;
    }

    static double GoldenSection(Func<double,double> phi, double a, double b, double tol) {
        double tau = (Math.Sqrt(5) - 1) / 2;
        double c = b - tau*(b-a), d = a + tau*(b-a);
        double fc = phi(c), fd = phi(d);
        while (b - a > tol) {
            if (fc < fd) { b=d; d=c; fd=fc; c=b - tau*(b-a); fc=phi(c); }
            else         { a=c; c=d; fc=fd; d=a + tau*(b-a); fd=phi(d); }
        }
        return 0.5*(a+b);
    }

    static void Main() {
        var x = new double[]{ 1.0, 2.0 };
        Console.WriteLine(string.Join(", ", GradAnalytic(x)));
        Console.WriteLine(string.Join(", ", GradCentral(x, 1e-6)));
        Console.WriteLine(GoldenSection(a => (a-2)*(a-2) + 1, 0.0, 5.0, 1e-6));
    }
}
```

!!! tip "Which language for what?"
    - **Python (NumPy + SciPy):** day-to-day optimization prototyping, `scipy.optimize.check_grad`, plotting.
    - **Java:** production-style implementation inside enterprise systems; also for teaching OOP-based algorithm design.
    - **C#:** Windows/.NET infrastructure, Unity simulations, and interop with MATLAB `Coder`-generated code.
    - **C/C++:** required when you must link the optimizer to a high-fidelity simulator.
    
    All three snippets above are mathematically equivalent — a good exercise is to *run all three* and compare the printed gradient vectors to machine precision.

### Cost comparison at a glance

| Derivative method | Cost to obtain \(\nabla f\) | Accuracy | Requires user code |
|---|---|---|---|
| Analytic (hand) | 1 evaluation | exact | Yes — derivation |
| Automatic diff. | ≤ 4× cost of \(f\) | machine precision | No (framework) |
| Central differences | \(2n\) evaluations of \(f\) | \(O(h^2)\) | No |
| Forward differences | \(n+1\) evaluations | \(O(h)\) | No |
| Adjoint (discrete or continuous) | 1 forward + 1 backward | exact (to ODE error) | Yes — adjoint code |

---

## 1.12 Chapter summary

1. **Local vs global vs Pareto.** Three independent notions; in engineering you usually need at least two of them (a Pareto front of local solutions, then a global choice).
2. **Gradient = slope, Hessian = curvature.** These two objects drive every algorithm from Chapter 2 onward.
3. **Linearity & quadraticity.** Linear problems have no interior minimum; quadratic problems are *solved in closed form* — the reference for what an iterative method should converge to.
4. **Optimality conditions.** FONC \(\nabla f = 0\), SOSC \(\nabla^2 f \succ 0\); with constraints, replace by **Lagrange/KKT** with multipliers and complementary slackness.
5. **Line search.** Bracket by expansion, refine by golden-section (or bisection/interpolation); a superlinear method with a poor line search can be worse than gradient descent with a good one.
6. **Derivatives.** Prefer analytic/AD; fall back to finite differences; use adjoints when the parameter dimension is large.

---

## 1.13 Quiz

### A. Multiple choice (choose the single best answer)

1. Which of the following is **necessary** for a differentiable \(f\) to have a local minimum at an interior point \(x^\star\)?
   (a) \(f(x^\star)=0\) (b) \(\nabla f(x^\star)=0\) (c) \(\nabla^2 f(x^\star)\succ 0\) (d) \(f\) is convex

2. If \(\nabla^2 f(x^\star) \succ 0\) and \(\nabla f(x^\star) = 0\), then \(x^\star\) is:
   (a) a strict local maximum (b) a strict local minimum (c) a saddle (d) undetermined

3. Which \(x\) is **Pareto-optimal** for \(\min (x^2, (x-2)^2)\) over \(\mathbb{R}\)?
   (a) \(x=-1\) (b) \(x=1\) (c) \(x=3\) (d) \(x=0.5\) only if \(x\le0\)

4. For \(f(x) = \tfrac12 x^\top Q x - c^\top x\) with \(Q\succ0\), the unique global minimizer is:
   (a) \(Qc\) (b) \(Q^{-1}c\) (c) \(-Q^{-1}c\) (d) \(c-Q\)

5. The gradient of a linear function \(f(x)=a^\top x + b\) is:
   (a) \(a\) (b) \(b\) (c) \(0\) (d) \(a x\)

6. Which finite-difference formula has \(O(h^2)\) truncation error?
   (a) forward (b) backward (c) central (d) all three

7. In golden-section search, the interval reduction factor per iteration is:
   (a) \(1/2\) (b) \((\sqrt 5 - 1)/2 \approx 0.618\) (c) \(1/\sqrt{2}\) (d) \(1/e\)

8. For \(\min x_1^2+x_2^2\) s.t. \(x_1+x_2=2\), the solution \(x^\star\) is:
   (a) \((2,0)\) (b) \((0,2)\) (c) \((1,1)\) (d) \((\sqrt2,\sqrt2)\)

9. In the KKT conditions for \(\min f\) s.t. \(g_i(x)\le 0\), the multipliers satisfy:
   (a) \(\mu_i \le 0\) (b) \(\mu_i = 0\) (c) \(\mu_i \ge 0\) (d) no sign restriction

10. For a quadratic with condition number \(\kappa(Q)\gg 1\), an iterative method will typically:
    (a) converge faster (b) converge more slowly (c) converge in one step (d) fail

11. Which derivative method is most efficient when there are **many** parameters and one scalar output?
    (a) forward finite differences (b) central finite differences (c) reverse-mode AD / adjoint (d) forward-mode AD

12. If \(\phi(\alpha) = f(x + \alpha d)\) with \(\nabla f(x)^\top d > 0\), then \(d\) is:
    (a) a descent direction (b) an ascent direction (c) orthogonal to \(\nabla f\) (d) not enough information

### B. Fill in the blank

1. The first-order necessary condition for an unconstrained minimizer is \(\nabla f(x^\star) = \_\_\_\_\).
2. The Hessian is the matrix of \_\_\_\_\_\_\_\_ partial derivatives.
3. A point that no other feasible point Pareto-dominates is called \_\_\_\_\_\_\_\_.
4. The second-order Taylor term is \(\tfrac12 d^\top \_\_\_\_ d\).
5. The central-difference formula for \(\partial f/\partial x_i\) is \((f(x+h e_i) - f(x-h e_i)) / \_\_\_\_\).
6. Golden-section uses \(\tau = (\sqrt 5 - 1)/2 \approx \_\_\_\_\).
7. For equality constraints \(Ax=b\), the FONC is \(\nabla f(x^\star) + \_\_\_\_ = 0\).
8. KKT complementary slackness reads \(\mu_i g_i(x^\star) = \_\_\_\_\).
9. The adjoint (costate) method's gradient cost is independent of the number of \_\_\_\_\_.
10. A symmetric matrix \(Q\) is positive definite if all eigenvalues are \_\_\_\_.

### C. Define / explain

1. Define **local minimum**, **strict local minimum**, and **global minimum**. Give a one-variable example illustrating the difference between local and global.
2. Explain **Pareto optimality** and the **Pareto front** with an example. Why is there usually no single optimum for a multi-objective problem?
3. State the **first-order** and **second-order** optimality conditions (both necessary and sufficient) for an unconstrained smooth problem. Give a counterexample showing that FONC alone is not sufficient.
4. Describe the **adjoint variable approach** for gradient computation. In two or three sentences, explain why it is preferable to finite differences for problems with many parameters.
5. Explain **bracketing** and **sectioning** in a line search. What is the reduction factor of golden-section, and what does the name mean?
6. For a quadratic form \(\tfrac12 x^\top Qx - c^\top x\), explain how the **eigenvalues of \(Q\)** determine (i) whether a unique minimum exists, and (ii) the shape of level sets. Give the formula for the optimal value.

---

## 1.14 Answer key

??? success "A. Multiple choice"
    1. **(b)** — \(\nabla f(x^\star)=0\) is necessary but not sufficient; \(\nabla^2 f\succ0\) is only sufficient.
    2. **(b)** — FONC + strict SOSC gives a strict local minimum.
    3. **(b)** — In the Pareto set \([0,2]\); the rest are dominated (see Example 1.4).
    4. **(b)** — \(Qx^\star = c\), so \(x^\star = Q^{-1}c\).
    5. **(a)** — \(\nabla (a^\top x+b) = a\), a constant.
    6. **(c)** — central differences are \(O(h^2)\); forward/backward are \(O(h)\).
    7. **(b)** — \(1/\varphi = (\sqrt 5-1)/2\approx 0.618\).
    8. **(c)** — \(x_1=x_2=1\) from symmetry and the constraint (Example 1.8).
    9. **(c)** — \(\mu_i\ge0\) for minimisation with \(g_i\le0\).
    10. **(b)** — large \(\kappa\) = slow convergence; \(\kappa=1\) would be ideal.
    11. **(c)** — reverse-mode AD / adjoint: one forward + one reverse pass.
    12. **(b)** — positive directional derivative means \(\phi\) initially increases.

??? success "B. Fill in the blank"
    1. \(0\)
    2. second (second-order / mixed)
    3. Pareto-optimal (Pareto-efficient)
    4. \(\nabla^2 f(x)\) (i.e. the Hessian)
    5. \(2h\)
    6. \(0.618\) (golden ratio conjugate)
    7. \(A^\top\lambda\) (or \(A^\top \lambda^\star\))
    8. \(0\)
    9. parameters (unknowns / control variables / inputs)
    10. positive (strictly greater than zero)

??? success "C. Define / explain — model answers"
    1. *Local minimum:* \(f(x^\star)\le f(x)\) for all feasible \(x\) in some neighborhood of \(x^\star\).
       *Strict local minimum:* the inequality is strict for \(x\ne x^\star\).
       *Global minimum:* the inequality holds over the whole feasible set.
       *Example:* \(f(x)=(x^2-1)^2\) has strict local minima at \(x=\pm1\) and both are also global; by contrast \(f(x)=0.5 x^2+3\sin(3x)\) on \([-3,3]\) has multiple local minima but only one is global.
    2. *Pareto dominance:* \(x\) dominates \(y\) if \(f_i(x)\le f_i(y)\) for all \(i\) and strictly for at least one \(i\). *Pareto-optimal:* not dominated by any feasible point. The *Pareto front* is the image of the Pareto set in objective space. There is generally no single optimum because objectives conflict — improving one worsens another. Example: \(\min(x^2,(x-2)^2)\); the Pareto set is \([0,2]\).
    3. FONC: \(\nabla f(x^\star)=0\). SOSC-necessary: \(\nabla^2 f(x^\star)\succeq0\). SOSC-sufficient: \(\nabla^2 f(x^\star)\succ0\) ⇒ strict local min. Counterexample to FONC alone: \(f(x)=x^3\) at \(x=0\) — \(\nabla f(0)=0\) but \(0\) is an inflection, not a minimum.
    4. The adjoint (costate) method augments the state with a dual variable \(\lambda\) obeying a backward-in-time (or backward-in-index) equation; a single forward pass yields the states and a single backward pass yields the costates. The gradient with respect to *all* parameters is then obtained by simple products, at a cost independent of the parameter dimension. This contrasts with finite differences, whose cost grows linearly with the number of parameters.
    5. *Bracketing* is the phase that finds three points \(a<b<c\) with \(\phi(b)<\phi(a)\) and \(\phi(b)<\phi(c)\) — a *bracket* that must contain a local minimum. *Sectioning* refines the bracket by placing one or two interior test points and discarding the region that cannot contain the minimum. Golden-section places the interior points at a ratio \(\tau=(\sqrt5-1)/2\approx0.618\) from each end; each iteration reduces the interval to a fraction \(\tau\) of its length while reusing one previously computed value — this is why it is called *golden-section*.
    6. (i) A unique global minimum of \(\tfrac12 x^\top Qx - c^\top x\) exists **iff \(Q\succ0\)** (all eigenvalues positive); if \(Q\) is indefinite, the quadratic is unbounded below along a negative-curvature direction. (ii) Level sets are ellipsoids whose axes are eigenvectors of \(Q\) with semi-axis lengths proportional to \(1/\sqrt{\lambda_i}\); the ratio of the longest to the shortest semi-axis is \(\sqrt{\kappa(Q)}\). Optimal value: \(f^\star = -\tfrac12 c^\top Q^{-1} c\).

---

## Further reading

- Boyd & Vandenberghe, *Convex Optimization*, §1–§5.
- Nocedal & Wright, *Numerical Optimization*, Chapters 1–3 and 10 (line search).
- Luenberger & Ye, *Linear and Nonlinear Programming*, Chapters 1–4.
- Press, Teukolsky, Vetterling & Flannery, *Numerical Recipes*, §10.
- Betts, *Practical Methods for Optimal Control and Estimation Using Nonlinear Programming*, Chapter 1.
- SciPy documentation: `scipy.optimize.check_grad`, `scipy.optimize.line_search`, `scipy.optimize.minimize`.
- CasADi / JAX tutorials on reverse-mode AD for gradient computation.

---

## Suggested tutorial problems for this chapter

1. **Unconstrained.** For \(f(x)=x^3-6x^2+9x+1\), find all stationary points, classify them, and identify the global minimum on \([0,4]\). Compare with a numerical golden-section result.
2. **Pareto.** Compute and sketch the Pareto front of \(\min(x^2+y^2, (x-1)^2+y^2)\) over \(\mathbb{R}^2\). Identify at least three Pareto-optimal points.
3. **Lagrange.** Use the Lagrange multiplier method to minimize \(x_1^2 + x_2^2 + x_3^2\) subject to \(x_1 + 2x_2 + 3x_3 = 6\). Verify the multiplier via the envelope theorem.
4. **KKT.** Solve \(\min (x-3)^2 + (y-2)^2\) s.t. \(x+y\le4\), \(x\ge0\), \(y\ge0\). Which constraints are active at the optimum, and what are the multipliers?
5. **Line search.** Implement golden-section search and use it to minimize \(\phi(\alpha)=\sin(3\alpha)+\alpha^2/4\) on \([0,2]\). Verify convergence against `scipy.optimize.minimize_scalar` (method `'golden'`).
6. **Adjoint.** Extend Example 1.13 to \(N=4\), \(x_{k+1}=x_k + u_k\), \(J = x_4^2 + \sum_{k=0}^3 u_k^2\). Derive the optimal \(u_0,\dots,u_3\) by both direct differentiation and the adjoint recursion, and confirm they agree.
7. **Gradient check.** For \(f(x_1,x_2,x_3) = \sin(x_1 x_2) + x_3^4 - x_1 e^{x_2 x_3}\), compute the analytic gradient at \((0.5, -1, 0.3)\) and compare with central-difference and reverse-mode AD (JAX) results.
