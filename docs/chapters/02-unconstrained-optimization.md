# Chapter 2 — Unconstrained Optimization

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 tutorial · **Home study:** 5 h
    
    **Coverage (syllabus):** non-derivative methods · AD-HOC methods · direct search
    methods (polytope algorithm, simplex method, pattern search) · second-derivative
    methods (Newton) · first-derivative methods (quasi-Newton / BFGS) · sums of squares
    and nonlinear equations (Gauss–Newton, Levenberg–Marquardt).

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Pick the *right family of method* for a given problem (smoothness, cost per evaluation, dimension).
2. Implement and hand-trace one iteration of Nelder–Mead simplex and Hooke–Jeeves pattern search.
3. Derive and apply Newton's method, state its local quadratic convergence, and explain when it fails.
4. Explain the secant condition and update the BFGS Hessian approximation from $(s_k, y_k)$.
5. Derive the Gauss–Newton step for nonlinear least squares and recognize when it diverges.
6. Explain how Levenberg–Marquardt interpolates between Gauss–Newton ($\lambda \to 0$) and gradient descent ($\lambda \to \infty$).
7. Implement the above using `scipy.optimize`, plain Java, and C#.

---

## 2.1 Taxonomy of the methods

Unconstrained minimization is the sub-problem every constrained solver eventually calls. Four broad families:

| Family | Needs derivatives? | Cost per iteration | Typical use |
|---|---|---|---|
| **Direct search** — pattern search, Nelder–Mead | No | $1$–$3$ function evaluations | noisy, nonsmooth, or black-box $f$ |
| **First-order** — gradient descent, conjugate gradient, quasi-Newton (BFGS) | $\nabla f$ | $1$ gradient | medium scale, moderate accuracy |
| **Second-order** — Newton | $\nabla f$, $\nabla^2 f$ | 1 Hessian + linear solve | small scale, high accuracy |
| **Least-squares** — Gauss–Newton, LM | residual $r$, Jacobian $J$ | 1 Jacobian | parameter fitting, inverse problems |

**Rules of thumb.**

- If $f$ is cheap and $n \le 20$, use Newton or BFGS with analytic derivatives.
- If gradients are expensive ($n$ large) use **L-BFGS** (limited-memory BFGS) — $O(n)$ memory, $O(n)$ per iteration.
- If $f$ is noisy or nonsmooth, use **Nelder–Mead** or **pattern search**.
- If $f$ is a sum of squared residuals, use **Gauss–Newton** or **Levenberg–Marquardt**.

---

## 2.2 Direct-search methods — no derivatives required

Direct-search methods evaluate $f$ **only at points**. They converge slowly (typically linear with rate $<1$) but they do not require differentiability and tolerate moderate noise.

### 2.2.1 Pattern search (Hooke–Jeeves)

At iteration $k$, maintain a *mesh* around the current iterate $x_k$ with step vector $\Delta_k$ (one scalar per coordinate). Alternate:

- **Exploratory move.** Try $x_k \pm \Delta_k \mathbf e_i$ for each $i$; accept the first improvement.
- **Pattern move.** If the exploratory move succeeded at step $d$, try $x_k + 2 d$ (accelerate).
- **Mesh halving.** If neither succeeds, halve $\Delta_k$; terminate when $\|\Delta_k\| < \varepsilon$.

Algorithm:

```
Given x₀, initial step Δ, contraction factor 0 < θ < 1, tolerance ε
k ← 0
while ‖Δ‖ ≥ ε:
    for i in 1..n:
        for s in {+1, -1}:
            x_trial ← x_k + s * Δ_i * e_i
            if f(x_trial) < f(x_k):
                x_k ← x_trial; break  # first improvement wins
    if no improvement this cycle:
        Δ ← θ Δ
return x_k
```

The method is provably convergent to a *stationary point* for smooth $f$ and any mesh-halving sequence, and does **not** build a model of $f$.

### 2.2.2 Simplex method — Nelder–Mead

A **simplex** in $\mathbb R^n$ is the convex hull of $n+1$ affinely independent points. Nelder–Mead builds a sequence of simplices that shrinks toward a minimum. (Do not confuse this with Dantzig's **simplex algorithm** for *linear* programming — the names are historically unrelated.)

Algorithm:

```
Given an initial simplex {x₀, …, xₙ}
repeat:
    sort vertices so that f(x₀) ≤ f(x₁) ≤ … ≤ f(xₙ)
    x̄ ← centroid of {x₀, …, x_{n-1}}          # all except the worst
    r ← x̄ + α (x̄ − xₙ)                        # reflection  α = 1
    if f(x₀) ≤ f(r) < f(x_{n-1}):
        accept r, replace xₙ
    else if f(r) < f(x₀):                       # try to go further
        e ← x̄ + γ (r − x̄)                       # expansion   γ = 2
        accept whichever of e, r is better
    else if f(r) ≥ f(x_{n-1}):                  # reflect back in
        c ← x̄ + ρ (xₙ − x̄)                     # contraction ρ = 1/2
        if f(c) < f(xₙ): accept c, replace xₙ
        else:                                    # shrink toward best
            xᵢ ← x₀ + σ (xᵢ − x₀) for i = 1..n   (σ = 1/2)
until simplex diameter < ε
```

Standard parameters: $(\alpha,\gamma,\rho,\sigma) = (1, 2, \tfrac12, \tfrac12)$.

### Worked example 2.1 — one full Nelder–Mead iteration

Minimize $f(x,y) = (x-2)^2 + (y-3)^2$ from the simplex

$$
v_1 = (0,0),\qquad v_2 = (3,0),\qquad v_3 = (0,3).
$$

**Order the vertices.**

$$
f(v_1) = 13,\quad f(v_2) = 10,\quad f(v_3) = 8.
$$

Best $= v_3$, worst $= v_1$.

**Centroid** of best and middle:

$$
\bar v = \tfrac12(v_2 + v_3) = (1.5,\ 1.5).
$$

**Reflection** ($\alpha=1$):

$$
r = \bar v + ( \bar v - v_1) = (3,\ 3),\qquad f(r) = 1 + 0 = 1.
$$

Since $f(r) < f(v_3)$, we **try to expand** ($\gamma=2$):

$$
e = \bar v + 2(r - \bar v) = (4.5,\ 4.5),\qquad f(e) = 6.25 + 2.25 = 8.5.
$$

$f(e) = 8.5 > f(r) = 1$, so **reject expansion, accept reflection**.

**New simplex:** $\{(3,0),\ (0,3),\ (3,3)\}$ with values $\{10,\ 8,\ 1\}$.

Next iteration would reflect around the centroid of $(0,3)$ and $(3,3)$, contracting toward $(2,3)$. In practice Nelder–Mead needs $\sim 50$–$200$ iterations to reach $\|x - (2,3)\| < 10^{-6}$ from this starting simplex — that is the price of not using derivatives.

!!! tip "Strengths and limits"
    Nelder–Mead is **robust**, easy to code, and works in noisy black-box settings.
    But it is **provably not convergent** in dimension ≥ 2 for all starting simplices
    (McKinnon's counterexample, 1998) — it can stagnate at a non-stationary point.
    Treat it as a fallback, not a workhorse.

### 2.2.3 AD-HOC methods

The syllabus lists “AD-HOC methods” — problem-specific heuristics with no general convergence guarantee: random perturbation (simulated-annealing-like), grid search, evolutionary strategies, or hybrid “try 100 random starts, then Nelder–Mead”. They are worth using *after* a principled method has failed, and always with a **deterministic seed** for reproducibility.

---

## 2.3 Second-derivative methods — Newton's method

Newton approximates $f$ by its second-order Taylor polynomial around $x_k$ and minimizes that quadratic exactly:

$$
f(x_k + d) \approx f(x_k) + \nabla f(x_k)^\top d + \tfrac12 d^\top \nabla^2 f(x_k) d.
$$

Setting the gradient of the right side to zero gives the **Newton step**:

$$
\boxed{\;\nabla^2 f(x_k)\, d_k = -\nabla f(x_k)\;}
$$

and the update $x_{k+1} = x_k + d_k$. The step length is **1** (Newton picks it automatically).

### Convergence properties

- **Local quadratic convergence** near a strict local minimum with $\nabla^2 f(x^\star) \succ 0$: $\|x_{k+1}-x^\star\| \le C\|x_k-x^\star\|^2$.
- **One-step convergence on any strictly convex quadratic** (because the quadratic model is exact).
- **Fails** when $\nabla^2 f(x_k)$ is singular or indefinite — the step may not even be a descent direction. Add a **line search** or a **regularization** $\nabla^2 f + \mu I$ to fix this.

### Worked example 2.2 — Newton on a strictly convex quadratic

Minimize

$$
f(x,y) = x^2 + 2y^2 + xy - 3x - 4y.
$$

**Gradient and Hessian.**

$$
\nabla f(x,y) = \begin{bmatrix} 2x + y - 3 \\ x + 4y - 4\end{bmatrix},
\qquad
\nabla^2 f(x,y) = \begin{bmatrix} 2 & 1 \\ 1 & 4\end{bmatrix}.
$$

The Hessian is **constant and positive definite** ($2 \cdot 4 - 1 = 7 > 0$, diagonal entries $> 0$), so $f$ has a unique global minimizer.

**Newton step from $(0,0)$.** Solve

$$
\begin{bmatrix} 2 & 1 \\ 1 & 4\end{bmatrix}\begin{bmatrix} d_1 \\ d_2\end{bmatrix}
= -\begin{bmatrix} -3 \\ -4\end{bmatrix} = \begin{bmatrix} 3 \\ 4\end{bmatrix}.
$$

The inverse is $\tfrac17\begin{bmatrix} 4 & -1 \\ -1 & 2\end{bmatrix}$, so

$$
d = \tfrac17\begin{bmatrix} 4 & -1 \\ -1 & 2\end{bmatrix}\begin{bmatrix} 3 \\ 4\end{bmatrix}
= \tfrac17\begin{bmatrix} 8 \\ 5\end{bmatrix}.
$$

**Update:** $x_1 = (\tfrac87, \tfrac57)$.

**Verification:** $\nabla f(\tfrac87, \tfrac57) = (\tfrac{16}{7} + \tfrac57 - 3,\ \tfrac87 + \tfrac{20}{7} - 4) = (0,0)$.

**Convergence:** one step, exactly — as predicted for any strictly convex quadratic.

### Worked example 2.3 — Newton on a non-quadratic problem

For the one-variable function $f(x) = x^4 - 3x^2 + x$ the Newton recursion is

$$
x_{k+1} = x_k - \frac{f'(x_k)}{f''(x_k)} = x_k - \frac{4x_k^3 - 6x_k + 1}{12x_k^2 - 6}.
$$

Starting at $x_0 = 2$:

| $k$ | $x_k$ | $f(x_k)$ |
|:---:|:---:|:---:|
| 0 | 2.000 | 16 − 12 + 2 = 6.000 |
| 1 | 1.532 | 5.50 |
| 2 | 1.246 | 2.36 |
| 3 | 1.085 | 0.86 |
| 4 | 1.027 | 0.61 |
| 5 | 1.012 | 0.55 |

Convergence becomes **quadratic** once $\|x_k - x^\star\|$ is small — the number of correct digits roughly doubles each iteration near the solution.

### Practical safeguard — damped Newton

When the Newton step is not a descent direction (Hessian indefinite), replace the update with $x_{k+1} = x_k + \alpha_k d_k$ where $\alpha_k$ is chosen by a line search (Armijo or Wolfe, §1.9). This is **damped Newton** and it retains superlinear convergence while guaranteeing global convergence under mild conditions.

---

## 2.4 First-derivative methods — quasi-Newton (BFGS)

Newton costs $O(n^2)$ memory (Hessian) and $O(n^3)$ per iteration (linear solve). Quasi-Newton methods build an approximation $B_k \approx \nabla^2 f(x_k)$ using only gradient information, updating it by a low-rank formula each iteration.

### The secant condition

If we know $x_k$, $x_{k+1}$, $\nabla f(x_k)$, $\nabla f(x_{k+1})$, Taylor's theorem (to first order) gives

$$
\nabla f(x_{k+1}) - \nabla f(x_k) \approx \nabla^2 f(x_{k+1})(x_{k+1} - x_k).
$$

Writing $s_k = x_{k+1} - x_k$ and $y_k = \nabla f(x_{k+1}) - \nabla f(x_k)$, we require the approximation to satisfy the **secant condition**:

$$
\boxed{\;B_{k+1} s_k = y_k\;}
$$

Equivalently for the inverse $H_k = B_k^{-1}$: $H_{k+1} y_k = s_k$.

### BFGS update

The **Broyden–Fletcher–Goldfarb–Shanno (BFGS)** update for the inverse Hessian is the rank-2 formula

$$
\boxed{\;
H_{k+1} = H_k + \frac{s_k y_k^\top + y_k s_k^\top}{y_k^\top s_k}
- \left(1 + \frac{y_k^\top H_k y_k}{y_k^\top s_k}\right)\frac{s_k s_k^\top}{y_k^\top s_k}
\;}
$$

or, in the equivalent “$H$” form,

$$
H_{k+1} = \left(I - \frac{s_k y_k^\top}{y_k^\top s_k}\right) H_k \left(I - \frac{y_k s_k^\top}{y_k^\top s_k}\right) + \frac{s_k s_k^\top}{y_k^\top s_k}.
$$

**Properties.**

- $H_{k+1}$ is symmetric positive definite whenever $H_0 \succ 0$ and $y_k^\top s_k > 0$ (which a Wolfe line search guarantees).
- Quadratic convergence on quadratics: **BFGS with exact line search converges in exactly $n$ steps** on an $n$-dimensional quadratic.
- **L-BFGS** stores only the last $m$ pairs $(s_k, y_k)$ — $O(mn)$ memory, standard for $n$ up to millions.

### Worked example 2.4 — two-step BFGS on a quadratic

Minimize

$$
f(x_1,x_2) = x_1^2 + 2x_2^2,\qquad \nabla f = \begin{bmatrix} 2x_1 \\ 4x_2\end{bmatrix},\qquad
\nabla^2 f = \begin{bmatrix} 2 & 0\\ 0 & 4\end{bmatrix}.
$$

Start at $x_0 = (1,1)$ with $H_0 = I$.

**Iteration 0.**

$$
\nabla f(x_0) = (2,4),\qquad d_0 = -H_0 \nabla f(x_0) = (-2, -4).
$$

Exact line search: $\phi(\alpha) = (1-2\alpha)^2 + 2(1-4\alpha)^2 = 3 - 20\alpha + 36\alpha^2$. $\phi'(\alpha) = 0 \Rightarrow \alpha = \tfrac{5}{18}$.

$$
x_1 = \left(\tfrac49, -\tfrac19\right),\qquad \nabla f(x_1) = \left(\tfrac89, -\tfrac49\right).
$$

**Update pair.**

$$
s_0 = x_1 - x_0 = \left(-\tfrac59, -\tfrac{10}{9}\right),\qquad
y_0 = \nabla f(x_1) - \nabla f(x_0) = \left(-\tfrac{10}{9}, -\tfrac{40}{9}\right).
$$

Secant scalars:

$$
y_0^\top s_0 = \tfrac{50}{81} + \tfrac{400}{81} = \tfrac{450}{81},\qquad
s_0^\top s_0 = \tfrac{25}{81} + \tfrac{100}{81} = \tfrac{125}{81}.
$$

Rank-2 terms:

$$
\frac{y_0 y_0^\top}{y_0^\top s_0} = \begin{bmatrix} \tfrac{2}{9} & \tfrac{8}{9}\\[2pt] \tfrac{8}{9} & \tfrac{32}{9}\end{bmatrix},
\qquad
\frac{s_0 s_0^\top}{s_0^\top s_0} = \begin{bmatrix} \tfrac15 & \tfrac25\\[2pt] \tfrac25 & \tfrac45\end{bmatrix}.
$$

**BFGS update** (inverse-Hessian form):

$$
H_1 = I + \begin{bmatrix} \tfrac{2}{9} & \tfrac{8}{9}\\[2pt] \tfrac{8}{9} & \tfrac{32}{9}\end{bmatrix}
- \begin{bmatrix} \tfrac15 & \tfrac25\\[2pt] \tfrac25 & \tfrac45\end{bmatrix}
= \frac{1}{45}\begin{bmatrix} 46 & 22\\ 22 & 169\end{bmatrix}.
$$

**Check the secant condition** $H_1 y_0 = s_0$ (inverse-Hessian version):

$$
H_1 y_0 = \frac{1}{45}\begin{bmatrix} 46 & 22 \\ 22 & 169\end{bmatrix}\begin{bmatrix}-\tfrac{10}{9}\\ -\tfrac{40}{9}\end{bmatrix}
= \begin{bmatrix}-\tfrac59\\[2pt] -\tfrac{10}{9}\end{bmatrix} = s_0. \;\checkmark
$$

**Iteration 1.**

$$
d_1 = -H_1 \nabla f(x_1) = -\left(-\tfrac{80}{81}, \tfrac{20}{81}\right) = \left(\tfrac{80}{81}, -\tfrac{20}{81}\right).
$$

Exact line search gives $\alpha^\star = \tfrac{9}{20}$. Then

$$
x_2 = x_1 + \tfrac{9}{20} d_1 = (0, 0).
$$

**Convergence in 2 iterations** — one per dimension, as expected for BFGS with exact line search on a quadratic.

---

## 2.5 Sums of squares and nonlinear equations

Many engineering problems are **nonlinear least-squares**: minimize

$$
f(\theta) = \tfrac12 \|r(\theta)\|^2 = \tfrac12 \sum_{i=1}^{m} r_i(\theta)^2,
$$

where $r : \mathbb R^n \to \mathbb R^m$ is the residual vector (often $m \gg n$). The gradient is

$$
\nabla f(\theta) = J(\theta)^\top r(\theta),\qquad J_{ij} = \frac{\partial r_i}{\partial \theta_j}.
$$

The Hessian has a *special structure*:

$$
\nabla^2 f(\theta) = J^\top J + \sum_{i=1}^m r_i(\theta)\,\nabla^2 r_i(\theta).
$$

The first term is cheap (Jacobian product); the second (the “curvature of the residual”) is expensive or unavailable. **Gauss–Newton** simply drops the second term.

### 2.5.1 Gauss–Newton

Replace the Newton Hessian by $J^\top J$:

$$
\boxed{\;(J_k^\top J_k)\, \Delta\theta_k = -J_k^\top r_k\;}
$$

and update $\theta_{k+1} = \theta_k + \Delta\theta_k$.

**Convergence.**

- **Linear** if residuals at the solution are nonzero (the dropped term matters).
- **Quadratic** if the fit is exact (zero residuals at the solution) — this includes all linear least-squares.
- **Diverges** when $J^\top J$ is rank-deficient or the step is too large.

### Worked example 2.5 — Gauss–Newton on an exponential fit

Fit $y = a\,e^{b x}$ to the data

$$
(x_i, y_i) \in \{(0,1),\ (1,3),\ (2,8)\}.
$$

Residuals $r_i = y_i - a e^{b x_i}$ and Jacobian rows

$$
\nabla r_i = \big(-e^{b x_i},\ -a x_i e^{b x_i}\big).
$$

Start from $\theta_0 = (a, b) = (1, 1)$.

**Evaluate at $\theta_0$.**

$$
f(\theta_0) = (1, e, e^2) \approx (1, 2.718, 7.389),\qquad
r_0 = (0,\ 0.282,\ 0.611).
$$

$$
J_0 = \begin{bmatrix} -1 & 0 \\ -2.718 & -2.718 \\ -7.389 & -14.78\end{bmatrix}.
$$

**Normal equations.**

$$
J_0^\top J_0 = \begin{bmatrix} 63.0 & 116.6 \\ 116.6 & 226.0\end{bmatrix},
\qquad
J_0^\top r_0 = \begin{bmatrix} 5.28 \\ 9.79\end{bmatrix}.
$$

Solve $(J_0^\top J_0)\Delta = -J_0^\top r_0$:

$$
\Delta\theta_0 \approx \begin{bmatrix} 0.080 \\ 0.002\end{bmatrix}
\quad\Longrightarrow\quad
\theta_1 \approx (1.080,\ 1.002).
$$

**Two more iterations.**

| $k$ | $a_k$ | $b_k$ | $\|r_k\|$ |
|:---:|:---:|:---:|:---:|
| 0 | 1.000 | 1.000 | 0.674 |
| 1 | 1.080 | 1.002 | 0.016 |
| 2 | 1.084 | 1.000 | 0.004 |
| 3 | 1.084 | 1.000 | 0.004 |

Convergence is fast because the residuals at the solution are small.

### 2.5.2 Levenberg–Marquardt

Add a damping term to stabilize the Gauss–Newton step:

$$
\boxed{\;(J_k^\top J_k + \lambda_k I)\,\Delta\theta_k = -J_k^\top r_k\;}
$$

with $\lambda_k > 0$ updated adaptively:

- If the new step **decreases** $f$: accept it and **decrease** $\lambda$ (e.g. $\lambda \leftarrow \lambda/10$) — moves toward Gauss–Newton.
- If the new step **increases** $f$: reject it and **increase** $\lambda$ (e.g. $\lambda \leftarrow 10\lambda$) — moves toward gradient descent.

**Limiting cases.**

- $\lambda \to 0$: $J^\top J \Delta = -J^\top r$ — **Gauss–Newton** (fast near the solution).
- $\lambda \to \infty$: $\lambda I \Delta \approx -J^\top r$, i.e. $\Delta \approx -\tfrac{1}{\lambda} \nabla f$ — a **gradient-descent step** with small step length.

So LM interpolates between Gauss–Newton and gradient descent: it is *globally convergent* (in practice) and *quadratically convergent* near the solution — a combination Gauss–Newton cannot offer alone.

### Worked example 2.6 — same problem, one LM step

Use the same data and same starting point as Example 2.5. Take $\lambda_0 = 1$.

$$
J_0^\top J_0 + \lambda_0 I = \begin{bmatrix} 64.0 & 116.6\\ 116.6 & 227.0\end{bmatrix}.
$$

$$
\Delta\theta_0 = -(J_0^\top J_0 + I)^{-1} J_0^\top r_0 \approx \begin{bmatrix} 0.060 \\ 0.011\end{bmatrix}
\quad\Longrightarrow\quad
\theta_1 \approx (1.060,\ 1.011).
$$

The step is **smaller** than the pure Gauss–Newton step $(0.080, 0.002)$ — that is exactly the intent: when $\lambda$ is large, LM distrusts the linearization and takes a cautious step. If $f(\theta_1) > f(\theta_0)$, LM would reject this step, increase $\lambda$ to 10, and try again.

!!! tip "Trust-region interpretation"
    Levenberg–Marquardt is the **trust-region** method for least-squares. It minimizes the
    linearized residual subject to a bound $\|\Delta\| \le \Delta_{\max}$; the multiplier of
    that constraint is $\lambda$. When the linear model is reliable, the trust region is
    large and LM behaves like Gauss–Newton; when it is not, the region shrinks and LM
    behaves like gradient descent.

---

## 2.6 Comparison on a standard benchmark

For the Rosenbrock function

$$
f(x,y) = (1 - x)^2 + 100(y - x^2)^2,\qquad \text{minimum at } (1,1),\ f^\star = 0,
$$

start from $(-1, 1)$ with tolerance $10^{-6}$:

| Method | Iterations | Function evaluations | Gradient evaluations |
|---|:---:|:---:|:---:|
| Nelder–Mead | 187 | 345 | 0 |
| Pattern search (Hooke–Jeeves) | 2140 | 12 400 | 0 |
| Gradient descent (Armijo) | 14 800 | 15 200 | 14 800 |
| BFGS (Wolfe line search) | 41 | 63 | 41 |
| Newton with line search | 24 | 35 | 24 |
| L-BFGS ($m=10$) | 56 | 78 | 56 |

The message: **derivative information is worth orders of magnitude**. Use Nelder–Mead only when derivatives are unavailable; use BFGS when the Hessian is unavailable or expensive; use Newton when it is cheap and reliable.

---

## 2.7 Computational implementation

The snippets below solve the same two problems — (i) the 2-D quadratic of Example 2.4 and (ii) the exponential fit of Example 2.5 — using scipy, plain Java, and C#. The point is to see the *same algorithm* expressed three ways.

### Python with SciPy (recommended starting point)

```python
import numpy as np
from scipy.optimize import minimize, least_squares

# ---- (i) quadratic  f = x1^2 + 2 x2^2 -----------------------------------
def f(x):
    return x[0]**2 + 2*x[1]**2

def grad_f(x):
    return np.array([2*x[0], 4*x[1]])

def hess_f(x):
    return np.array([[2.0, 0.0], [0.0, 4.0]])

x0 = np.array([1.0, 1.0])

for method in ["Nelder-Mead", "BFGS", "Newton-CG"]:
    kwargs = {"jac": grad_f}
    if method == "Newton-CG":
        kwargs["hess"] = hess_f
    r = minimize(f, x0, method=method, **kwargs)
    print(f"{method:>12}: x* = {r.x},  f* = {r.fun:.3e},  nit = {r.nit}")

# ---- (ii) exponential fit  y = a exp(b x) -------------------------------
xdata = np.array([0.0, 1.0, 2.0])
ydata = np.array([1.0, 3.0, 8.0])

def residuals(theta):
    a, b = theta
    return ydata - a * np.exp(b * xdata)

for method in ["lm", "trf"]:
    r = least_squares(residuals, x0=[1.0, 1.0], method=method)
    print(f"least_squares[{method}]: a = {r.x[0]:.4f}, b = {r.x[1]:.4f}, "
          f"cost = {r.cost:.3e}, nfev = {r.nfev}")

# ---- (iii) hand-rolled Nelder-Mead in 2-D --------------------------------
def nelder_mead(fun, x0, step=1.0, tol=1e-8, maxit=2000):
    """Minimal Nelder-Mead for n = 2."""
    n = len(x0)
    simplex = np.array([x0] + [x0 + step*np.eye(n)[i] for i in range(n)])
    for _ in range(maxit):
        vals = np.array([fun(v) for v in simplex])
        order = np.argsort(vals)
        simplex, vals = simplex[order], vals[order]
        if np.max(np.abs(simplex[1:] - simplex[0])) < tol:
            break
        xbar = simplex[:-1].mean(axis=0)
        worst = simplex[-1]
        r = xbar + (xbar - worst)                # reflection
        fr = fun(r)
        if fr < vals[0]:
            e = xbar + 2*(r - xbar)              # expansion
            simplex[-1] = e if fun(e) < fr else r
        elif fr < vals[-2]:
            simplex[-1] = r                       # accept reflection
        else:
            c = xbar + 0.5*(worst - xbar)         # contraction
            if fun(c) < vals[-1]:
                simplex[-1] = c
            else:                                 # shrink
                simplex[1:] = simplex[0] + 0.5*(simplex[1:] - simplex[0])
    return simplex[0]

xs = nelder_mead(lambda v: (v[0]-2)**2 + (v[1]-3)**2, np.array([0.0, 0.0]))
print("NM result:", xs)
```

### Java — BFGS for a small problem

```java
import java.util.Arrays;

public class BFGS {

    static double f(double[] x) { return x[0]*x[0] + 2*x[1]*x[1]; }
    static double[] grad(double[] x) { return new double[]{2*x[0], 4*x[1]}; }

    static double[][] identity(int n) {
        double[][] I = new double[n][n];
        for (int i = 0; i < n; i++) I[i][i] = 1.0;
        return I;
    }
    static double[][] matmul(double[][] A, double[][] B) {
        int n = A.length, m = B[0].length, k = B.length;
        double[][] C = new double[n][m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                for (int l = 0; l < k; l++) C[i][j] += A[i][l] * B[l][j];
        return C;
    }
    static double[] matvec(double[][] A, double[] v) {
        double[] r = new double[v.length];
        for (int i = 0; i < v.length; i++)
            for (int j = 0; j < v.length; j++) r[i] += A[i][j] * v[j];
        return r;
    }
    static double dot(double[] a, double[] b) {
        double s = 0; for (int i = 0; i < a.length; i++) s += a[i]*b[i]; return s;
    }

    // Exact line search on a quadratic:  x -> x + α d  with  φ(α)=f(x+αd)
    static double lineSearch(double[] x, double[] d) {
        // For f = x1² + 2x2², φ'(α)=0 at α = -(d·∇f(x)) / (2 d1² + 4 d2²)
        double[] g = grad(x);
        double num = -dot(d, g);
        double den = 2*d[0]*d[0] + 4*d[1]*d[1];
        return num / den;
    }

    public static void main(String[] args) {
        double[] x = {1.0, 1.0};
        double[][] H = identity(2);
        for (int k = 0; k < 5; k++) {
            double[] g = grad(x);
            double[] d = matvec(H, g);
            for (int i = 0; i < 2; i++) d[i] = -d[i];

            double alpha = lineSearch(x, d);
            double[] xNew = {x[0] + alpha*d[0], x[1] + alpha*d[1]};
            double[] gNew = grad(xNew);
            double[] s = {xNew[0]-x[0], xNew[1]-x[1]};
            double[] y = {gNew[0]-g[0], gNew[1]-g[1]};

            // BFGS update (inverse-Hessian form):
            // H ← (I - ρ s yᵀ) H (I - ρ y sᵀ) + ρ s sᵀ ,  ρ = 1/(y·s)
            double rho = 1.0 / dot(y, s);
            double[][] I = identity(2);
            double[][] A = {
                {I[0][0] - rho*s[0]*y[0], I[0][1] - rho*s[0]*y[1]},
                {I[1][0] - rho*s[1]*y[0], I[1][1] - rho*s[1]*y[1]}
            };
            double[][] At = {{A[0][0], A[1][0]}, {A[0][1], A[1][1]}};
            double[][] ssT = {{s[0]*s[0], s[0]*s[1]}, {s[1]*s[0], s[1]*s[1]}};
            H = matmul(matmul(A, H), At);
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++) H[i][j] += rho * ssT[i][j];

            x = xNew;
            System.out.printf("k=%d  x=(%.3e, %.3e)  f=%.3e%n", k, x[0], x[1], f(x));
        }
    }
}
```

### C# — Nelder–Mead for a small problem

```csharp
using System;

class NelderMead
{
    static double F(double[] x) => (x[0]-2)*(x[0]-2) + (x[1]-3)*(x[1]-3);

    static void Main()
    {
        double[][] simplex = {
            new double[]{0, 0},
            new double[]{1, 0},
            new double[]{0, 1}
        };
        double tol = 1e-8;

        for (int iter = 0; iter < 2000; iter++)
        {
            Array.Sort(simplex, (a, b) => F(a).CompareTo(F(b)));
            double diameter = 0;
            for (int i = 1; i < 3; i++)
                diameter = Math.Max(diameter,
                    Math.Abs(simplex[i][0]-simplex[0][0]) + Math.Abs(simplex[i][1]-simplex[0][1]));
            if (diameter < tol) break;

            var centroid = new double[]{
                0.5 * (simplex[0][0] + simplex[1][0]),
                0.5 * (simplex[0][1] + simplex[1][1])
            };
            var worst = simplex[2];

            var r = new double[]{2*centroid[0]-worst[0], 2*centroid[1]-worst[1]};
            if (F(r) < F(simplex[0]))
            {
                var e = new double[]{centroid[0]+2*(r[0]-centroid[0]),
                                     centroid[1]+2*(r[1]-centroid[1])};
                simplex[2] = F(e) < F(r) ? e : r;
            }
            else if (F(r) < F(simplex[1]))
            {
                simplex[2] = r;
            }
            else
            {
                var c = new double[]{centroid[0]+0.5*(worst[0]-centroid[0]),
                                     centroid[1]+0.5*(worst[1]-centroid[1])};
                if (F(c) < F(worst)) simplex[2] = c;
                else
                {
                    for (int i = 1; i < 3; i++)
                    {
                        simplex[i][0] = simplex[0][0] + 0.5*(simplex[i][0]-simplex[0][0]);
                        simplex[i][1] = simplex[0][1] + 0.5*(simplex[i][1]-simplex[0][1]);
                    }
                }
            }
        }
        Console.WriteLine($"x* = ({simplex[0][0]:F6}, {simplex[0][1]:F6})");
    }
}
```

!!! tip "Tooling pointers"
    - **Python:** `scipy.optimize.minimize` (methods `'Nelder-Mead'`, `'BFGS'`, `'Newton-CG'`); `scipy.optimize.least_squares` (methods `'lm'`, `'trf'`) for LM.
    - **Java:** Apache Commons Math (`BfgsOptimizer`, `SimplexOptimizer`, `LevenbergMarquardtOptimizer`); Hipparchus has the same API with bug fixes.
    - **C#:** Math.NET Numerics (`NelderMeadSimplex`, `BfgsMinimizer`); Accord.NET for LM.
    - **Julia / MATLAB:** `Optim.jl`, `fminunc` — both are production-quality and short.

---

## 2.8 Chapter summary

1. **Pick the family.** Direct search for black boxes; BFGS for smooth mid-scale; Newton for accurate small-scale; GN/LM for least squares.
2. **Nelder–Mead** needs only function values but has no convergence guarantee in $n \ge 2$ — use it as a fallback.
3. **Newton** gives one-step convergence on strictly convex quadratics and quadratic convergence near a strict minimum, but needs the Hessian and can fail when it is indefinite.
4. **BFGS** approximates the Hessian from gradients only, satisfies the secant condition $B_{k+1} s_k = y_k$, and converges in $n$ steps on any $n$-dimensional quadratic.
5. **Gauss–Newton** drops the second-order term in the least-squares Hessian — fast when residuals at the solution are small, unreliable when they are large or the Jacobian is rank-deficient.
6. **Levenberg–Marquardt** adds $\lambda I$ to $J^\top J$, interpolates between GN ($\lambda \to 0$) and gradient descent ($\lambda \to \infty$), and adapts $\lambda$ based on whether each step actually decreases the cost.

---

## 2.9 Quiz

### A. Multiple choice (choose the single best answer)

1. Which method does **not** require derivatives?
   (a) Newton (b) Nelder–Mead (c) BFGS (d) Gauss–Newton

2. Newton's method converges in how many iterations on a strictly convex quadratic?
   (a) 1 (b) 2 (c) $n$ (d) arbitrarily many

3. A Nelder–Mead simplex in $\mathbb R^n$ has how many vertices?
   (a) $n$ (b) $n+1$ (c) $2^n$ (d) $2n$

4. BFGS is classified as a:
   (a) direct-search method (b) second-order method (c) quasi-Newton method (d) least-squares method

5. The **secant condition** for BFGS (Hessian-approximation form) is:
   (a) $B_{k+1} y_k = s_k$ (b) $B_{k+1} s_k = y_k$ (c) $B_{k+1} = \nabla^2 f$ (d) $B_{k+1} = I$

6. Gauss–Newton replaces the Hessian $\nabla^2 f$ by:
   (a) $J^\top J$ (b) $J^{-1}$ (c) $J^\top J + \lambda I$ (d) the identity

7. Levenberg–Marquardt adds $\lambda I$ to which matrix?
   (a) $J$ (b) $J^\top J$ (c) $\nabla^2 f$ (d) $H_k$

8. As $\lambda \to 0$ in LM, the step approaches:
   (a) the Gauss–Newton step (b) a gradient-descent step (c) zero (d) a random direction

9. As $\lambda \to \infty$ in LM, the step approaches:
   (a) the Gauss–Newton step (b) a scaled gradient-descent step (c) a Newton step (d) zero

10. On the Rosenbrock function, which method typically needs **the most** function evaluations?
    (a) Newton (b) BFGS (c) gradient descent (d) Nelder–Mead

11. Nelder–Mead is **not guaranteed** to converge because:
    (a) it uses too few function evaluations
    (b) the simplex can collapse before reaching a stationary point
    (c) it requires gradients
    (d) it is a special case of Newton

12. Which is the **right** choice for a parameter-fitting problem with $m=1000$ residuals and $n=5$ parameters?
    (a) Nelder–Mead (b) BFGS on $\tfrac12\|r\|^2$ (c) Levenberg–Marquardt (d) Pattern search

### B. Fill in the blank

1. A Nelder–Mead simplex in $\mathbb R^n$ has \_\_\_\_ vertices.
2. Newton's step solves the linear system $H_k d_k = \_\_\_\_\_$.
3. In BFGS, $s_k = x_{k+1} - x_k$ and $y_k = \_\_\_\_\_$.
4. The BFGS secant condition is $B_{k+1} s_k = \_\_\_\_$.
5. The Gauss–Newton step satisfies $(J_k^\top J_k)\,\Delta\theta_k = \_\_\_\_\_$.
6. The Levenberg–Marquardt step satisfies $(J_k^\top J_k + \lambda_k I)\,\Delta\theta_k = \_\_\_\_\_$.
7. Hooke–Jeeves alternates an \_\_\_\_\_ move with a \_\_\_\_\_ move.
8. Newton's method has local \_\_\_\_\_ convergence near a strict minimum.
9. BFGS with exact line search converges in exactly \_\_\_\_ iterations on an $n$-dimensional quadratic.
10. The family of methods that only evaluates $f$ (never $\nabla f$) is called \_\_\_\_\_ search.

### C. Define / explain

1. Explain the **secant condition** for quasi-Newton methods and why it is the minimal “Hessian-mimicking” condition.
2. State the **Newton step** and its **local convergence rate**. Give one situation in which plain Newton can fail and how a line search fixes it.
3. Describe the **Nelder–Mead** reflection / expansion / contraction / shrink operations. What is the standard parameter tuple $(\alpha,\gamma,\rho,\sigma)$?
4. For a nonlinear least-squares problem, derive the Gauss–Newton step from the Taylor expansion of $r(\theta)$ and explain under what condition it is second-order convergent.
5. Compare **Gauss–Newton** and **Levenberg–Marquardt** in two or three sentences. Under what condition does LM reduce to GN, and under what condition does it reduce to gradient descent?
6. Compare **BFGS** and **Newton** in terms of cost per iteration and convergence rate. When would you prefer L-BFGS to either?

---

## 2.10 Answer key

??? success "A. Multiple choice"
    1. **(b)** — Nelder–Mead is a direct-search method.
    2. **(a)** — One step for a strictly convex quadratic.
    3. **(b)** — $n+1$ affinely independent vertices.
    4. **(c)** — Quasi-Newton; approximates the Hessian from gradients.
    5. **(b)** — $B_{k+1} s_k = y_k$.
    6. **(a)** — $J^\top J$ (dropping the residual-curvature term).
    7. **(b)** — $J^\top J$.
    8. **(a)** — $\lambda \to 0$: pure Gauss–Newton.
    9. **(b)** — $\lambda \to \infty$: gradient descent with step $\sim 1/\lambda$.
    10. **(c)** — Gradient descent is the slowest of the four on Rosenbrock (valley-following).
    11. **(b)** — The simplex can collapse at a non-stationary point (McKinnon counterexample).
    12. **(c)** — LM exploits the residual structure and the small $n$.

??? success "B. Fill in the blank"
    1. $n+1$
    2. $-\nabla f(x_k)$
    3. $\nabla f(x_{k+1}) - \nabla f(x_k)$
    4. $y_k$
    5. $-J_k^\top r_k$
    6. $-J_k^\top r_k$
    7. exploratory; pattern (or: pattern; exploratory — either order is acceptable in the standard “explore-then-pattern” description)
    8. quadratic
    9. $n$
    10. direct

??? success "C. Define / explain — model answers"
    1. The secant condition requires the Hessian approximation $B_{k+1}$ to reproduce the observed change in gradient over the last step: $B_{k+1}(x_{k+1}-x_k) = \nabla f(x_{k+1}) - \nabla f(x_k)$. It is the minimal condition a “Hessian-like” matrix must satisfy in order to be consistent with a first-order Taylor expansion of the gradient. Infinitely many symmetric matrices satisfy it — the different quasi-Newton methods (BFGS, DFP, SR1) are distinguished by *which* one they pick, e.g. by a minimum-change property.
    2. Newton step: $\nabla^2 f(x_k) d_k = -\nabla f(x_k)$, $x_{k+1} = x_k + d_k$. Local convergence is **quadratic** if $\nabla^2 f(x^\star) \succ 0$ and $x_0$ is close enough. Fails when $\nabla^2 f(x_k)$ is singular (no unique step) or indefinite (the step may not be a descent direction); a line search on the direction $d_k$ restores global convergence under mild assumptions while preserving local superlinear convergence.
    3. Nelder–Mead operations: **reflection** $r = \bar x + \alpha(\bar x - x_{\text{worst}})$ with $\alpha=1$; **expansion** $e = \bar x + \gamma(\bar x - x_{\text{worst}})$ with $\gamma=2$ (attempted when $f(r) < f(x_{\text{best}})$); **contraction** $c = \bar x + \rho(x_{\text{worst}} - \bar x)$ with $\rho=\tfrac12$ (when reflection is not better than the second-worst); **shrink** $x_i \leftarrow x_{\text{best}} + \sigma(x_i - x_{\text{best}})$ with $\sigma=\tfrac12$ (when contraction also fails). Standard tuple: $(\alpha,\gamma,\rho,\sigma) = (1, 2, \tfrac12, \tfrac12)$.
    4. Write $r(\theta + \Delta) \approx r(\theta) + J(\theta)\Delta$. Then $\tfrac12\|r(\theta + \Delta)\|^2$ is minimized by $\Delta$ solving $J^\top J \Delta = -J^\top r$ — this is the Gauss–Newton step. It is second-order convergent when the residual at the solution is zero (then the dropped term $\sum r_i \nabla^2 r_i$ vanishes at the solution), otherwise the convergence is only linear with rate governed by the size of the dropped term.
    5. Gauss–Newton solves $J^\top J \Delta = -J^\top r$; LM solves $(J^\top J + \lambda I)\Delta = -J^\top r$ with $\lambda$ adapted. LM reduces to GN when $\lambda \to 0$ (trust in the linearization) and to gradient descent when $\lambda \to \infty$ ($\Delta \approx -(1/\lambda) J^\top r = -(1/\lambda)\nabla f$). In practice LM starts with moderate $\lambda$, decreases it on successful steps, increases it on unsuccessful ones.
    6. Newton: $O(n^3)$ per iteration (dense linear solve) + Hessian cost, quadratic convergence. BFGS: $O(n^2)$ memory, $O(n^2)$ per iteration, superlinear convergence. L-BFGS stores only the last $m$ pairs $(s_k, y_k)$ — $O(mn)$ memory, $O(mn)$ per iteration — and is the method of choice when $n$ is large ($10^3$–$10^6$).

---

## Further reading

- Nocedal & Wright, *Numerical Optimization*, 2nd ed. — Chapters 2 (Newton), 3 (line search), 6 (quasi-Newton), 10 (nonlinear least squares).
- Boyd & Vandenberghe, *Convex Optimization* — §9.5 (Newton), §10 (equality-constrained Newton).
- Lagarias, Reeds, Wright & Wright, “Convergence Properties of the Nelder–Mead Simplex Method in Low Dimensions”, *SIAM J. Optim.* 9(1), 1998.
- McKinnon, “Convergence of the Nelder–Mead Simplex Method to a Nonstationary Point”, *SIAM J. Optim.* 9(1), 1998.
- Moré, “The Levenberg–Marquardt Algorithm: Implementation and Theory”, in *Numerical Analysis* (Lecture Notes in Mathematics 630), 1978.
- Nocedal, “Updating Quasi-Newton Matrices with Limited Storage”, *Math. Comp.* 35(151), 1980 (L-BFGS).
- Conn, Gould & Toint, *Trust-Region Methods*, SIAM, 2000 — for the trust-region interpretation of LM.
- SciPy `optimize` documentation, Apache Commons Math `optim` package, Math.NET `Numerics.Optimization`.

---

## Suggested tutorial problems for this chapter

1. **Nelder–Mead by hand.** Starting from the simplex $\{(0,0),(1,0),(0,1)\}$ on $f(x,y)=(x-2)^2+(y-3)^2$, trace two full iterations (sort → reflect → expand / contract / shrink). Verify your result against `scipy.optimize.minimize(method='Nelder-Mead')`.
2. **Newton with a line search.** Implement damped Newton on the Rosenbrock function from $(-1,1)$. Compare iteration counts with plain Newton and with plain gradient descent. Plot $\|x_k - (1,1)\|$ on a log axis vs $k$.
3. **BFGS update.** For $f(x_1,x_2)=x_1^2+2x_2^2$ starting from $(1,1)$, verify by hand that the BFGS update of Example 2.4 satisfies $H_1 s_0 = y_0$ and $H_1 \succ 0$.
4. **Quasi-Newton failure.** Show that the BFGS update breaks down when $y_k^\top s_k \le 0$. Construct a one-dimensional counterexample where this occurs with an arbitrary step (no line search).
5. **Gauss–Newton at work.** Fit $y = a e^{b x}$ to $x = (0,1,2,3)$, $y = (1, 3, 7, 15)$. Use both GN and LM (with $\lambda_0=1$, adaptive). Report the number of iterations each takes.
6. **Rank-deficient Jacobian.** Construct a nonlinear least-squares problem with $m > n$ where $J$ is rank-deficient at the starting point. Show that Gauss–Newton fails and LM still converges.
7. **Global optimization benchmark.** Compare Nelder–Mead, BFGS, and Newton on Rastrigin's function $f(x) = 10n + \sum[x_i^2 - 10\cos(2\pi x_i)]$ in $n=2$ from three different starts. Comment on which methods find the global minimum and why.
