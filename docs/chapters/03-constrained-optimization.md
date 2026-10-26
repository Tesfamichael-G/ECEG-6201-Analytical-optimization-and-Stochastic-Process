# Chapter 3 — Constrained Optimization

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 tutorial · **Home study:** 5 h
    
    **Coverage (syllabus):** penalty functions · barrier functions · augmented
    Lagrangian · quadratic programming — focus on QP with equality constraints,
    with an active-set preview for inequalities.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Formulate the constraint set $g_i(x) \le 0$, $h_j(x) = 0$ and identify feasible, strictly feasible, and boundary points.
2. Build **exterior penalty functions**, prove that the solution of the penalized problem converges to the constrained optimum as the penalty parameter grows, and quantify the error $O(1/\mu)$.
3. Build **logarithmic barrier functions**, explain why the barrier requires a strictly feasible starting point, and show that the error is $O(\mu)$ as $\mu \to 0^+$.
4. Derive and iterate the **augmented Lagrangian** (method of multipliers), and explain why it avoids the ill-conditioning of pure penalties.
5. Solve a **quadratic program with equality constraints** by forming and solving the KKT linear system; state existence and uniqueness conditions.
6. Describe the **active-set method** for QP with inequality constraints.
7. Implement each of the above in Python, Java, or C#.

---

## 3.1 Problem statement

Minimize $f : \mathbb{R}^n \to \mathbb{R}$ subject to

$$
\begin{aligned}
& \text{minimize}   && f(x) \\
& \text{subject to} && g_i(x) \le 0, \quad i = 1,\dots,p, \\
&                   && h_j(x) = 0, \quad j = 1,\dots,q.
\end{aligned}
$$

The **feasible set** is $\Omega = \{ x : g_i(x) \le 0,\; h_j(x) = 0 \}$.

| Concept | Definition |
|---|---|
| Feasible point | $x \in \Omega$ |
| Strictly feasible point | $g_i(x) < 0$ for all $i$, $h_j(x) = 0$ for all $j$ |
| Active inequality | $g_i(x) = 0$ at a feasible point |
| Inactive inequality | $g_i(x) < 0$ |
| Local minimum | $f(x^\star) \le f(x)$ for all feasible $x$ near $x^\star$ |

In Chapter 1 you saw the **KKT conditions**:

$$
\nabla f(x^\star) + \sum_i \mu_i^\star \nabla g_i(x^\star) + \sum_j \lambda_j^\star \nabla h_j(x^\star) = 0,
$$

$$
\mu_i^\star \ge 0,\qquad \mu_i^\star\, g_i(x^\star) = 0,\qquad g_i(x^\star) \le 0,\qquad h_j(x^\star) = 0.
$$

The methods in this chapter all converge to a KKT point — the difference is *how* they enforce feasibility.

---

## 3.2 Penalty functions — exterior methods

### Idea

Add a penalty for violating the constraints and minimize the penalized function **without constraints**:

$$
\boxed{\;P(x, \mu) = f(x) + \frac{\mu}{2} \sum_{i=1}^{p} \bigl[\max(0,\ g_i(x))\bigr]^2 + \frac{\mu}{2} \sum_{j=1}^{q} h_j(x)^2\;}
$$

with penalty parameter $\mu > 0$. Iterate $\mu \to \infty$: the penalised minimiser $x_\mu$ converges to the constrained optimum $x^\star$.

**Why “exterior”.** Because the iterates stay **outside** the feasible set (or on its boundary); the penalty pulls them back in as $\mu$ increases.

### Worked example 3.1 — exterior penalty on an inequality

Minimize $f(x) = (x-2)^2$ subject to $x \le 1$.

The unconstrained minimiser is $x = 2$, which is **infeasible** ($2 \le 1$ is false). So the constraint is active at the solution.

**Penalized problem.** With $g(x) = x-1 \le 0$:

$$
P(x, \mu) = (x-2)^2 + \frac{\mu}{2}\bigl[\max(0, x-1)\bigr]^2.
$$

For $x \le 1$: $P = (x-2)^2$, decreasing in $x$ ⇒ best is $x = 1$.
For $x > 1$:

$$
P = (x-2)^2 + \tfrac{\mu}{2}(x-1)^2.
$$

**First-order condition on $x > 1$:**

$$
\frac{dP}{dx} = 2(x-2) + \mu(x-1) = 0.
$$

Solving:

$$
2x - 4 + \mu x - \mu = 0 \;\Longrightarrow\; x_\mu = \frac{4+\mu}{2+\mu}.
$$

**Check convergence:** $x_\mu \to 1$ as $\mu \to \infty$. ✓ The constrained optimum is $x^\star = 1$ with $f^\star = 1$.

| $\mu$ | $x_\mu$ | $\|x_\mu - x^\star\|$ |
|:---:|:---:|:---:|
| 1 | 1.6667 | 6.67·10⁻¹ |
| 10 | 1.1667 | 1.67·10⁻¹ |
| 10² | 1.0196 | 1.96·10⁻² |
| 10³ | 1.0020 | 2.00·10⁻³ |
| 10⁴ | 1.0002 | 2.00·10⁻⁴ |

The error is $O(1/\mu)$ — to halve the error, double $\mu$ — which is why exterior penalties are *slow* near the solution and get numerically fragile at $\mu \sim 10^{10}$.

### Worked example 3.2 — exterior penalty on an equality

Minimize $f(x,y) = x^2 + y^2$ subject to $x + 2y = 4$.

**KKT (reference answer).**

$$
\mathcal L = x^2 + y^2 + \lambda(x + 2y - 4),\qquad
2x = -\lambda,\qquad 2y = -2\lambda,\qquad x + 2y = 4.
$$

The first two equations give $y = 2x$; substituting into the constraint, $x + 4x = 4$, so $x^\star = 0.8$, $y^\star = 1.6$, $\lambda^\star = -1.6$, $f^\star = 0.64 + 2.56 = 3.2$.

**Penalized problem.**

$$
P(x, y, \mu) = x^2 + y^2 + \frac{\mu}{2}(x + 2y - 4)^2.
$$

**First-order conditions.**

$$
\frac{\partial P}{\partial x} = 2x + \mu(x + 2y - 4) = 0,\qquad
\frac{\partial P}{\partial y} = 2y + 2\mu(x + 2y - 4) = 0.
$$

From the first equation, $\mu(x + 2y - 4) = -2x$. Substituting into the second:

$$
2y = -2\bigl(-2x\bigr) = 4x \;\Longrightarrow\; y = 2x.
$$

With $y = 2x$ in the first equation:

$$
2x + \mu(5x - 4) = 0 \;\Longrightarrow\; x_\mu = \frac{4\mu}{5\mu + 2}.
$$

**As $\mu \to \infty$:** $x_\mu \to 4/5 = 0.8$, $y_\mu \to 1.6$ — recovering the KKT solution. ✓

| $\mu$ | $x_\mu$ | $y_\mu$ | $\|(x_\mu, y_\mu) - (0.8, 1.6)\|$ |
|:---:|:---:|:---:|:---:|
| 1 | 0.571 | 1.143 | 5.11·10⁻¹ |
| 10 | 0.769 | 1.538 | 7.02·10⁻² |
| 10² | 0.797 | 1.594 | 7.34·10⁻³ |
| 10³ | 0.800 | 1.599 | 7.38·10⁻⁴ |

Again $O(1/\mu)$ convergence.

### Properties and pitfalls

!!! warning "Ill-conditioning"
    The Hessian of $P$ at the penalised minimiser has eigenvalues that scale like
    $\mu$ in directions normal to the active constraints — so as $\mu \to \infty$,
    the condition number grows like $\mu$. Simple unconstrained solvers (Newton,
    BFGS) will struggle to solve the inner problem accurately at large $\mu$.

!!! tip "When the constraint is inactive"
    If the unconstrained minimiser of $f$ already satisfies $g_i(x) \le 0$, then
    for any $\mu$ the penalised minimiser equals it. **Complementary slackness
    automatically** — no penalty is paid.

---

## 3.3 Barrier functions — interior methods

### Idea

Instead of penalising violations from outside, add a term that goes to $+\infty$ as the boundary is approached, keeping iterates **strictly inside** the feasible set:

$$
\boxed{\;B(x, \mu) = f(x) - \mu \sum_{i=1}^{p} \log\bigl(-g_i(x)\bigr)\;}
$$

with barrier parameter $\mu > 0$, and $\mu \to 0^+$ as the method converges. This is the **logarithmic barrier**, the basis of modern interior-point methods.

**Why “interior”.** Because the barrier is defined only when every $g_i(x) < 0$ strictly. The iterates stay inside the feasible set — the opposite of the penalty method.

### Worked example 3.3 — logarithmic barrier on a simple inequality

Minimize $f(x) = x$ subject to $x \ge 1$ (equivalently $g(x) = 1 - x \le 0$).

Unconstrained minimum would be $x \to -\infty$; the constraint forces $x^\star = 1$ with $f^\star = 1$.

**Barrier problem.**

$$
B(x, \mu) = x - \mu \log(x - 1), \qquad x > 1.
$$

**First-order condition.**

$$
B'(x) = 1 - \frac{\mu}{x-1} = 0 \;\Longrightarrow\; x_\mu = 1 + \mu.
$$

**Convergence.** $x_\mu \to 1$ as $\mu \to 0^+$; the gap is exactly $\mu$. ✓

**Optimal barrier value.**

$$
B(x_\mu, \mu) = (1 + \mu) - \mu \log \mu \;\xrightarrow{\mu \to 0^+}\; 1.
$$

| $\mu$ | $x_\mu$ | $\|x_\mu - x^\star\|$ |
|:---:|:---:|:---:|
| 1 | 2.00 | 1.00 |
| 0.5 | 1.50 | 0.50 |
| 0.1 | 1.10 | 0.10 |
| 0.01 | 1.01 | 0.01 |

The **duality gap** between the barrier value and the true optimum is $O(\mu)$ — so to reach tolerance $\varepsilon$, use $\mu \sim \varepsilon$. Contrast with the exterior penalty, which needs $\mu \sim 1/\varepsilon$.

### Starting point and Slater's condition

!!! note "Requirement"
    The barrier is well-defined only at points where $g_i(x) < 0$ **for all** $i$ — a **strictly feasible** point. If none exists, the barrier method cannot start.

!!! abstract "Slater's condition"
    There exists a point $x_0$ with $g_i(x_0) < 0$ for all $i$ and $h_j(x_0) = 0$ for all $j$.  
    Slater's condition guarantees strong duality (the KKT multipliers exist and the
    duality gap vanishes in the limit); it is the standard assumption for interior-point
    methods on convex problems.

### Convergence proof sketch

Differentiate $B$ at the barrier minimiser $x_\mu$:

$$
0 = \nabla B(x_\mu, \mu) = \nabla f(x_\mu) - \mu \sum_{i=1}^{p} \frac{1}{-g_i(x_\mu)} \nabla g_i(x_\mu).
$$

Compare with the KKT stationarity condition, $\nabla f + \sum \mu_i \nabla g_i = 0$. This suggests identifying

$$
\mu_i(\mu) := \frac{\mu}{-g_i(x_\mu)} > 0.
$$

The **complementary-slackness gap** at the barrier minimiser is then

$$
\sum_i \mu_i(\mu)\, g_i(x_\mu) = \sum_i \frac{\mu\, g_i(x_\mu)}{-g_i(x_\mu)} = -\mu\, p,
$$

i.e. the total duality gap is $p\,\mu$ — hence $O(\mu)$. That is the mathematical content behind “take $\mu \to 0$ to reach the constrained optimum.”

---

## 3.4 Augmented Lagrangian — method of multipliers

### Motivation

Exterior penalty methods are ill-conditioned as $\mu \to \infty$; barrier methods need a strictly feasible start. The **augmented Lagrangian** (also called the **method of multipliers**) fixes both problems by adding a *Lagrange multiplier term* alongside the quadratic penalty:

$$
\boxed{\;\mathcal L_\rho(x, \lambda) = f(x) + \sum_{j=1}^{q} \lambda_j h_j(x) + \frac{\rho}{2} \sum_{j=1}^{q} h_j(x)^2\;}
$$

with penalty parameter $\rho > 0$ and multiplier estimate $\lambda$.

### Algorithm (equality-constrained case)

```
Choose x₀, λ₀, ρ > 0, tolerance ε
repeat:
    x_{k+1} ← argmin_x  L_ρ(x, λ_k)      (unconstrained, or box)
    λ_{k+1} ← λ_k + ρ · h(x_{k+1})       (multiplier update)
until ‖h(x_{k+1})‖ < ε
```

### Why it works — and why it is better than pure penalty

Consider the augmented Lagrangian and differentiate with respect to $x$:

$$
\nabla_x \mathcal L_\rho = \nabla f(x) + \sum_j \bigl[\lambda_j + \rho\, h_j(x)\bigr] \nabla h_j(x).
$$

Compare with the KKT condition $\nabla f(x^\star) + \sum \lambda_j^\star \nabla h_j(x^\star) = 0$. If we set the effective multiplier to be

$$
\lambda_j^{\text{eff}} := \lambda_j + \rho\, h_j(x),
$$

then stationarity of the augmented Lagrangian **is exactly the KKT stationarity condition** — provided $\lambda_j^{\text{eff}} = \lambda_j^\star$. The multiplier update $\lambda \leftarrow \lambda + \rho\, h(x)$ is a fixed-point iteration that drives $\lambda \to \lambda^\star$ and, simultaneously, $h(x) \to 0$.

!!! tip "The key advantage"
    Once $\lambda = \lambda^\star$ is reached, the constrained minimiser is
    the **exact unconstrained minimiser** of $\mathcal L_\rho$ for *any* fixed
    $\rho > 0$ — **no need to drive $\rho \to \infty$**. This is why the
    augmented Lagrangian avoids the ill-conditioning of pure penalties.

### Worked example 3.4 — augmented Lagrangian on a one-dimensional problem

Minimize $f(x) = x^2$ subject to $x = 3$.

**Constrained optimum.** $x^\star = 3$, $\lambda^\star = -2x^\star = -6$, $f^\star = 9$.

**Augmented Lagrangian.**

$$
\mathcal L_\rho(x, \lambda) = x^2 + \lambda(x - 3) + \frac{\rho}{2}(x - 3)^2.
$$

**Inner minimisation.** For fixed $\lambda$:

$$
\frac{\partial \mathcal L_\rho}{\partial x} = 2x + \lambda + \rho(x - 3) = 0
\;\Longrightarrow\; x = \frac{3\rho - \lambda}{2 + \rho}.
$$

**Multiplier update.** $\lambda \leftarrow \lambda + \rho(x - 3)$.

**Iterations** ($\rho = 1$, $\lambda_0 = 0$):

| $k$ | $x_k$ | $\lambda_k$ | $\|x_k - 3\|$ |
|:---:|:---:|:---:|:---:|
| 0 | 1.0000 | 0.0000 | 2.000 |
| 1 | 1.6667 | −2.0000 | 1.333 |
| 2 | 2.1111 | −3.3333 | 0.889 |
| 3 | 2.4074 | −4.2222 | 0.593 |
| 4 | 2.6049 | −4.8148 | 0.395 |
| 5 | 2.7366 | −5.2100 | 0.263 |

**Convergence analysis.** With the closed form $x_k = (3\rho - \lambda_k)/(2 + \rho)$, one verifies by algebra:

$$
\lambda_{k+1} - \lambda^\star = \frac{2}{2+\rho}\,(\lambda_k - \lambda^\star).
$$

So the multiplier converges **geometrically with rate** $2/(2+\rho)$ — for $\rho = 1$, rate $= 2/3$. With a larger $\rho$, faster convergence but a more ill-conditioned inner problem; with a smaller $\rho$, slower convergence but a well-conditioned inner problem. The trick is to choose $\rho$ (or adapt it) to balance these.

### Combined equality + inequality

For inequalities $g_i(x) \le 0$, add a barrier-like penalty on the *positive part* of the multiplier-adjusted slack:

$$
\mathcal L_\rho(x, \lambda, \mu) = f(x) + \sum_j \lambda_j h_j(x) + \frac{\rho}{2}\sum_j h_j(x)^2 + \frac{1}{2\rho} \sum_i \Bigl[\max\bigl(0,\ \mu_i + \rho\, g_i(x)\bigr)^2 - \mu_i^2\Bigr],
$$

with updates

$$
\lambda \leftarrow \lambda + \rho\, h(x), \qquad \mu_i \leftarrow \max\bigl(0,\ \mu_i + \rho\, g_i(x)\bigr).
$$

This is the **ALGENCAN / LANCELOT** formulation and is what practical solvers such as `Ipopt`, `Algencan`, and `MINOS` build on.

---

## 3.5 Quadratic programming (QP)

### Standard form

$$
\begin{aligned}
& \text{minimize}   && \tfrac{1}{2} x^\top Q x + c^\top x \\
& \text{subject to} && A x = b, \\
&                   && G x \le h.
\end{aligned}
$$

with $Q = Q^\top \in \mathbb{R}^{n \times n}$, $c \in \mathbb{R}^n$, $A \in \mathbb{R}^{m \times n}$, $G \in \mathbb{R}^{r \times n}$.

**Key fact.** If $Q \succ 0$ and the feasible set is non-empty, QP is a **convex** problem — any KKT point is a global minimum, and solvers converge reliably.

### 3.5.1 QP with equality constraints only

Minimize $\tfrac12 x^\top Q x + c^\top x$ subject to $A x = b$.

**Lagrangian.**

$$
\mathcal L(x, \lambda) = \tfrac12 x^\top Q x + c^\top x + \lambda^\top(Ax - b).
$$

**KKT conditions (linear in $(x, \lambda)$).**

$$
Q x + c + A^\top \lambda = 0, \qquad A x = b.
$$

Stack into a single linear system — the **KKT matrix**:

$$
\boxed{\;
\begin{bmatrix} Q & A^\top \\ A & 0 \end{bmatrix}
\begin{bmatrix} x \\ \lambda \end{bmatrix}
= \begin{bmatrix} -c \\ b \end{bmatrix}
\;}
$$

**Existence and uniqueness.**

| Case | Solution set |
|---|---|
| $Q \succ 0$, $A$ full row rank | **Unique** $(x^\star, \lambda^\star)$ |
| $Q \succ 0$, $A$ rank-deficient | Unique $x^\star$; $\lambda$ not unique |
| $Q$ indefinite, $A$ full row rank | Unique stationary point; **not necessarily** a minimum |
| $Q$ indefinite, problem unbounded | No solution |

!!! tip "Practical solvers"
    For small problems, just solve the KKT matrix directly (`numpy.linalg.solve`).  
    For larger problems, exploit structure: solve the Schur complement
    $A Q^{-1} A^\top \lambda = A Q^{-1} c - b$ first, then recover
    $x = Q^{-1}(-c - A^\top \lambda)$. This reduces the dense system from
    $(n+m)\times(n+m)$ to $m\times m$.

### Worked example 3.5 — QP with an equality constraint, solved by KKT

Minimize

$$
f(x, y) = x^2 + y^2 + xy - 2x + y \quad \text{subject to} \quad x + y = 3.
$$

**Standard form.** Write as $\tfrac12 z^\top Q z + c^\top z$ with $z = (x, y)$:

$$
Q = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix},\qquad
c = \begin{bmatrix} -2 \\ 1 \end{bmatrix},\qquad
A = \begin{bmatrix} 1 & 1 \end{bmatrix},\qquad
b = 3.
$$

Since $Q \succ 0$ (eigenvalues $1, 3$) and $A$ has rank $1$, the solution is unique.

**KKT system.**

$$
\begin{bmatrix}
2 & 1 & 1 \\
1 & 2 & 1 \\
1 & 1 & 0
\end{bmatrix}
\begin{bmatrix} x \\ y \\ \lambda \end{bmatrix}
= \begin{bmatrix} 2 \\ -1 \\ 3 \end{bmatrix}.
$$

**Row reduction.** Row 1 − Row 2 gives $x - y = 3$. With Row 3 ($x + y = 3$), we get $x = 3$, $y = 0$. Then Row 1 gives $\lambda = 2 - 2 \cdot 3 - 0 = -4$.

**Verification.** $x^\star = (3, 0)$, $\lambda^\star = -4$, $f^\star = 9 + 0 + 0 - 6 + 0 = 3$.

**Direct check** by eliminating the constraint: with $y = 3 - x$,

$$
f(x) = x^2 + (3-x)^2 + x(3-x) - 2x + (3-x) = x^2 - 6x + 12.
$$

$df/dx = 2x - 6 = 0 \Rightarrow x = 3$, $y = 0$, $f = 3$. ✓

### 3.5.2 QP with inequality constraints — active-set preview

The QP with both equalities and inequalities is handled by identifying the **active set** $\mathcal W$ (indices of constraints satisfied as equalities at the solution) and solving the equality-constrained KKT system on $\mathcal W$.

**Active-set algorithm (conceptual).**

```
Choose a feasible starting point x₀ and an initial working set W
repeat:
    solve the equality-constrained QP with Ax = b, G_W x = h_W
        -> candidate step d
    if d = 0 (current point is optimal for W):
        compute multipliers μ for constraints in W
        if all μ_i ≥ 0:  STOP (KKT satisfied)
        else:  remove the most negative μ_i from W
    else:
        take a step in direction d, stopping at the first
        constraint that becomes active; add that constraint to W
until convergence
```

This is the classical **Wolfe / Goldfarb–Idnani** algorithm used in production QP solvers (`quadprog`, `OSQP` uses a different interior-point approach).

### Worked example 3.6 — QP with an inequality, solved by inspection

Minimize

$$
f(x,y) = (x-2)^2 + (y-2)^2 \quad \text{subject to} \quad x + y \le 3,\ x \ge 0,\ y \ge 0.
$$

**Step 1 — guess the active set.** The unconstrained minimiser is $(2, 2)$, but $2 + 2 = 4 > 3$ ⇒ the inequality $x + y \le 3$ is violated ⇒ it is active. Both non-negativity constraints hold at the solution (to be verified).

So $\mathcal W = \{x + y \le 3\}$.

**Step 2 — solve the equality-constrained QP** on $x + y = 3$. Substitute $y = 3 - x$:

$$
f(x) = (x - 2)^2 + (1 - x)^2,\qquad
\frac{df}{dx} = 2(x-2) - 2(1-x) = 4x - 6 = 0 \;\Longrightarrow\; x = \tfrac{3}{2}.
$$

So $x = 1.5$, $y = 1.5$, $f = 0.25 + 0.25 = 0.5$.

**Step 3 — verify KKT multipliers.**

$$
\nabla f(1.5, 1.5) = (2(1.5-2),\ 2(1.5-2)) = (-1, -1).
$$

For the constraint $x + y - 3 \le 0$, $\nabla g = (1, 1)$. Stationarity:

$$
\nabla f + \mu \nabla g = 0 \;\Longrightarrow\; (-1,-1) + \mu(1,1) = 0 \;\Longrightarrow\; \mu = 1 \ge 0. \;\checkmark
$$

For $x \ge 0$: at $x = 1.5 > 0$ the constraint is inactive, so its multiplier is $0$ (complementary slackness). Similarly for $y \ge 0$. All KKT conditions hold.

**Answer.** $x^\star = (1.5, 1.5)$, $f^\star = 0.5$, $\mu^\star_{x+y} = 1$, all other multipliers zero.

---

## 3.6 Computational implementation

The examples below solve the same three problems — (i) Example 3.4's augmented Lagrangian, (ii) Example 3.5's QP with equality constraint, and (iii) Example 3.6's inequality-constrained QP — using SciPy, Java, and C#. The point is to show the same algorithm expressed three ways.

### Python with SciPy

```python
import numpy as np
from scipy.optimize import minimize, minimize_scalar

# ---- (i) augmented Lagrangian for min x^2 s.t. x = 3 --------------------
def aug_lag_demo(rho=1.0, maxiter=30, tol=1e-8):
    x, lam = 1.0, 0.0
    for k in range(maxiter):
        # inner: unconstrained in x, fixed lam and rho
        f_inner = lambda x: x**2 + lam*(x - 3) + 0.5*rho*(x - 3)**2
        x = minimize_scalar(f_inner, bounds=(-50, 50), method='bounded').x
        gap = abs(x - 3)
        print(f"k={k:2d}  x={x:.6f}  λ={lam:+.6f}  |x-3|={gap:.2e}")
        if gap < tol:
            break
        lam += rho * (x - 3)

print("=== augmented Lagrangian ===")
aug_lag_demo()

# ---- (ii) QP with equality, KKT linear system --------------------------
Q = np.array([[2., 1.], [1., 2.]])
c = np.array([-2., 1.])
A = np.array([[1., 1.]])
b = np.array([3.])

n, m = Q.shape[0], A.shape[0]
KKT = np.block([[Q, A.T], [A, np.zeros((m, m))]])
rhs = np.concatenate([-c, b])
sol = np.linalg.solve(KKT, rhs)
xstar, lamstar = sol[:n], sol[n:]
print("\n=== QP with equality (KKT) ===")
print("x* =", xstar, "  λ* =", lamstar,
      "  f* =", 0.5*xstar@Q@xstar + c@xstar)

# ---- (iii) QP with inequality — scipy SLSQP -----------------------------
def f_qp(z):
    return (z[0]-2)**2 + (z[1]-2)**2
def df_qp(z):
    return np.array([2*(z[0]-2), 2*(z[1]-2)])

cons = [
    {'type': 'ineq', 'fun': lambda z: 3 - z[0] - z[1],
                     'jac': lambda z: np.array([-1., -1.])},
    {'type': 'ineq', 'fun': lambda z: z[0], 'jac': lambda z: np.array([1., 0.])},
    {'type': 'ineq', 'fun': lambda z: z[1], 'jac': lambda z: np.array([0., 1.])},
]
r = minimize(f_qp, x0=np.array([0.1, 0.1]), jac=df_qp,
             constraints=cons, method='SLSQP', options={'ftol': 1e-12})
print("\n=== QP with inequality (SLSQP) ===")
print("x* =", r.x, "  f* =", r.fun)

# ---- (iv) exterior penalty, educational --------------------------------
print("\n=== exterior penalty on min (x-2)^2 s.t. x <= 1 ===")
for mu in [1., 10., 100., 1000.]:
    P = lambda x: (x-2)**2 + 0.5*mu*max(0, x-1)**2
    x = minimize_scalar(P, bounds=(-5, 5), method='bounded').x
    print(f"μ={mu:6.0f}  x*={x:.6f}  exact={(4+mu)/(2+mu):.6f}")
```

**Expected output (abridged).**

```
=== augmented Lagrangian ===
k= 0  x=1.000000  λ=+0.000000  |x-3|=2.00e+00
k= 1  x=1.666667  λ=-2.000000  |x-3|=1.33e+00
...
k=12  x=2.999323  λ=-5.998645  |x-3|=6.77e-04
=== QP with equality (KKT) ===
x* = [3. 0.]   λ* = [-4.]   f* = 3.0
=== QP with inequality (SLSQP) ===
x* = [1.5 1.5]   f* = 0.5
```

### Java — augmented Lagrangian

```java
public class AugmentedLagrangian {

    // min  f(x) = x^2  s.t.  h(x) = x - 3 = 0
    static double innerMin(double lambda, double rho) {
        // d/dx [ x² + λ(x-3) + (ρ/2)(x-3)² ] = 0
        //  2x + λ + ρ(x-3) = 0  =>  x = (3ρ - λ)/(2+ρ)
        return (3*rho - lambda) / (2 + rho);
    }

    public static void main(String[] args) {
        double rho = 1.0, lam = 0.0, tol = 1e-8;
        for (int k = 0; k < 50; k++) {
            double x = innerMin(lam, rho);
            double gap = Math.abs(x - 3);
            System.out.printf("k=%2d  x=%.8f  λ=%+.6f  |x-3|=%.2e%n",
                              k, x, lam, gap);
            if (gap < tol) break;
            lam += rho * (x - 3);
        }
    }
}
```

### C# — exterior penalty + QP by KKT

```csharp
using System;

class ConstrainedOptimization
{
    // ---- (i) exterior penalty for min (x-2)^2 s.t. x <= 1 --------
    static double PenaltyMin(double mu) {
        // P(x) = (x-2)² + (μ/2) max(0, x-1)²
        // For x > 1: x = (4+μ)/(2+μ); for x ≤ 1: x = 1.
        double x = (4 + mu) / (2 + mu);
        return x > 1 ? x : 1.0;
    }

    // ---- (ii) QP with equality: solve KKT by Gaussian elimination ----
    static (double[] x, double[] lam) SolveQP(double[,] Q, double[] c,
                                              double[,] A, double[] b)
    {
        int n = Q.GetLength(0), m = A.GetLength(0);
        int N = n + m;
        var M = new double[N, N + 1];

        // top-left Q, top-right Aᵀ
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) M[i, j] = Q[i, j];
            for (int j = 0; j < m; j++) M[i, n + j] = A[j, i];
            M[i, N] = -c[i];
        }
        // bottom-left A, bottom-right 0
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) M[n + i, j] = A[i, j];
            for (int j = 0; j < m; j++) M[n + i, n + j] = 0.0;
            M[n + i, N] = b[i];
        }

        // Gauss-Jordan elimination
        for (int col = 0; col < N; col++) {
            int piv = col;
            for (int r = col + 1; r < N; r++)
                if (Math.Abs(M[r, col]) > Math.Abs(M[piv, col])) piv = r;
            for (int j = col; j <= N; j++)
                (M[col, j], M[piv, j]) = (M[piv, j], M[col, j]);

            double d = M[col, col];
            for (int j = col; j <= N; j++) M[col, j] /= d;
            for (int r = 0; r < N; r++) {
                if (r == col) continue;
                double factor = M[r, col];
                for (int j = col; j <= N; j++) M[r, j] -= factor * M[col, j];
            }
        }

        var sol = new double[N];
        for (int i = 0; i < N; i++) sol[i] = M[i, N];
        var x = new double[n]; var lam = new double[m];
        Array.Copy(sol, 0, x, 0, n);
        Array.Copy(sol, n, lam, 0, m);
        return (x, lam);
    }

    static void Main() {
        Console.WriteLine("=== exterior penalty ===");
        foreach (double mu in new[]{1.0, 10.0, 100.0, 1000.0})
            Console.WriteLine($"μ={mu,6:F0}  x*={PenaltyMin(mu):F6}  " +
                              $"exact={(4+mu)/(2+mu):F6}");

        Console.WriteLine("\n=== QP with equality ===");
        double[,] Q = { {2, 1}, {1, 2} };
        double[]  c = { -2, 1 };
        double[,] A = { {1, 1} };
        double[]  b = { 3 };
        var (x, lam) = SolveQP(Q, c, A, b);
        Console.WriteLine($"x* = ({x[0]:F4}, {x[1]:F4}),  λ* = {lam[0]:F4}");
    }
}
```

!!! tip "Production solvers you can use directly"
    - **Python:** `scipy.optimize.minimize` (`'SLSQP'`, `'trust-constr'`),
      `scipy.optimize.quadprog`-compatible `quadprog`, `cvxpy` for modeling.
    - **Java:** Apache Commons Math `QuadraticMultivariateFunction` + `QuadraticOptimizer`;
      Hipparchus; `ojalgo` (`ConvexSolver`).
    - **C#:** Math.NET Numerics has `BFGS` but no native QP; use `Cplex`/`Gurobi`
      via their .NET APIs, or `Google OR-Tools` (`Google.OrTools.LinearSolver`)
      for LPs and QPs.
    - **MATLAB / Julia:** `quadprog`, `fmincon`, `JuMP + Gurobi / OSQP`.

---

## 3.7 Summary — choosing between methods

| Method | Needs feasible start? | Sub-problem per iteration | Convergence | Best for |
|---|:---:|---|---|---|
| Exterior penalty | No (infeasible ok) | Unconstrained | $O(1/\mu)$ | Simple problems, teaching |
| Barrier (interior point) | **Yes** (strictly feasible) | Unconstrained + barrier | $O(\mu)$ | Large convex problems |
| Augmented Lagrangian | No | Unconstrained | Superlinear in $\lambda$ | Equality & mixed, production |
| QP with equality | Any | Linear solve (KKT) | Exact (one step) | QP with only equalities |
| QP active-set | Any | Equality-constrained QP | Finite | QP with inequalities |

**Rule of thumb.**

- If the problem is a **QP with only equality constraints**, solve the KKT linear system once.
- If it is a **general smooth nonlinear problem**, use an **augmented Lagrangian** (equalities) or a **logarithmic barrier** (inequalities, convex case) — avoid plain quadratic penalties in production.
- If it is a **QP with inequalities**, use an **active-set** or **interior-point** QP solver, not a penalty loop.

---

## 3.8 Quiz

### A. Multiple choice (choose the single best answer)

1. Which method requires a **strictly feasible** starting point?
   (a) exterior penalty (b) logarithmic barrier (c) augmented Lagrangian (d) QP active-set

2. In the exterior penalty method, the parameter $\mu$ should tend to:
   (a) $0^+$ (b) $1$ (c) $+\infty$ (d) any value

3. In the logarithmic barrier method, the parameter $\mu$ should tend to:
   (a) $0^+$ (b) $1$ (c) $+\infty$ (d) any value

4. The main purpose of the **augmented Lagrangian** is to avoid:
   (a) KKT conditions (b) ill-conditioning of the pure penalty method (c) complementarity (d) QP structure

5. The KKT system for a QP with $n$ variables and $m$ equality constraints has dimension:
   (a) $n \times n$ (b) $m \times m$ (c) $(n+m) \times (n+m)$ (d) $2n \times 2n$

6. For a QP $\min \tfrac12 x^\top Q x + c^\top x$ s.t. $Ax = b$, a unique solution exists iff:
   (a) $A$ is square (b) $Q \succ 0$ and $A$ has full row rank (c) $Q$ is positive semi-definite (d) $Q = I$

7. In the exterior penalty method with finite $\mu$, the constraint $g_i(x) \le 0$ is:
   (a) exactly satisfied (b) approximately satisfied, with error $O(1/\mu)$ (c) always violated (d) ignored

8. In the augmented Lagrangian $\mathcal L_\rho(x,\lambda) = f(x) + \lambda h(x) + \tfrac{\rho}{2}h(x)^2$, the update $\lambda \leftarrow \lambda + \rho\,h(x)$ drives the pair $(x, \lambda)$ toward:
   (a) the unconstrained minimum of $f$ (b) a KKT point of the constrained problem (c) the starting point (d) infinity

9. The barrier term $-\mu \log(-g_i(x))$ is defined only when:
   (a) $g_i(x) \ge 0$ (b) $g_i(x) < 0$ strictly (c) $g_i(x) = 0$ (d) no condition

10. Which is **not** a property of the exterior penalty method?
    (a) iterates may be infeasible (b) ill-conditioned as $\mu \to \infty$ (c) requires strictly feasible start (d) converges from outside the feasible set

11. In the active-set method for QP, an active constraint at the solution has a KKT multiplier satisfying:
    (a) $\mu_i = 0$ (b) $\mu_i \ge 0$ (c) $\mu_i \le 0$ (d) any sign

12. **Slater's condition** for the log-barrier method on a convex problem asserts that:
    (a) all constraints are equalities (b) there exists $x_0$ with $g_i(x_0) < 0$ for all $i$ (c) the objective is linear (d) $Q \succ 0$

### B. Fill in the blank

1. Exterior penalty: $P(x,\mu) = f(x) + \tfrac{\mu}{2}\sum_i [\max(0,\ \_\_\_\_)]^2$.
2. Barrier: $B(x,\mu) = f(x) - \mu \sum_i \log(\_\_\_\_)$.
3. Augmented Lagrangian (equality): $\mathcal L_\rho(x,\lambda) = f(x) + \lambda^\top h(x) + \_\_\_\_\,\|h(x)\|^2$.
4. Augmented-Lagrangian multiplier update: $\lambda_{k+1} = \lambda_k + \rho\,\_\_\_\_$.
5. As $\mu \to$ \_\_\_, the exterior penalty minimiser converges to the constrained optimum.
6. As $\mu \to$ \_\_\_, the barrier minimiser converges to the constrained optimum.
7. The QP KKT block matrix is $\begin{bmatrix} Q & A^\top \\ A & \_\_\_\_\end{bmatrix}$.
8. The barrier method needs a \_\_\_\_ starting point.
9. The exterior penalty method allows \_\_\_\_ iterates.
10. The augmented Lagrangian avoids the \_\_\_\_ of pure penalty methods.

### C. Define / explain

1. Compare **exterior penalty** and **log-barrier** methods in three respects: starting-point requirement, side of the boundary the iterates live on, and asymptotic error as the parameter is driven to its limit.
2. Explain the **ill-conditioning** of the exterior penalty method and how the augmented Lagrangian overcomes it. Include the effective-multiplier identity $\lambda^{\text{eff}} = \lambda + \rho\,h(x)$.
3. **Derive** the KKT linear system for a QP with equality constraints. State the existence/uniqueness condition and sketch the Schur-complement reduction.
4. Describe the **augmented-Lagrangian iteration** on the equality-constrained problem. Explain, in one sentence, why the multiplier update $\lambda \leftarrow \lambda + \rho\,h(x)$ is a fixed-point iteration toward the correct KKT multiplier.
5. State **Slater's condition** and explain why it is the standard assumption for interior-point methods. What can go wrong if it fails?
6. Describe the **active-set method** for QP with inequalities. In what sense does it generalise the pure equality-constrained QP of §3.5.1?

---

## 3.9 Answer key

??? success "A. Multiple choice"
    1. **(b)** — the log barrier is defined only where $g_i(x) < 0$ strictly.
    2. **(c)** — exterior penalty: $\mu \to \infty$.
    3. **(a)** — barrier: $\mu \to 0^+$.
    4. **(b)** — AL removes the need for $\mu \to \infty$.
    5. **(c)** — $n + m$.
    6. **(b)** — $Q \succ 0$ and $A$ full row rank.
    7. **(b)** — penalty error is $O(1/\mu)$, not zero for finite $\mu$.
    8. **(b)** — fixed point of the iteration is the KKT pair.
    9. **(b)** — strictly negative $g_i(x)$ inside the barrier.
    10. **(c)** — that is a property of the *barrier*, not the penalty.
    11. **(b)** — active inequalities have $\mu_i \ge 0$ at a minimum.
    12. **(b)** — strict feasibility is exactly Slater's condition for inequality-only problems.

??? success "B. Fill in the blank"
    1. $g_i(x)$
    2. $-g_i(x)$
    3. $\tfrac{\rho}{2}$
    4. $h(x_k)$
    5. $+\infty$
    6. $0^+$
    7. $0$
    8. strictly feasible (interior)
    9. infeasible
    10. ill-conditioning

??? success "C. Define / explain — model answers"
    1. *Starting point.* Barrier needs a strictly feasible $x_0$; penalty does not. *Side of the boundary.* Penalty iterates stay outside or on the boundary; barrier iterates stay strictly inside. *Asymptotic error.* Penalty error is $O(1/\mu)$ as $\mu \to \infty$; barrier error is $O(\mu)$ as $\mu \to 0^+$. Consequently the barrier reaches a given tolerance with a much smaller parameter excursion — one reason interior-point methods displaced penalty methods in the 1990s.
    2. The Hessian of the quadratic penalty has eigenvalues that scale like $\mu$ in constraint-normal directions, so its condition number grows like $\mu$ as $\mu \to \infty$; simple solvers then fail to solve the inner unconstrained problem accurately. The augmented Lagrangian adds a Lagrange-multiplier term $\lambda^\top h(x)$ so that at stationarity $\nabla f + (\lambda + \rho\,h(x))^\top \nabla h = 0$. If $\lambda$ is updated to equal the correct KKT multiplier $\lambda^\star$, the minimiser of $\mathcal L_\rho$ is exactly the constrained optimum for any fixed $\rho$ — hence no need to drive $\rho \to \infty$, hence no ill-conditioning. The identity $\lambda^{\text{eff}} = \lambda + \rho\,h(x)$ is what makes stationarity of $\mathcal L_\rho$ coincide with KKT stationarity once $\lambda = \lambda^\star$.
    3. Lagrangian $\mathcal L = \tfrac12 x^\top Qx + c^\top x + \lambda^\top(Ax-b)$. Gradient in $x$: $Qx + c + A^\top\lambda = 0$. Gradient in $\lambda$: $Ax - b = 0$. Stack into $\begin{bmatrix}Q & A^\top \\ A & 0\end{bmatrix}\begin{bmatrix}x\\\lambda\end{bmatrix} = \begin{bmatrix}-c\\b\end{bmatrix}$. Unique solution iff $Q \succ 0$ and $A$ has full row rank. Schur complement: eliminate $x$ via $x = Q^{-1}(-c - A^\top\lambda)$, substitute into $Ax = b$: $A Q^{-1}(A^\top\lambda + c) = -b$, so $\lambda$ is obtained from an $m\times m$ linear system, then $x$ is recovered.
    4. Choose $x_0, \lambda_0, \rho$. Inner step: $x_{k+1} = \arg\min_x \mathcal L_\rho(x, \lambda_k)$. Multiplier update: $\lambda_{k+1} = \lambda_k + \rho\,h(x_{k+1})$. At the fixed point, $h(x^\star) = 0$ (feasibility) and stationarity of $\mathcal L_\rho$ gives the KKT stationarity condition with multiplier $\lambda^\star$ — so the iteration converges exactly to the constrained KKT pair. The multiplier update is a gradient ascent on the dual function (equivalently, a fixed-point iteration driven by the primal residual $h(x_k)$).
    5. Slater's condition asserts existence of a strictly feasible point $x_0$ with $g_i(x_0) < 0$ for all $i$ (and $h_j(x_0) = 0$ for all equalities). It guarantees that the KKT multipliers exist and the duality gap vanishes in the limit of the interior-point method, so the barrier path converges to the true constrained optimum. If Slater's condition fails (no strictly feasible point), the barrier is undefined — no $x_0$ exists from which to start — and the duality gap may not converge to zero.
    6. The active-set method solves a sequence of *equality-constrained QPs*, one per working set $\mathcal W$ of constraints. At each step it either (i) drops a constraint whose multiplier is negative (leaving the set), or (ii) adds a constraint that is violated on the current step (entering the set). When a step finds $d = 0$ with all $\mu_i \ge 0$ for $i \in \mathcal W$, KKT is satisfied and the algorithm terminates. It generalises §3.5.1 by turning inequality constraints into equalities *only when they matter*, solving the KKT system of §3.5.1 on the reduced set — hence “active-set”.

---

## Further reading

- Bertsekas, *Constrained Optimization and Lagrange Multiplier Methods*, Academic Press, 1982 — the definitive reference on augmented Lagrangians.
- Bertsekas, *Nonlinear Programming*, 2nd ed., Athena Scientific — Chapters 4 (penalty/barrier), 5 (duality, AL).
- Nocedal & Wright, *Numerical Optimization*, 2nd ed. — §16 (penalty/AL), §17 (interior-point), §16.5 (QP active set).
- Boyd & Vandenberghe, *Convex Optimization* — §4.2 (linear/quadratic/QP), §5.3 (Lagrange duality), §11 (interior-point).
- Conn, Gould & Toint, *Trust-Region Methods*, SIAM 2000 — §14 (AL).
- Fletcher, *Practical Methods of Optimization*, 2nd ed. — §12 (QP, active set), §12.3 (Goldfarb–Idnani).
- Wright, *Primal-Dual Interior-Point Methods*, SIAM, 1997.
- Birge & Louveaux, *Introduction to Stochastic Programming* — for QP applied to stochastic control.
- SciPy `optimize` (`SLSQP`, `trust-constr`) and `cvxpy` documentation; Apache Commons Math and `ojalgo` documentation; Gurobi / OSQP documentation.

---

## Suggested tutorial problems for this chapter

1. **Penalty iteration by hand.** Minimize $f(x) = (x-3)^2$ subject to $x \ge 2$, using the exterior penalty with $\mu \in \{1, 10, 10^2, 10^3\}$. Tabulate $x_\mu$ and confirm $x_\mu \to 2$ as $\mu \to \infty$. Compare with the barrier $B(x, \mu) = (x-3)^2 - \mu\log(x-2)$ and note the sign and rate of convergence.
2. **Barrier starting point.** Solve $\min\, x_1 + x_2$ subject to $x_1^2 + x_2^2 \le 1$ using the log barrier. Find a strictly feasible starting point. What happens at $\mu = 1$, $\mu = 0.1$, $\mu = 0.01$?
3. **Augmented Lagrangian with both signs.** Minimize $f(x) = x^2 - 2x$ subject to $x = 1$. Take $\rho = 1$, $\lambda_0 = 0$ and iterate; compute the analytic convergence rate from the formula $(\lambda_{k+1} - \lambda^\star) = \frac{2}{2+\rho}(\lambda_k - \lambda^\star)$.
4. **Mixed constraints.** Minimize $x^2 + y^2$ subject to $x + y = 2$ and $x \le 0.5$. Write the augmented Lagrangian (equality + inequality penalty), set up the KKT system by hand, and verify that the inequality is active with multiplier $\mu \ge 0$.
5. **QP with three equalities.** Min $\tfrac12 x^\top Q x + c^\top x$ subject to $A x = b$ with
   $Q = I_4$, $c = (1, -2, 3, -4)$, and $A x = b$ given by $x_1 + x_2 = 1$, $x_3 + x_4 = 0$, $x_1 - x_3 = 1$. Solve via the KKT matrix and confirm the optimum satisfies the KKT conditions.
6. **Active-set exercise.** Solve by hand the QP $\min\ (x_1-1)^2 + (x_2-2)^2$ subject to $x_1 + 2x_2 \le 2$, $x_1, x_2 \ge 0$. Enumerate the active-set possibilities and show that the optimum lies on $x_1 + 2x_2 = 2$ with $x_2 = 0$.
7. **Solver cross-check.** For any of the above, solve numerically in Python (`scipy.optimize.minimize(method='SLSQP')`) and confirm that the KKT multipliers reported by `method='trust-constr'` match your hand-derived values.
