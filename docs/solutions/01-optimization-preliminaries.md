# Solutions to Chapter 1 — Optimization Preliminaries

[:octicons-arrow-left-24: Back to Chapter 1](../chapters/01-optimization-preliminaries.md)

Key steps and final answers. See the corresponding worked examples in the chapter
for full derivations.

---

## 1.1 Unconstrained minimum of a cubic

$f(x) = x^3 - 6x^2 + 9x + 1$ on $[0,4]$.

- $f'(x) = 3(x-1)(x-3)$ → stationary points $x = 1, 3$.
- $f''(x) = 6(x-2)$: $f''(1) = -6$ (local max, $f=5$); $f''(3) = 6$ (local min, $f=1$).
- Boundary values: $f(0) = 1$, $f(4) = 5$.

**Global min:** $x=0$ and $x=3$, $f^\star = 1$. **Global max:** $x=1$ and $x=4$, $f=5$.

Golden-section on $[0,4]$ finds $x \approx 3$ (the interior min); the other global min at $x=0$ is a boundary point that golden-section cannot reach.

---

## 1.2 Pareto front of two objectives

$\min (x^2+y^2,\ (x-1)^2+y^2)$ over $\mathbb{R}^2$.

For any fixed $x$, $y=0$ Pareto-dominates any $y \ne 0$ (both objectives strictly smaller). Along $y=0$:

$$
f_1 = x^2,\quad f_2 = (x-1)^2.
$$

For $x<0$: moving to $x=0$ reduces both. For $x>1$: moving to $x=1$ reduces both. So Pareto set $= \{(x,0): 0 \le x \le 1\}$, Pareto front $= \{(x^2, (x-1)^2): 0 \le x \le 1\}$ — a segment from $(0,1)$ to $(1,0)$.

**Three Pareto-optimal points:** $(0,0) \to (0,1)$; $(0.5,0) \to (0.25, 0.25)$; $(1,0) \to (1,0)$.

---

## 1.3 Lagrange multiplier with equality constraint

$\min x_1^2+x_2^2+x_3^2$ s.t. $x_1+2x_2+3x_3=6$.

Lagrangian $\mathcal L = \sum x_i^2 + \lambda(x_1+2x_2+3x_3-6)$. From $\partial/\partial x_i = 0$:

$$
x_1 = -\tfrac{\lambda}{2},\quad x_2 = -\lambda,\quad x_3 = -\tfrac{3\lambda}{2}.
$$

Constraint: $-\tfrac{\lambda}{2} - 2\lambda - \tfrac{9\lambda}{2} = -7\lambda = 6 \Rightarrow \lambda = -\tfrac{6}{7}$.

**$x^\star = \bigl(\tfrac{3}{7}, \tfrac{6}{7}, \tfrac{9}{7}\bigr)$, $f^\star = \tfrac{135}{49} \approx 2.755$, $\lambda^\star = -\tfrac{6}{7}$.**

Envelope theorem check: $df^\star/db = \lambda = -6/7$ — raising $b$ by $\epsilon$ reduces $f^\star$ by $\epsilon \cdot 6/7$.

---

## 1.4 KKT with three constraints

$\min (x-3)^2+(y-2)^2$ s.t. $x+y \le 4$, $x \ge 0$, $y \ge 0$.

Unconstrained minimiser $(3,2)$ violates $x+y \le 4$ (5 > 4). On $x+y = 4$: $y = 4-x$,

$$
f(x) = (x-3)^2 + (2-x)^2 = 2x^2 - 10x + 13.
$$

$f'(x) = 4x - 10 = 0 \Rightarrow x = 2.5$, $y = 1.5$, $f = 0.5$. Both $x, y > 0$ → non-negativity inactive. KKT: $\nabla f(2.5, 1.5) = (-1,-1)$, $\nabla g = (1,1)$ → $-1 + \mu = 0$ → $\mu = 1 \ge 0$ ✓.

**$x^\star = (2.5, 1.5)$, $f^\star = 0.5$, $\mu^\star = 1$, other multipliers $0$.**

---

## 1.5 Golden-section search

$\phi(\alpha) = \sin(3\alpha) + \alpha^2/4$ on $[0,2]$, tolerance $10^{-4}$.

Numerically, the global minimum is **$\alpha^\star \approx 0.14$, $\phi^\star \approx 0.42$**. But golden-section is a *local* method — depending on the bracket, it can converge to the local minimum near $\alpha \approx 1.5$, which is deeper in some runs. This is a cautionary example: **a global-scope claim for golden-section is wrong**; verify against `scipy.optimize.minimize_scalar(method='bounded')` and, if needed, multi-start.

---

## 1.6 Adjoint for $N=4$

$x_{k+1} = x_k + u_k$, $x_0 = 2$, $J = x_4^2 + \sum_{k=0}^3 u_k^2$.

**Direct.** $x_4 = 2 + \sum u_k$; $\partial J/\partial u_k = 2x_4 + 2u_k = 0$ ⇒ $u_k = -x_4$ for all $k$. Then $x_4 = 2 - 4x_4$ ⇒ $x_4 = 2/5$.

**$u^\star = (-0.4, -0.4, -0.4, -0.4)$, $J^\star = 0.8$.**

**Adjoint check.** $\lambda_4 = 2x_4$; $\lambda_k = \lambda_{k+1}$ (costates constant); $\partial J/\partial u_k = 2u_k + \lambda_{k+1}$. At optimum $u_k = -\lambda/2 = -x_4$ — identical to direct. ✓

---

## 1.7 Gradient check

$f = \sin(x_1x_2) + x_3^4 - x_1 e^{x_2 x_3}$.

$$
\frac{\partial f}{\partial x_1} = x_2\cos(x_1x_2) - e^{x_2 x_3},\quad
\frac{\partial f}{\partial x_2} = x_1\cos(x_1x_2) - x_1 x_3 e^{x_2 x_3},\quad
\frac{\partial f}{\partial x_3} = 4x_3^3 - x_1 x_2 e^{x_2 x_3}.
$$

At $(0.5, -1, 0.3)$: $\nabla f \approx (-1.511,\ -0.328,\ -0.262)$.

Central differences with $h=10^{-6}$ match to 6 digits. `jax.grad` matches to machine precision.