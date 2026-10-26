# Solutions to Chapter 3 — Constrained Optimization

[:octicons-arrow-left-24: Back to Chapter 3](../chapters/03-constrained-optimization.md)

---

## 3.1 Penalty iteration

$\min (x-2)^2$ s.t. $x \ge 1$. Penalty $P_\mu(x) = (x-2)^2 + \tfrac{\mu}{2}\max(0, 1-x)^2$.

From the worked example, $x_\mu = (4+\mu)/(2+\mu)$. As $\mu \to \infty$, $x_\mu \to 1$ with rate $O(1/\mu)$.

| $\mu$ | $x_\mu$ | $\|x_\mu - 1\|$ |
|:---:|:---:|:---:|
| 1 | 1.667 | 0.667 |
| 10 | 1.182 | 0.182 |
| $10^2$ | 1.020 | 0.020 |
| $10^3$ | 1.002 | 0.002 |

---

## 3.2 Barrier starting point

$\min x_1+x_2$ s.t. $x_1^2+x_2^2 \le 1$.

A strictly feasible starting point: $(0,0)$ (inside the unit disk).

Barrier: $B(x,\mu) = x_1+x_2 - \mu\log(1-x_1^2-x_2^2)$.

**First-order conditions** (by symmetry $x_1=x_2=:t$):

$$
1 - \mu\cdot\frac{2t}{1-2t^2} = 0 \Rightarrow 1-2t^2 = 2\mu t \Rightarrow t^2 + \mu t - \tfrac12 = 0
$$

$$
t = \frac{-\mu + \sqrt{\mu^2 + 2}}{2}.
$$

| $\mu$ | $t$ | $(x_1, x_2)$ | $f$ |
|:---:|:---:|:---:|:---:|
| 1 | $-0.366$ | $(-0.366, -0.366)$ | $-0.732$ |
| 0.1 | $-0.606$ | $(-0.606, -0.606)$ | $-1.212$ |
| 0.01 | $-0.700$ | $(-0.700, -0.700)$ | $-1.400$ |

As $\mu \to 0$: $t \to -\sqrt{1/2}$, $f \to -\sqrt 2 \approx -1.414$. The constrained optimum is indeed $(-\sqrt{1/2}, -\sqrt{1/2})$.

---

## 3.3 Augmented Lagrangian

$\min x^2-2x$ s.t. $x = 1$.

$\mathcal L_\rho(x,\lambda) = x^2-2x + \lambda(x-1) + \tfrac{\rho}{2}(x-1)^2$.

**Inner minimisation.** $\partial/\partial x = 2x-2+\lambda+\rho(x-1) = 0$ ⇒ $x = (2-\lambda+\rho)/(2+\rho)$.

**Update.** $\lambda \leftarrow \lambda + \rho(x-1)$.

The convergence rate is $\frac{2}{2+\rho}$ per iteration: for $\rho=1$ the rate is $2/3$.

---

## 3.4 Mixed constraints

$\min x^2+y^2$ s.t. $x+y=2$, $x \le 0.5$.

Try active $x = 0.5$: then $y = 1.5$. $\nabla f = (1, 3)$, $\nabla h = (1,1)$, $\nabla g = (1,0)$.

KKT: $(1,3) + \lambda(1,1) + \mu(1,0) = 0$. From the second component, $\lambda = -3$; then $1 - 3 + \mu = 0 \Rightarrow \mu = 2 > 0$ ✓.

**$x^\star = (0.5, 1.5)$, $f^\star = 2.5$, $\lambda^\star = -3$, $\mu^\star = 2$.**

---

## 3.5 QP with three equalities

$Q = I_4$, $c = (1,-2,3,-4)$, constraints $x_1+x_2=1$, $x_3+x_4=0$, $x_1-x_3=1$.

KKT system is $7 \times 7$. Solving:

**$x^\star = (0.5, 0.5, -0.5, 0.5)$, $\lambda^\star = (-1.5, 0.5, -1)$.**

Check: $\nabla f + A^\top\lambda = (1,-2,3,-4) + A^\top\lambda = 0$ with these multipliers ✓.

---

## 3.6 Active-set for QP

$\min (x_1-1)^2+(x_2-2)^2$ s.t. $x_1+2x_2 \le 2$, $x_1, x_2 \ge 0$.

**Try $x_1+2x_2 = 2$ active, both $x_i > 0$.** $x_1 = 2 - 2x_2$. $f(x_2) = (1-2x_2)^2 + (x_2-2)^2 = 5x_2^2 - 8x_2 + 5$. $f' = 10x_2 - 8 = 0$ ⇒ $x_2 = 0.8$, $x_1 = 0.4$, $f = 1.8$. Both $x_i > 0$ ✓.

**KKT.** $\nabla f(0.4, 0.8) = (-1.2, -2.4)$, $\nabla g = (1, 2)$. $-1.2 + \mu = 0$ ⇒ $\mu = 1.2 > 0$ ✓.

**$x^\star = (0.4, 0.8)$, $f^\star = 1.8$, $\mu^\star = 1.2$.**

---

## 3.7 Solver cross-check

`scipy.optimize.minimize(..., method='SLSQP')` gives the same $x^\star$. `method='trust-constr'` returns the multipliers in the result object (`.v` for inequality, `.multipliers` in newer SciPy), which match the hand-computed values in 3.4, 3.5, and 3.6.