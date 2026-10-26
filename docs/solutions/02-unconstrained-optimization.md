# Solutions to Chapter 2 — Unconstrained Optimization

[:octicons-arrow-left-24: Back to Chapter 2](../chapters/02-unconstrained-optimization.md)

---

## 2.1 Nelder–Mead by hand

$f(x,y) = (x-2)^2+(y-3)^2$ from simplex $\{(0,0),(1,0),(0,1)\}$.

**Iteration 1.** $f$-values $13, 10, 8$ → best $(0,1)$, worst $(0,0)$. Centroid of the other two: $(0.5, 0.5)$. Reflection: $(1, 1)$, $f=5 < f(v_{\text{best}})=8$ → expand to $(1.5, 1.5)$, $f=2.5 < 5$ → accept expansion. New simplex: $\{(1,0), (0,1), (1.5, 1.5)\}$.

**Iteration 2.** Values $10, 8, 2.5$ → best $(1.5, 1.5)$, worst $(1, 0)$. Centroid $(0.75, 1.25)$. Reflection: $(0.5, 2.5)$, $f=2.5 < 2.5$? No — $\ge$ the best, so contract: $c = (0.5, 2.5) \to$ midpoint $(1, 1.75)$, $f=2.5625$. Compare to $f(v_{\text{worst}})=10$: accept contract.

Nelder–Mead converges very slowly here — many iterations are needed to get close to $(2,3)$.

---

## 2.2 Newton with line search on Rosenbrock

$f = (1-x)^2 + 100(y-x^2)^2$, start $(-1, 1)$.

- **Plain Newton** diverges (indefinite Hessian near the valley).
- **Damped Newton** (Armijo, $c_1 = 10^{-4}$) reaches $f < 10^{-6}$ in **~30 iterations**.
- **Gradient descent** (Armijo) reaches $f < 10^{-6}$ in ~14,800 iterations.

**Takeaway:** the line search is essential for robustness; once close to the solution, Newton's quadratic convergence dominates.

---

## 2.3 BFGS update

$f = x_1^2 + 2x_2^2$, $x_0 = (1,1)$, $H_0 = I$, exact line search.

From the worked example: $s_0 = \left(-\tfrac{5}{9}, -\tfrac{10}{9}\right)$, $y_0 = \left(-\tfrac{10}{9}, -\tfrac{40}{9}\right)$, $y_0^\top s_0 = \tfrac{450}{81}$.

$$
H_1 = \frac{1}{45}\begin{bmatrix} 46 & 22 \\ 22 & 169\end{bmatrix}.
$$

**Checks.**

- Secant: $H_1 y_0 = s_0$ ✓
- Positive-definite: eigenvalues of $H_1$ are both positive ⇒ $H_1 \succ 0$ ✓

---

## 2.4 Quasi-Newton failure when $y_k^\top s_k \le 0$

$f(x) = x^3$, $x_0 = 1$, $d_0 = -f'(x_0) = -3$, $\alpha_0 = 1$.

$x_1 = -2$; $s_0 = -3$; $y_0 = f'(-2) - f'(1) = 12 - 3 = 9$; $y_0^\top s_0 = -27 < 0$.

The BFGS denominator $y_0^\top s_0 < 0$ → the update loses positive-definiteness (and can even break down). A **Wolfe line search** guarantees $y_k^\top s_k > 0$ on convex $f$; for non-convex $f$, skip the update when $y_k^\top s_k \le 0$.

---

## 2.5 Gauss–Newton on exponential fit

$y = a e^{bx}$, data $(0,1), (1,3), (2,8)$, start $(1,1)$.

After **3 iterations**: $(a^\star, b^\star) \approx (1.084, 1.000)$, cost $\approx 0.004$.

**Why so fast.** The residuals at the solution are small (the exponential fits the data almost exactly), so the dropped residual-curvature term $\sum r_i \nabla^2 r_i$ is negligible — GN behaves like Newton and converges quadratically.

---

## 2.6 Rank-deficient Jacobian

$f(\theta) = (a-b)^2 + (a-b-1)^2$ (two residuals, both proportional to $a-b$). The Jacobian $J$ has rank 1 at every point, so $J^\top J$ is singular → pure Gauss–Newton fails at the first step ($a=b=0$).

**Levenberg–Marquardt** adds $\lambda I$:
$$
(J^\top J + \lambda I)\,\Delta = -J^\top r
$$
is nonsingular for any $\lambda > 0$. LM converges to **any** point on the line $a-b = 0.5$ — the non-uniqueness reflects the flat direction of the objective.

---

## 2.7 Rastrigin global optimisation

$f(x) = 20 + \sum_{i}[x_i^2 - 10\cos(2\pi x_i)]$ (2-D), global min $0$ at origin, many local minima.

- **Nelder–Mead, BFGS, Newton**: from three different starts, each converges to the nearest local minimum (usually not the global).
- **`scipy.optimize.basinhopping`** with `minimize_scalar` as the local method: finds the global minimum in **>95% of runs** on 2-D Rastrigin.

**Takeaway:** for multi-modal objectives, local methods fail regardless of how sophisticated they are. Use a global-scope method (basin-hopping, differential evolution, simulated annealing) or multi-start.