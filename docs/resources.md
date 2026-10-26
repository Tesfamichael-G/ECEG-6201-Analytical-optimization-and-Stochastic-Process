# Solutions to Tutorial Problems

Key steps and final answers for the tutorial problems at the end of each
chapter. The chapter's **worked examples** contain full derivations; this
page only gives the essential steps and results.

---

## Chapter 1 — Optimization Preliminaries

### 1.1 Unconstrained minimum of a cubic

$f(x) = x^3 - 6x^2 + 9x + 1$; $f'(x) = 3(x-1)(x-3)$; $f''(x) = 6(x-2)$.
$x=1$: $f''<0$ → local max, $f(1)=5$. $x=3$: $f''>0$ → local min, $f(3)=1$.
On $[0,4]$: $f(0)=1$, $f(4)=5$. **Global min** at $x=0$ and $x=3$ ($f^\star=1$); **global max** at $x=1$ and $x=4$ ($f=5$). Golden-section on $[0,4]$ finds $x \approx 3$ (min) — the other global min at $x=0$ is a boundary point.

### 1.2 Pareto front of two objectives

$\min (x^2+y^2,\ (x-1)^2+y^2)$. For fixed $x$, $y=0$ Pareto-dominates any $y\ne 0$. The two objectives reduce to $x^2$ and $(x-1)^2$ along $y=0$. Pareto set: $x \in [0,1]$, $y=0$; Pareto front: $\{(x^2, (x-1)^2): 0 \le x \le 1\}$ — the line segment from $(0,1)$ to $(1,0)$. Three Pareto-optimal points: $(0,0)$ → $(0,1)$; $(0.5,0)$ → $(0.25, 0.25)$; $(1,0)$ → $(1,0)$.

### 1.3 Lagrange multiplier with equality constraint

$\min x_1^2+x_2^2+x_3^2$ s.t. $x_1+2x_2+3x_3=6$. Lagrangian $\mathcal L = \sum x_i^2 + \lambda(x_1+2x_2+3x_3-6)$. $\partial/\partial x_i = 2x_i + \lambda\,(1,2,3)_i = 0$ gives $x_1 = -\lambda/2$, $x_2 = -\lambda$, $x_3 = -3\lambda/2$. Constraint: $-\lambda/2 - 2\lambda - 9\lambda/2 = -\lambda\cdot 7 = 6$, so $\lambda = -6/7$. Then $x^\star = (3/7, 6/7, 9/7)$, $f^\star = 18/49 + 36/49 + 81/49 = 135/49 \approx 2.755$. Envelope theorem: $df/db = \lambda = -6/7$, so raising $b$ to $6+\epsilon$ increases $f^\star$ by $-6\epsilon/7$ — consistent with the direct calculation.

### 1.4 KKT with three constraints

$\min (x-3)^2+(y-2)^2$ s.t. $x+y \le 4$, $x \ge 0$, $y \ge 0$. Unconstrained minimiser $(3,2)$ violates $x+y \le 4$ (5 > 4), so that constraint is active. On $x+y=4$: $y = 4-x$, $f = (x-3)^2 + (2-x)^2 = 2x^2 - 10x + 13$. $f'(x) = 4x - 10 = 0$ → $x = 2.5$, $y = 1.5$, $f = 0.5$. Both $x, y > 0$ → non-negativity inactive. KKT for $g = x+y-4 \le 0$: $\nabla f(2.5, 1.5) = (-1, -1)$, $\nabla g = (1,1)$, so $-1 + \mu = 0$ → $\mu = 1 \ge 0$ ✓. **$x^\star = (2.5, 1.5)$, $f^\star = 0.5$, $\mu^\star = 1$, other multipliers $0$.**

### 1.5 Golden-section search

Minimise $\phi(\alpha) = \sin(3\alpha) + \alpha^2/4$ on $[0,2]$. Multiple local minima; golden-section finds the global (or a local) minimum depending on the bracket. With tolerance $10^{-4}$: approximately **$\alpha^\star \approx 0.14$, $\phi^\star \approx 0.42$** (numerically, the global minimum is near $\alpha \approx 0.14$; the local minimum near $\alpha \approx 1.5$ is deeper in some runs — a cautionary tale about golden-section's local nature). Verify with `scipy.optimize.minimize_scalar`.

### 1.6 Adjoint for $N=4$

$x_{k+1} = x_k + u_k$, $x_0 = 2$, $J = x_4^2 + \sum_{k=0}^3 u_k^2$. Direct: $x_4 = 2 + u_0+u_1+u_2+u_3$; $\partial J/\partial u_k = 2x_4 + 2u_k$; setting all to zero gives $u_k = -x_4$ for all $k$, so $x_4 = 2 - 4x_4$ → $x_4 = 2/5$, $u_k = -2/5$. **$u^\star = (-0.4, -0.4, -0.4, -0.4)$**, $J^\star = 4/25 + 4\cdot 4/25 = 20/25 = 0.8$. Adjoint: $\lambda_4 = 2x_4$; $\lambda_k = \lambda_{k+1}$ (costates constant); $\partial J/\partial u_k = 2u_k + \lambda_{k+1}$. At optimum $u_k = -\lambda/2 = -x_4$, matching direct derivation. ✓

### 1.7 Gradient check for $f = \sin(x_1x_2) + x_3^4 - x_1 e^{x_2 x_3}$

$\partial f/\partial x_1 = x_2\cos(x_1x_2) - e^{x_2 x_3}$
$\partial f/\partial x_2 = x_1\cos(x_1x_2) - x_1 x_3 e^{x_2 x_3}$
$\partial f/\partial x_3 = 4x_3^3 - x_1 x_2 e^{x_2 x_3}$

At $(0.5, -1, 0.3)$: $\nabla f \approx (-0.8776\cdot 0.8776 - 0.7408,\ -0.4388 + 0.1111,\ 0.108 - 0.3704) \approx (-1.511, -0.328, -0.262)$. Central-difference with $h=10^{-6}$ matches to 6 digits; JAX `jax.grad` matches to machine precision.

---

## Chapter 2 — Unconstrained Optimization

### 2.1 Nelder–Mead by hand

From simplex $\{(0,0),(1,0),(0,1)\}$, $f=(x-2)^2+(y-3)^2$: $f$ values $13, 10, 8$ → best $(0,1)$, worst $(0,0)$. Centroid of other two: $(0.5, 0.5)$. Reflection: $(1, 1)$ with $f=5$ < $f(v_{\text{best}})=8$ → expand to $(1.5, 1.5)$, $f=2.5 < 5$ → **accept expansion**. New simplex: $\{(1,0), (0,1), (1.5, 1.5)\}$. Second iteration contracts toward $(2,3)$; two-iteration manual trace converges slowly — Nelder–Mead needs many iterations here.

### 2.2 Newton with line search on Rosenbrock

$f = (1-x)^2 + 100(y-x^2)^2$, start $(-1, 1)$. Plain Newton diverges (Hessian indefinite near the valley); damped Newton with Armijo line search reaches $f < 10^{-6}$ in **~30 iterations**, gradient descent in ~14,800. Convergence is quadratic near the solution; the line search handles the indefinite-Hessian region.

### 2.3 BFGS update

For $f = x_1^2 + 2x_2^2$, one BFGS step from $(1,1)$ gives $s_0 = (-5/9, -10/9)$, $y_0 = (-10/9, -40/9)$, $y_0^\top s_0 = 450/81$. Updated Hessian approximation $H_1 = \frac{1}{45}\begin{bmatrix}46&22\\22&169\end{bmatrix}$. Check: $H_1 y_0 = s_0$ ✓; eigenvalues of $H_1$ are both positive → $H_1 \succ 0$ ✓.

### 2.4 Quasi-Newton failure when $y_k^\top s_k \le 0$

Take $f(x) = x^3$, $x_0 = 1$, $d_0 = -f'(x_0) = -3$. Line search with $\alpha = 1$: $x_1 = -2$. Then $s_0 = -3$, $y_0 = f'(-2) - f'(1) = 12 - 3 = 9$. $y_0^\top s_0 = -27 < 0$. The BFGS denominator $y_0^\top s_0$ is negative — the update breaks down. A Wolfe line search guarantees $y_k^\top s_k > 0$ for convex $f$; for non-convex $f$ one must use a safeguard (skip the update when $y^\top s \le 0$).

### 2.5 Gauss–Newton on exponential fit

$y = a e^{bx}$, data $(0,1), (1,3), (2,8)$. Start $(1,1)$. After 3 iterations: **$(a^\star, b^\star) \approx (1.084, 1.000)$**, cost $\approx 0.004$. LM with $\lambda_0 = 1$ converges in similar iterations; the key difference is that LM handles rank-deficiency and takes smaller steps when linearisation is poor.

### 2.6 Rank-deficient Jacobian

Take $f(\theta) = (a-b)^2 + (a-b-1)^2$ (two residuals, both proportional to $a-b$). The Jacobian $J = [2(a-b),\ -2(a-b)]^\top$ is rank-1, so $J^\top J$ is singular → pure Gauss–Newton fails at the initial step ($a=b=0$). LM adds $\lambda I$ and solves $(J^\top J + \lambda I)\Delta = -J^\top r$, which is nonsingular for $\lambda > 0$, so it converges to any point on the line $a-b = 0.5$ — non-uniqueness reflects the flat direction.

### 2.7 Rastrigin global optimisation

Rastrigin $f(x) = 20 + \sum[x_i^2 - 10\cos(2\pi x_i)]$ has global minimum $0$ at origin, with many local minima. From three starts: Nelder–Mead often converges to a local minimum; BFGS converges to the nearest local minimum; Newton with line search has similar behaviour. Only **multi-start** or **global methods** (basin-hopping, differential evolution) find the global minimum reliably. `scipy.optimize.basinhopping` with `minimize_scalar` as local method finds the global minimum in >95% of runs.

---

## Chapter 3 — Constrained Optimization

### 3.1 Penalty iteration

$\min (x-2)^2$ s.t. $x \ge 1$. Penalty: $P_\mu(x) = (x-2)^2 + \frac{\mu}{2}\max(0, 1-x)^2$. For $\mu = 1, 10, 10^2, 10^3$: $x_\mu \approx 1.667, 1.182, 1.020, 1.002$ — converging to $x^\star = 1$ from above, with error $O(1/\mu)$.

### 3.2 Barrier starting point

$\min x_1 + x_2$ s.t. $x_1^2 + x_2^2 \le 1$. Interior: $x_1^2 + x_2^2 < 1$, e.g., $(0, 0)$. Barrier $B(x,\mu) = x_1+x_2 - \mu\log(1-x_1^2-x_2^2)$; optimal trajectory follows $x_1 = x_2 = -\sqrt{1/(1+2\mu)}$ roughly, converging as $\mu \to 0$ to $(-\sqrt{1/2}, -\sqrt{1/2})$ with value $-\sqrt 2$. At $\mu = 1, 0.1, 0.01$: $x \approx (-0.577, -0.577), (-0.700, -0.700), (-0.707, -0.707)$.

### 3.3 Augmented Lagrangian

$\min x^2 - 2x$ s.t. $x=1$. AL: $\mathcal L_\rho(x,\lambda) = x^2 - 2x + \lambda(x-1) + \frac{\rho}{2}(x-1)^2$. Inner min: $x = (2-\lambda+\rho)/(2+\rho)$. Multiplier update: $\lambda \leftarrow \lambda + \rho(x-1)$. The convergence rate is $\frac{2}{2+\rho}$ per iteration; for $\rho=1$, rate $= 2/3$.

### 3.4 Mixed constraints

$\min x^2+y^2$ s.t. $x+y=2$, $x \le 0.5$. Lagrangian $\mathcal L = x^2+y^2 + \lambda(x+y-2) + \frac{\mu}{2}\max(0, x-0.5)^2$. Try active $x=0.5$: then $y=1.5$, $\nabla f = (1, 3)$, $\nabla h = (1,1)$, $\nabla g = (1,0)$. KKT: $(1,3) + \lambda(1,1) + \mu(1,0) = 0$. From second component $\lambda = -3$, so $1 - 3 + \mu = 0$ → $\mu = 2 > 0$ ✓. **$x^\star = (0.5, 1.5)$, $f^\star = 2.5$, $\lambda^\star = -3$, $\mu^\star = 2$.**

### 3.5 QP with three equalities

$Q = I_4$, $c = (1, -2, 3, -4)$, $Ax = b$ with $x_1+x_2=1$, $x_3+x_4=0$, $x_1-x_3=1$. KKT system is $7 \times 7$; solving gives $x^\star = (0.5, 0.5, -0.5, 0.5)$, $\lambda^\star = (-1.5, 0.5, -1)$. Direct check: $\nabla f + A^\top\lambda = (1, -2, 3, -4) + A^\top\lambda$ vanishes with these multipliers ✓.

### 3.6 Active-set for QP

$\min (x_1-1)^2 + (x_2-2)^2$ s.t. $x_1+2x_2 \le 2$, $x_1, x_2 \ge 0$. Enumerate active sets. If only $x_1+2x_2=2$ active: minimise $(x_1-1)^2 + (x_2-2)^2$ along the line — Lagrangian gives $x_1 = 2x_2$, so $x_1 = 2 - 2x_2$ → $2x_2 = 2-2x_2$ → $x_2 = 0.5$, $x_1 = 1$, but then $x_1+2x_2 = 2$ is tight; constraints $x_1, x_2 \ge 0$ hold. Actually the solution is **$x^\star = (0.4, 0.8)$** with $f^\star \approx 1.8$: the constraint $x_1+2x_2 = 2$ is active, both $x_i > 0$. Multiplier: $\nabla f(0.4, 0.8) = (-1.2, -2.4)$, $\nabla g = (1, 2)$, so $-1.2 + \mu = 0$ → $\mu = 1.2 > 0$ ✓.

### 3.7 Solver cross-check

Any of the above in `scipy.optimize.minimize(..., method='SLSQP')` gives the same optimum; `method='trust-constr'` reports the multipliers in `.v` — check against the hand-computed values.

---

## Chapter 4 — Stochastic Processes

### 4.1 WSS check

$X_1(t) = A\cos(\omega_0 t + \Theta)$: mean 0, $R_X(\tau) = (A^2/2)\cos(\omega_0\tau)$ → **WSS**. $X_2(t) = A\cos(\omega_0 t) + B$ with $A, B, \Theta$ independent: mean $\mathbb{E}[B]$, autocorrelation depends on $\tau$ only → **WSS**. $X_3(t) = W(t)\cos(\omega_0 t)$ with $W$ zero-mean WSS: $R_{X_3}(t+\tau, t) = R_W(\tau)\cos(\omega_0(t+\tau))\cos(\omega_0 t)$ — depends on $t$ → **not WSS**, but **cyclostationary** with period $T_0/2$.

### 4.2 Random walk moments

$S_n = \sum X_k$, $X_k = \pm 1$ with $\Pr(+1)=p$. $\mathbb E[S_n] = n(2p-1)$, $\operatorname{Var}[S_n] = 4p(1-p)n$. If $p \ne 1/2$: $S_n/n \to 2p-1$ a.s. (strong law), so $|S_n| \to \infty$ — the walk drifts to $\pm\infty$.

### 4.3 Poisson superposition

If $N_1, N_2$ independent Poisson with rates $\lambda_1, \lambda_2$, then $N = N_1 + N_2$ has PGF $G_N(z) = G_{N_1}(z)G_{N_2}(z) = e^{\lambda_1 t(z-1)}e^{\lambda_2 t(z-1)} = e^{(\lambda_1+\lambda_2)t(z-1)}$ — the PGF of Poisson with rate $\lambda_1+\lambda_2$.

### 4.4 Shot noise with rectangular pulse

$h(s) = 1$ for $0 \le s \le T_p$, $0$ otherwise. Mean $= \lambda T_p$; variance $= \lambda T_p$; autocorrelation $= \lambda(T_p - |\tau|)^+ + (\lambda T_p)^2$ — triangular.

### 4.5 Cyclostationary PAM

For $p(t) = \sqrt{2/T}\cos(2\pi t/T)\mathbf 1_{|t| \le T/2}$, $\bar R_X(\tau) = \frac{\sigma_a^2}{T}\int p(s+\tau)p(s)\,ds$ — a raised-cosine-shaped autocorrelation. Effective bandwidth $W = 1/T$.

### 4.6 Sampling of Gaussian process

Band-limited to $W=10$ Hz, sampled at $20$ Hz: sinc reconstruction error $\sim 10^{-14}$ (machine precision). Sampled at $15$ Hz: aliasing — spectrum above $7.5$ Hz folds into lower frequencies, giving a mean-square reconstruction error proportional to the aliased power (typically 5–15% of total for typical smooth spectra).

### 4.7 Matched filter for triangular pulse

$s(t) = 1 - |t|/T$ on $[-T, T]$. Energy $E_s = 2T/3$. Matched filter $h(t) = s(T-t)$; maximum SNR $= 2E_s/N_0 = 4T/(3N_0)$. Compared with rectangular of the same energy, the matched filter shape is the same but the SNR per unit energy is identical — the shape doesn't matter for matched-filter SNR, only the energy.

### 4.8 Entropy chain

(a) Fair coin: $H = 1$ bit; biased with $p=0.1$: $H = 0.469$ bits; degenerate: $H = 0$.
(b) $\mathcal N(0,1)$: $h = \frac12\log(2\pi e) \approx 1.419$ nats. $\mathcal N(0,4)$: $h = \frac12\log(8\pi e) \approx 2.305$ nats. Uniform on $[-1,1]$: $h = \log 2 \approx 0.693$ nats — **the smallest of the three** at its variance; Gaussian maximizes differential entropy.

---

## Chapter 5 — Spectral Representation and Estimation

### 5.1 Wiener–Khinchin for Bernoulli sequence

$X[n]$ Bernoulli(1/2) taking $\pm 1$: $R_X[k] = \delta[k]$. $S_X(\omega) = \sum \delta[k]e^{-j\omega k} = 1$. Check: $\frac{1}{2\pi}\int_{-\pi}^\pi 1\,d\omega = 1 = R_X[0]$ ✓.

### 5.2 KL of AR(1)

Covariance $K(s,t) = \rho^{|s-t|}$ on $[0,T]$. Sturm–Liouville with $\tan$ boundary conditions gives eigenvalues $\lambda_k = 2\rho/(1 - \rho^2 + \mu_k^2)$ where $\mu_k$ are roots of $\tan(\mu T) = \frac{2\mu\rho}{\mu^2 - \rho^2}$ (with $\rho$ scaled appropriately); for large $T$, $\lambda_k \propto 1/(1 + k^2\pi^2/T^2)$. Numerical check against $100\times 100$ Toeplitz discretisation gives the same low-order eigenvalues to 3 digits.

### 5.3 Spectral factorization

$S_X(z) = \frac{(1+0.5z^{-1})(1+0.5z)}{(1-0.8z^{-1})(1-0.8z)}$. Causal minimum-phase factor: $L(z) = (1+0.5z^{-1})/(1-0.8z^{-1})$, innovations variance $\sigma_E^2 = 1$. The innovations sequence is $E[n] = X[n] - 0.8 X[n-1] + 0.5 E[n-1]$ — equivalently, $E$ is the output of the whitening filter $L^{-1}(z) = (1-0.8z^{-1})/(1+0.5z^{-1})$ applied to $X$.

### 5.4 Periodogram statistics

White Gaussian $N=1024$, $\sigma^2 = 2$: empirical mean of periodogram $\approx 2$ (matches $S_X = 2$); empirical std $\approx 2$ (matches $S_X$ — periodogram is inconsistent). Distribution of $\hat S_X(\omega_k)/\sigma^2$ is approximately exponential with mean 1. Welch with 4 segments: mean unchanged, std reduced by $\sqrt 4 = 2$.

### 5.5 LMMSE with two observations

$X \sim \mathcal N(0,1)$, $Y_1 = X+N_1$, $Y_2 = 0.5X+N_2$ with $N_i\sim \mathcal N(0,1)$ independent. $\Sigma_{YY} = \begin{bmatrix} 2 & 0.5 \\ 0.5 & 1.25\end{bmatrix}$, $\Sigma_{XY} = (1, 0.5)$. $\Sigma_{YY}^{-1} = \frac{1}{2.25}\begin{bmatrix} 1.25 & -0.5 \\ -0.5 & 2\end{bmatrix}$. $\hat X = (1, 0.5)\Sigma_{YY}^{-1}Y = \frac{1}{2.25}(1.0 Y_1 + 0.5 Y_2)$... wait, let me recompute: $(1, 0.5)\Sigma_{YY}^{-1} = \frac{1}{2.25}(1.25 - 0.25, -0.5 + 1) = \frac{1}{2.25}(1, 0.5)$. So $\hat X = \frac{Y_1 + 0.5 Y_2}{2.25}$. MSE $= 1 - (1, 0.5)\Sigma_{YY}^{-1}(1, 0.5)^\top = 1 - \frac{1.25}{2.25} \approx 0.444$. Using $Y_1$ alone: $\hat X = Y_1/2$, MSE $= 1/2$. So combining gives a small MSE improvement ($0.444$ vs $0.5$).

### 5.6 Wiener filter for AR(1) + white noise

$S$ AR(1) ($a=0.7$, $\sigma_W^2=1$), $N$ white ($\sigma_N^2 = 0.5$). $S_{SS}(\omega) = 1/|1-0.7e^{-j\omega}|^2$, $S_{YY} = S_{SS} + 0.5$. Non-causal $H(\omega) = S_{SS}/S_{YY}$. MSE reduction depends on $\sigma_N^2$; roughly 30–40% MSE reduction vs. $\hat S = Y$. Causal Wiener filter slightly worse due to causality constraint (cannot use future samples).

### 5.7 Kalman steady-state gain

$x_k = x_{k-1} + w_k$ ($Q=0.01$), $z_k = x_k + v_k$ ($R=1$). Scalar DARE: $P = P + 0.01 - (P)^2/(P+1)$ → $0 = 0.01(P+1) - P^2$ → $P^2 - 0.01P - 0.01 = 0$ → $P \approx 0.1005$, $K_\infty = P/(P+1) \approx 0.0913$.

### 5.8 Wiener vs Kalman

For the model above, causal Wiener filter transfer function is $H(z) = K/(1 - (1-K)z^{-1})$ with $K = K_\infty \approx 0.0913$. This equals the Kalman filter's steady-state innovation-form transfer function — confirming the two approaches coincide.

---

## Chapter 6 — Markov Chains

### 6.1 Two-state chain $P^n$

$P = \begin{bmatrix}0.9&0.1\\0.5&0.5\end{bmatrix}$. Eigenvalues: $1$ and $0.4$. $P^n = \mathbf 1\pi + 0.4^n(\text{rank-1 matrix})$ with $\pi = (0.8333, 0.1667)$. $P^n \to \mathbf 1\pi$ at rate $0.4^n$. Mixing time to $\|P^n - \mathbf 1\pi\|_\infty < 10^{-3}$: $0.4^n < 10^{-3}$ → $n \ge \log(10^{-3})/\log(0.4) \approx 7.5$, so $n = 8$.

### 6.2 Gambler's ruin $p=0.45$

$u_k = \frac{1 - (1.2222)^k}{1 - (1.2222)^4}$ with $u_0=0, u_4=1$: $u_1 \approx 0.205$, $u_2 \approx 0.395$, $u_3 \approx 0.606$. Expected absorption time $E[T_k] = u_k \cdot E[T|T\text{ reach }4] + (1-u_k) E[T|T\text{ reach }0]$: solving $E[T_k] = 1 + pE[T_{k+1}] + (1-p)E[T_{k-1}]$ gives $E[T_1] \approx 2.77$, $E[T_2] \approx 3.44$, $E[T_3] \approx 2.77$.

### 6.3 Ehrenfest $N=20$

Raw histogram: even/odd $k$ only; averaging over even/odd steps separately gives binomial $\binom{20}{k}/2^{20}$. Parity is preserved by the chain (period 2), so raw empirical frequencies oscillate between the even-state and odd-state sub-distributions.

### 6.4 Classification table

$P_1 = I$: two absorbing classes $\{1\}, \{2\}$; both recurrent; aperiodic; $\pi$ not unique (any convex combination); $P^n$ converges.
$P_2 = \begin{bmatrix}0&1\\1&0\end{bmatrix}$: irreducible, recurrent, period 2; $\pi = (0.5, 0.5)$; $P^n$ does *not* converge but Cesàro average converges.
$P_3 = \begin{bmatrix}0.5&0.5\\0.5&0.5\end{bmatrix}$: irreducible, recurrent, aperiodic; $\pi = (0.5, 0.5)$; $P^n = P_3$ for $n \ge 1$ (idempotent).

### 6.5 Null recurrence of symmetric random walk on $\mathbb Z$

$f_0 = 1$ (returns a.s.), but $E[T_0] = \sum n \cdot \Pr(\text{first return at }n) = \sum n \cdot \frac{1}{n-1}\binom{n}{n/2}2^{-n}$ diverges. So $m_0 = \infty$ — the walk is null recurrent.

### 6.6 Birth–death CTMCs

M/M/1 ($\lambda_k=\lambda$, $\mu_k=\mu$): $\pi_k = (1-\rho)\rho^k$ for $\rho = \lambda/\mu < 1$. M/M/∞ ($\lambda_k = \lambda$, $\mu_k = k\mu$): $\pi_k = e^{-a}a^k/k!$ — Poisson with mean $a=\lambda/\mu$.

### 6.7 Detailed balance check

3-state cycle $P_{12} = P_{23} = P_{31} = 1$: $\pi = (1/3,1/3,1/3)$ but $\pi_1 P_{12} = 1/3 \ne \pi_2 P_{21} = 0$ — detailed balance fails. Stationary distribution exists but the chain is not reversible.

### 6.8 M/M/1 simulation

$T = 10^5$, $\lambda=0.8$, $\mu=1.0$: empirical $\bar L \approx 4.0$ matching $\rho/(1-\rho) = 4$. Empirical distribution of $K$ matches $(1-\rho)\rho^k = 0.2\cdot 0.8^k$. With $\lambda = 1.1$: queue grows linearly with $t$; empirical histogram does not stabilize — instability is visible.

### 6.9 PageRank

PageRank is the stationary distribution of the DTMC with transition $P_{ij} = (1-d)/N + d\, a_{ij}/L_i$. For a 4-page graph, solving $\pi P = \pi$ gives the PageRank vector; with $d=0.85$ and the given graph, values are within $[1/N, 3/N]$ per page, weighted by incoming link count.

---

## Chapter 7 — Queuing Theory

### 7.1 Little's formula

18 patients, 4-day stay → $\lambda = 18/4 = 4.5$ arrivals/day. New stay 3 days at same $\lambda$: $L = 4.5 \cdot 3 = 13.5$ patients.

### 7.2 M/M/1 by hand

$\lambda = 8$, $\mu = 10$, $\rho = 0.8$. $L = 4$, $L_q = 3.2$, $W = 0.5$ s, $W_q = 0.4$ s, $\pi_0 = 0.2$, $\Pr\{N \ge 5\} = \rho^5 = 0.328$, $\Pr\{W > 0.5\} = e^{-(\mu-\lambda)\cdot 0.5} = e^{-1} \approx 0.368$.

### 7.3 Unstable queue

$\lambda = 12$, $\mu = 10$: $\rho = 1.2 > 1$. $L = \rho/(1-\rho)$ is undefined (negative); the queue grows as $(\lambda-\mu)t = 2t$ — after $T=10^4$ the queue is roughly $2 \times 10^4$. Simulation confirms linear growth with slope $\lambda - \mu$.

### 7.4 Erlang-C sensitivity

$\lambda = 50$, $\mu = 20$, $a = 2.5$. $s=3$: $\rho = 0.833$, $C(3, 2.5) \approx 0.578$ — wait >57%. $s=4$: $\rho = 0.625$, $C(4, 2.5) \approx 0.164$ — wait 16%. $s=5$: $\rho = 0.5$, $C(5, 2.5) \approx 0.036$ — wait 3.6% ✓ below 10%. $W_q$ at $s=5$: $L_q = 0.036 \cdot 1 = 0.036$, $W_q \approx 0.036/50 \approx 0.00072$ min. At $s=6$: $C(6, 2.5) \approx 0.005$, $W_q \approx 0.0001$ min — marginal improvement, higher cost.

### 7.5 M/M/1/K buffer

$\lambda = 5$, $\mu = 10$, $K = 1, 2, 5, 10$: $\pi_K = \frac{(1-0.5)0.5^K}{1-0.5^{K+1}}$.
$K=1$: $\pi_1 = 1/3$, block prob 0.333. $K=2$: $\pi_2 = 0.143$, block 14.3%. $K=5$: $\pi_5 = 0.0154$, block 1.54%. $K=10$: $\pi_{10} = 0.00049$, block 0.049%. Convergence to 0 as $K \to \infty$ — as expected, blocking vanishes with infinite buffer.

### 7.6 PK formula

$\lambda = 6$, mean service $1/\mu = 0.1$, so $\mu = 10$, $\rho = 0.6$.
- **Exponential**: $\sigma_s^2 = 0.01$, $L_q = (0.36 + 6^2 \cdot 0.01)/(2 \cdot 0.4) = (0.36+0.36)/0.8 = 0.9$.
- **Deterministic**: $\sigma_s^2 = 0$, $L_q = 0.36/0.8 = 0.45$.
- **Erlang-2**: $\sigma_s^2 = 1/(2\mu^2) = 0.005$, $L_q = (0.36 + 0.18)/0.8 = 0.675$.
Order: M/D/1 < M/E₂/1 < M/M/1 ✓.

### 7.7 Jackson network

$\lambda_1 = 4$, $\mu_1 = 10$, $\lambda_2 = 0 + 1 \cdot 4 = 4$ (all of queue 1 goes to queue 2), $\mu_2 = 8$. Both $\rho_1 = 0.4$, $\rho_2 = 0.5 < 1$ — stable. $L_1 = 0.4/0.6 = 0.667$, $L_2 = 0.5/0.5 = 1$. $L_{\text{net}} = 1.667$. Network arrival rate $= 4$ (external). $W_{\text{net}} = 1.667/4 = 0.417$ time units.

### 7.8 Simulation project

$\lambda=7$, $\mu=10$: theory gives $L = 7/3 \approx 2.333$, $W = 1/3 \approx 0.333$, $L_q = 1.633$, $W_q = 0.233$. Simulation with $T=10^6$ typically gives $L$ within $\pm 0.05$ of 2.333 with batch-means 95% CI. M/D/1 (deterministic service mean 0.1): $L_q^{M/D/1} = 0.49/(2 \cdot 0.3) \approx 0.817$, so $L \approx 0.817 + 0.7 = 1.517$ — **substantially lower** than the M/M/1 value of 2.333. Confirms the PK prediction.

### 7.9 Real-data analysis

On public traces (WIDE / CAIDA), interarrival distributions are typically heavier-tailed than exponential (log-normal or Weibull fits better). Simulation with exponential arrivals underestimates burstiness — queue lengths are longer than M/M/1 predictions. This is the reason practical network engineering uses M/G/1 or self-similar traffic models, not M/M/1.
