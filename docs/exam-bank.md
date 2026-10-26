# Exam Bank — ECEG-6201 Analytical Optimization and Stochastic Processes

40 questions spanning all seven chapters. Use for self-assessment and exam
preparation. **Answers and grading rubric at the end.**

**Instructions (for a mock exam):**

- Time: **3 hours**
- Closed book; **one A4 handwritten formula sheet** allowed
- Non-programmable calculator permitted
- Total marks: **100**
- Marks per question in brackets

---

## Section A — Optimization Preliminaries (Ch 1) — 15 marks

### A1. Multiple choice [2]

Which of the following is a **necessary** condition for a differentiable $f$
to have a local minimum at an interior point $x^\star$?

(a) $f(x^\star) = 0$
(b) $\nabla f(x^\star) = 0$
(c) $\nabla^2 f(x^\star) \succ 0$
(d) $f$ is convex

### A2. Multiple choice [2]

In multi-objective optimization, a point $x^\star$ is **Pareto-optimal** if:

(a) It minimizes a weighted sum of objectives
(b) No feasible point improves some objective without worsening another
(c) It is a local minimum of each individual objective
(d) It satisfies the KKT conditions

### A3. Short answer [3]

Define the **gradient** and **Hessian** for a function $f: \mathbb R^n \to \mathbb R$.
State the geometric interpretation of the Hessian's eigenvalues.

### A4. Short answer [3]

State the **first-order necessary condition** (FONC) and the **second-order
sufficient condition** (SOSC) for a strict local minimum. Give an example of a
function that satisfies FONC but not SOSC at a point that is not a minimum.

### A5. Derivation [5]

Starting from the **line-search subproblem** $\min_{\alpha > 0} \phi(\alpha) = f(x_k + \alpha d_k)$,
derive the **golden-section search** update rule. State the shrinkage factor
per iteration and the number of iterations needed to reach tolerance $\epsilon$
from initial bracket length $L_0$.

---

## Section B — Unconstrained Optimization (Ch 2) — 15 marks

### B1. Multiple choice [2]

Newton's method converges in how many iterations on a strictly convex quadratic?

(a) 1 (b) 2 (c) $n$ (d) arbitrarily many

### B2. Multiple choice [2]

The **secant condition** for the BFGS Hessian approximation is:

(a) $B_{k+1} y_k = s_k$ (b) $B_{k+1} s_k = y_k$ (c) $B_{k+1} = I$ (d) $B_{k+1} = \nabla^2 f$

### B3. Short answer [4]

Describe the **Nelder–Mead simplex algorithm**. State its reflection,
expansion, contraction, and shrink steps with their standard parameters.

### B4. Derivation [4]

Derive the **Newton step** for unconstrained minimization from the second-order
Taylor expansion. State the local convergence rate and one reason plain Newton
can fail.

### B5. Short answer [3]

State the **Gauss–Newton** step for nonlinear least squares. Under what
condition does it achieve quadratic convergence?

---

## Section C — Constrained Optimization (Ch 3) — 15 marks

### C1. Multiple choice [2]

Which method requires a **strictly feasible** starting point?

(a) Exterior penalty (b) Logarithmic barrier
(c) Augmented Lagrangian (d) QP active-set

### C2. Multiple choice [2]

In the exterior penalty method, the penalty parameter $\mu$ must tend to:

(a) $0^+$ (b) $1$ (c) $+\infty$ (d) any value

### C3. Short answer [3]

State the **KKT conditions** for $\min f(x)$ subject to $g_i(x) \le 0$,
$h_j(x) = 0$. Include complementary slackness.

### C4. Derivation [5]

For a quadratic program $\min \tfrac12 x^\top Q x + c^\top x$ subject to
$Ax = b$, derive the **KKT linear system**. State the existence and
uniqueness condition.

### C5. Short answer [3]

Compare the **augmented Lagrangian** and the **exterior penalty method**.
Why does the augmented Lagrangian avoid the ill-conditioning of pure penalties?

---

## Section D — Stochastic Processes & Random Walks (Ch 4) — 15 marks

### D1. Multiple choice [2]

A **WSS** process has:

(a) Constant variance, arbitrary mean
(b) Constant mean, autocorrelation depending only on lag
(c) Gaussian distribution
(d) Finite energy

### D2. Multiple choice [2]

The variance of a shot-noise process $\sum_i h(t - t_i)$ (Poisson arrivals, rate $\lambda$) is:

(a) $\lambda \int h(s)\, ds$
(b) $\lambda \int h^2(s)\, ds$
(c) $\lambda^2 (\int h)^2$
(d) $0$

### D3. Short answer [3]

State **Campbell's theorem**. Give two engineering applications of shot noise.

### D4. Derivation [4]

For $X(t) = \sum_k a_k\, p(t - kT)$ with i.i.d. zero-mean $a_k$, derive
$R_X(t + \tau, t)$. Show that $X$ is cyclostationary with period $T$.

### D5. Short answer [4]

Derive the **matched filter** for detecting a known signal $s(t)$ in white
Gaussian noise. State the maximum output SNR.

---

## Section E — Spectral Representation & Estimation (Ch 5) — 15 marks

### E1. Multiple choice [2]

The **Wiener–Khinchine theorem** states that the PSD is:

(a) The variance of the autocorrelation
(b) The Fourier transform of the autocorrelation
(c) The reciprocal of the autocorrelation
(d) The mean of the signal

### E2. Multiple choice [2]

The **periodogram** of white noise is:

(a) A flat line
(b) Unbiased but has variance equal to its mean
(c) Consistent
(d) Zero at all frequencies

### E3. Short answer [3]

State the **Spectral Factorization Theorem**. Why must the factor $L(z)$ be
**causal, monic, and minimum-phase**?

### E4. Derivation [4]

For an observed process $Y = S + N$ with $S \perp N$, derive the **non-causal
Wiener filter** $H_{\text{opt}}(\omega)$ that minimizes the mean-square error.

### E5. Short answer [4]

State the **Kalman filter predict and update equations**. Interpret the
Kalman gain $K_k$ in words.

---

## Section F — Markov Chains (Ch 6) — 13 marks

### F1. Multiple choice [2]

A probability vector $\pi$ is a **stationary distribution** of a DTMC with
transition matrix $P$ if:

(a) $\pi P = \pi$ and $\sum \pi_i = 1$
(b) $P\pi = \pi$ and $\sum P = 1$
(c) $\pi = P^T$
(d) $P = \pi \pi^\top$

### F2. Short answer [4]

Classify states as **recurrent / transient** and **periodic / aperiodic**.
State the **ergodic theorem** for finite irreducible aperiodic chains.

### F3. Derivation [4]

For a birth–death DTMC with birth probabilities $p_k$ and death probabilities
$q_k$, use **detailed balance** to derive the stationary distribution $\pi_k$.
State the positive-recurrence condition.

### F4. Short answer [3]

Write the **generator matrix** $Q$ for a CTMC. State the **Kolmogorov forward
equation** and how the stationary distribution is obtained from $Q$.

---

## Section G — Queuing Theory (Ch 7) — 12 marks

### G1. Multiple choice [2]

Little's formula $L = \lambda W$ holds:

(a) Only for M/M/1 (b) Only for Poisson arrivals
(c) For any stable queuing system (d) Never

### G2. Short answer [3]

State the **stability condition** for an M/M/1 queue. What happens to the
stationary distribution when the condition is violated?

### G3. Derivation [4]

For an M/M/1 queue with arrival rate $\lambda$ and service rate $\mu$,
derive the **stationary distribution** $\pi_n$ and the mean number in
system $L$. Assume $\lambda < \mu$.

### G4. Short answer [3]

State the **Erlang-C formula**. Give one practical application in call-centre
staffing.

---

## Answer key

??? success "Click to reveal answers"

    ### Section A

    **A1.** (b) — $\nabla f(x^\star) = 0$ is necessary; SOSC $\succ 0$ is sufficient.

    **A2.** (b) — Pareto-optimality: no feasible point Pareto-dominates it.

    **A3.** Gradient: $\nabla f = (\partial f/\partial x_1, \ldots, \partial f/\partial x_n)^\top$. Hessian: $\nabla^2 f = [\partial^2 f/\partial x_i \partial x_j]_{i,j}$. Geometric: eigenvalues of the Hessian are the **principal curvatures** of $f$ at that point; eigenvectors are the directions of maximum and minimum curvature. Positive eigenvalues ⇒ bowl up (convex); negative ⇒ bowl down; mixed ⇒ saddle.

    **A4.** **FONC:** $\nabla f(x^\star) = 0$ if $x^\star$ is a local min and $f$ is differentiable. **SOSC:** $\nabla f(x^\star) = 0$ and $\nabla^2 f(x^\star) \succ 0$ ⇒ $x^\star$ is a **strict** local min. **Counterexample:** $f(x) = x^3$ at $x = 0$ has $f'(0) = 0$ but $f''(0) = 0$ — the point is an inflection, not a minimum.

    **A5.** Golden-section places two interior points $c = b - \tau(b-a)$ and $d = a + \tau(b-a)$ with $\tau = (\sqrt 5 - 1)/2 \approx 0.618$. Evaluate $\phi(c), \phi(d)$; keep the subinterval that contains the smaller value. Each iteration reduces the bracket by $\tau$ with **one new evaluation**. Number of iterations to reach tolerance $\epsilon$ from initial bracket $L_0$: $k = \lceil \log(\epsilon / L_0)/\log\tau \rceil$.

    ### Section B

    **B1.** (a) — one step for a strictly convex quadratic.

    **B2.** (b) — $B_{k+1} s_k = y_k$ with $s_k = x_{k+1} - x_k$, $y_k = \nabla f(x_{k+1}) - \nabla f(x_k)$.

    **B3.** **Nelder–Mead:** maintain a simplex of $n+1$ points. Sort by $f$ value. Compute centroid $\bar v$ of all but worst. **Reflection** $\alpha = 1$: $r = \bar v + (\bar v - v_{\text{worst}})$. **Expansion** $\gamma = 2$: $e = \bar v + 2(r - \bar v)$ if reflection is better than the best. **Contraction** $\rho = \tfrac12$: $c = \bar v + \tfrac12(v_{\text{worst}} - \bar v)$ if reflection is worse than second-worst. **Shrink** $\sigma = \tfrac12$: $v_i \leftarrow v_{\text{best}} + \tfrac12(v_i - v_{\text{best}})$ if contraction also fails.

    **B4.** Taylor: $f(x_k + d) \approx f(x_k) + \nabla f(x_k)^\top d + \tfrac12 d^\top \nabla^2 f(x_k)\, d$. Set gradient in $d$ to zero: $\nabla f(x_k) + \nabla^2 f(x_k)\, d = 0$, giving $d_k = -(\nabla^2 f(x_k))^{-1}\nabla f(x_k)$. **Convergence:** local quadratic if Hessian at the solution is positive-definite. **Fails:** when the Hessian is singular (no unique step) or indefinite (not a descent direction); a line search or regularization $H + \mu I$ fixes this.

    **B5.** For least squares $\min \tfrac12 \|r(\theta)\|^2$, the Gauss–Newton step solves $(J^\top J)\Delta = -J^\top r$, where $J = \partial r/\partial \theta$. **Quadratic convergence** iff the residuals at the solution are zero (dropped term $\sum r_i \nabla^2 r_i = 0$); otherwise only **linear**.

    ### Section C

    **C1.** (b) — logarithmic barrier requires strictly feasible interior point.

    **C2.** (c) — $\mu \to \infty$.

    **C3.** KKT for $\min f$ s.t. $g_i \le 0$, $h_j = 0$:
    (i) $\nabla f + \sum_i \mu_i \nabla g_i + \sum_j \lambda_j \nabla h_j = 0$;
    (ii) $g_i(x) \le 0$, $h_j(x) = 0$ (feasibility);
    (iii) $\mu_i \ge 0$ (dual feasibility);
    (iv) $\mu_i g_i(x) = 0$ for each $i$ (complementary slackness).

    **C4.** Lagrangian $\mathcal L = \tfrac12 x^\top Q x + c^\top x + \lambda^\top(Ax - b)$. Gradient in $x$: $Qx + c + A^\top\lambda = 0$. Gradient in $\lambda$: $Ax - b = 0$. Combine: $\begin{bmatrix}Q & A^\top \\ A & 0\end{bmatrix}\begin{bmatrix}x \\ \lambda\end{bmatrix} = \begin{bmatrix}-c \\ b\end{bmatrix}$. **Existence & uniqueness:** $Q \succ 0$ and $A$ full row rank.

    **C5.** Exterior penalty drives $\mu \to \infty$, which makes the Hessian eigenvalues grow like $\mu$ in constraint-normal directions — condition number grows like $\mu$; simple solvers fail. **Augmented Lagrangian** $\mathcal L_\rho = f + \lambda^\top h + \tfrac{\rho}{2}\|h\|^2$ uses a multiplier update $\lambda \leftarrow \lambda + \rho h(x)$. When $\lambda$ converges to the true KKT multiplier $\lambda^\star$, the constrained optimum is the exact minimizer of $\mathcal L_\rho$ for **any fixed** $\rho > 0$. No need to drive $\rho \to \infty$ ⇒ no ill-conditioning.

    ### Section D

    **D1.** (b).

    **D2.** (b) — $\lambda\int h^2(s)\,ds$.

    **D3.** Campbell's theorem: for Poisson arrivals rate $\lambda$ and response $h$, $\mathbb E[X(t)] = \lambda\int h$, $\operatorname{Var}[X(t)] = \lambda\int h^2$, $R_X(\tau) = \lambda\int h(s) h(s+\tau)\,ds + \lambda^2 (\int h)^2$. Applications: photodetector current (photon arrivals), neural spike trains (synaptic filtering), interference modeling in wireless.

    **D4.** $R_X(t+\tau,t) = \sum_{k,m}\mathbb E[a_k a_m] p(t+\tau-kT)p(t-mT) = \sigma_a^2\sum_k p(t+\tau-kT)p(t-kT)$. The right side is periodic in $t$ with period $T$ — so $X$ is **cyclostationary**. The time-averaged autocorrelation $\bar R_X(\tau) = \frac{\sigma_a^2}{T}\int p(s+\tau)p(s)\,ds$ is the WSS-equivalent.

    **D5.** Output SNR at time $T$: $\text{SNR} = |\int s(\tau) h(T-\tau)\,d\tau|^2 / \left[(N_0/2)\int h^2(u)\,du\right]$. By Cauchy–Schwarz, $\text{SNR} \le 2E_s/N_0$ with equality iff $h(t) \propto s(T-t)$. **Matched filter:** $h(t) = s(T-t)$, maximum output SNR $= 2 E_s / N_0$.

    ### Section E

    **E1.** (b).

    **E2.** (b).

    **E3.** Any rational PSD with no zeros on the unit circle and $\int\log S > -\infty$ can be factored uniquely as $S_X(z) = \sigma_e^2 L(z)L(z^{-1})$ with $L$ causal, monic, minimum-phase. **Why:** (i) causality makes $1/L$ causal (whitening realizable); (ii) monic fixes the scale so $\sigma_e^2$ is uniquely the innovations variance; (iii) minimum-phase ensures smallest $\sigma_e^2$ — any non-minimum-phase factor would have larger innovation energy.

    **E4.** $H_{\text{opt}}(\omega) = S_{SY}(\omega)/S_{YY}(\omega) = S_{SS}/(S_{SS} + S_{NN})$. High SNR → $H \approx 1$ (pass); low SNR → $H \approx 0$ (reject); equal → $H = 0.5$. The filter shapes the passband by the local SNR.

    **E5.** Predict: $\hat x_{k|k-1} = A\hat x_{k-1|k-1}$, $P_{k|k-1} = AP_{k-1|k-1}A^\top + Q$. Gain: $K_k = P_{k|k-1}H^\top(HP_{k|k-1}H^\top + R)^{-1}$. Update: $\hat x_{k|k} = \hat x_{k|k-1} + K_k(z_k - H\hat x_{k|k-1})$, $P_{k|k} = (I - K_kH)P_{k|k-1}$. $K_k$ is the **weight on the new information** — large when measurements are reliable (small $R$) or prediction is uncertain (large $P$), small in the opposite case.

    ### Section F

    **F1.** (a).

    **F2.** *Recurrent:* return probability $f_i = 1$. *Transient:* $f_i < 1$. *Periodic:* period $d(i) = \gcd\{n \ge 1 : P^n_{ii} > 0\} > 1$; aperiodic if period 1. **Ergodic theorem:** for a finite irreducible aperiodic DTMC with stationary $\pi$, $\lim_{n\to\infty} P^n = \mathbf 1\pi$ and $\frac{1}{N}\sum_{n=1}^N f(X_n) \to \sum_j \pi_j f(j)$ a.s.

    **F3.** Detailed balance for birth–death chain: $\pi_k p_k = \pi_{k+1} q_{k+1}$, so $\pi_k = \pi_0\prod_{j=0}^{k-1} p_j/q_{j+1}$. Normalize: $\pi_0 = 1/(1 + \sum_{k=1}^\infty \prod_{j=0}^{k-1} p_j/q_{j+1})$. **Positive recurrence:** the infinite sum converges.

    **F4.** Generator $Q$: $Q_{ij} \ge 0$ for $i \ne j$ (transition rate $i \to j$), $Q_{ii} = -\sum_{j \ne i} Q_{ij}$ (so each row sums to 0). Kolmogorov forward: $\dot P(t) = P(t) Q$, solution $P(t) = e^{Qt}$. Stationary distribution solves $\pi Q = 0$ with $\sum \pi_i = 1$.

    ### Section G

    **G1.** (c).

    **G2.** $\rho = \lambda/\mu < 1$. If $\rho \ge 1$, the queue is **not positive recurrent** — no stationary distribution; queue length grows unboundedly.

    **G3.** Detailed balance: $\pi_n \lambda = \pi_{n+1}\mu$, so $\pi_n = \pi_0 \rho^n$ with $\rho = \lambda/\mu$. Normalize: $\pi_0 = 1 - \rho$. **Mean:** $L = \sum n\pi_n = \rho/(1-\rho)$.

    **G4.** $C(s, a) = \Pr\{N \ge s\}$ for an $s$-server queue with offered load $a = \lambda/\mu$: $C(s,a) = \pi_0 a^s / (s!(1 - \rho))$ with $\rho = a/s$. **Application:** given target wait probability (say 20%), find the minimum $s$ such that $C(s,a) \le 0.20$ — the standard staffing formula used in call centres.

---

## Grading rubric (for instructors)

| Band | Marks | Criteria |
|---|---|---|
| **A** | 90–100 | All correct or with trivial slips; rigorous derivations; clear interpretation |
| **B** | 75–89 | Most correct; occasional gaps; derivations mostly sound |
| **C** | 60–74 | Basic understanding; some errors in derivations or interpretation |
| **F** | < 60 | Major conceptual gaps; significant derivations missing |

**Point allocation.** MCQ all-or-nothing. Short answers partial credit per correct step. Derivations partial credit per correct line.

**Common error checks (and deductions):**

- **Wrong sign** in KKT or Lagrange multiplier: −1 to −2.
- **FONC/SOSC confusion** (claiming FONC ⇒ minimum): −2.
- **Golden-section factor** wrong (0.5 instead of 0.618): −1.
- **Confusing WSS and stationarity** conditions: −1.
- **Wiener–Khinchine sign convention** (using $+j\omega$ instead of $-j\omega$): −1.
- **Uniqueness condition wrong** for QP (missing full row rank): −1.
- **Little's formula misapplied** to unstable queue: −2.
- **Stationary distribution of M/M/1 wrong** ($\pi_n = \rho^n$ instead of $(1-\rho)\rho^n$): −2.
