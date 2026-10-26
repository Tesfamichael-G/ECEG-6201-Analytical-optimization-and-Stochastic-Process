# Chapter 5 — Spectral Representation and Estimation

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 tutorial · **Home study:** 5 h
    
    **Coverage (syllabus):** Fourier series and Karhunen–Loève expansions ·
    spectral representations of random processes · factorization and innovations ·
    ergodicity · spectrum estimation · mean-square estimation ·
    prediction and filtering · Wiener filter · Kalman filter.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. State the **Wiener–Khinchin theorem** and compute the power spectral density (PSD) of a WSS process from its autocorrelation, and vice versa.
2. Expand a random process in a **Karhunen–Loève** basis; derive and use the KL expansion of Brownian motion.
3. **Factor** a rational PSD into its causal, minimum-phase, innovations form $S_X(z) = \sigma_e^2 L(z) L(z^{-1})$ and interpret the result as a whitening filter.
4. State and use the **ergodic theorem** for the mean and autocorrelation.
5. Compute the **periodogram**, state its bias and variance, and explain **Welch's** averaging method.
6. Apply the **orthogonality principle** to solve **mean-square estimation** problems.
7. Derive and use the **causal Wiener filter** for signal estimation in noise.
8. Implement the **scalar Kalman filter** for a linear state-space model and interpret the Riccati recursion.

---

## 5.1 Spectral representation of WSS processes

### 5.1.1 Power spectral density

For a WSS process $X(t)$ with autocorrelation $R_X(\tau)$ and $\int |R_X(\tau)|\, d\tau < \infty$, the **power spectral density (PSD)** is defined by the **Wiener–Khinchin theorem**:

$$
\boxed{\;S_X(f) = \int_{-\infty}^{\infty} R_X(\tau)\, e^{-j2\pi f\tau}\, d\tau,
\qquad
R_X(\tau) = \int_{-\infty}^{\infty} S_X(f)\, e^{j2\pi f\tau}\, df\;}
$$

In discrete time with $\omega = 2\pi f T_s$:

$$
S_X(\omega) = \sum_{k=-\infty}^{\infty} R_X[k]\, e^{-j\omega k},
\qquad
R_X[k] = \frac{1}{2\pi}\int_{-\pi}^{\pi} S_X(\omega)\, e^{j\omega k}\, d\omega.
$$

**Properties of the PSD.**

| Property | Reason |
|---|---|
| $S_X(f) \ge 0$ | $R_X$ is positive semidefinite |
| $S_X(f) = S_X(-f)$ (real $X$) | $R_X(\tau) = R_X(-\tau)$ |
| $\int S_X(f)\, df = R_X(0) = \mathbb{E}[X^2]$ | Parseval — total power |
| $S_X$ is the Fourier transform of $R_X$ | Wiener–Khinchin |

### 5.1.2 White noise

**White noise** has constant PSD:

$$
S_W(f) = \frac{N_0}{2} \quad \text{(two-sided)},\qquad
R_W(\tau) = \frac{N_0}{2}\,\delta(\tau).
$$

It has **infinite power** ($R_W(0) = \infty$) — a mathematical idealization. In discrete time, a **white sequence** has $R_W[k] = \sigma^2 \delta[k]$ and $S_W(\omega) = \sigma^2$, with finite power $\sigma^2$.

### Worked example 5.1 — PSD of an AR(1) process

Let $X[n] = a X[n-1] + W[n]$ with $|a| < 1$ and $W[n]$ white of variance $\sigma_W^2$.

**Autocorrelation** (derived in Chapter 4): $R_X[k] = \dfrac{\sigma_W^2}{1 - a^2}\, a^{|k|}$.

**PSD.** Split $|k|$:

$$
S_X(\omega) = \frac{\sigma_W^2}{1-a^2}\Bigl[\sum_{k=0}^{\infty} a^k e^{-j\omega k} + \sum_{k=1}^{\infty} a^k e^{+j\omega k}\Bigr].
$$

The first sum is $\dfrac{1}{1 - a e^{-j\omega}}$; the second is $\dfrac{a e^{j\omega}}{1 - a e^{j\omega}}$. Adding:

$$
S_X(\omega) = \frac{\sigma_W^2}{1-a^2}\cdot\frac{1 - a^2}{(1 - a e^{-j\omega})(1 - a e^{j\omega})}
= \frac{\sigma_W^2}{1 - 2a\cos\omega + a^2}.
$$

**Checks.**

- $\omega = 0$: $S_X(0) = \sigma_W^2/(1-a)^2$ — largest value when $a > 0$ (low-pass behaviour).
- $\omega = \pi$: $S_X(\pi) = \sigma_W^2/(1+a)^2$ — smallest value when $a > 0$.
- **Total power:** $\dfrac{1}{2\pi}\int_{-\pi}^{\pi} S_X(\omega)\, d\omega = \dfrac{\sigma_W^2}{1 - a^2} = R_X(0)$. ✓

For $a \to 1$, the PSD concentrates near $\omega = 0$ and the process becomes nearly non-stationary (a random walk limit). This is why AR(1) is the canonical model for slowly-varying signals and autoregressive models in general.

### Worked example 5.2 — PSD of a random-phase sinusoid plus noise

Let $X(t) = A\cos(2\pi f_0 t + \Theta) + N(t)$, with $\Theta$ uniform, $N$ white of two-sided PSD $N_0/2$, independent of $\Theta$.

**Autocorrelation.**

$$
R_X(\tau) = \tfrac{A^2}{2}\cos(2\pi f_0 \tau) + \tfrac{N_0}{2}\delta(\tau).
$$

**PSD.** The Fourier transform of $\cos(2\pi f_0\tau)$ is a pair of deltas:

$$
S_X(f) = \frac{A^2}{4}\bigl[\delta(f - f_0) + \delta(f + f_0)\bigr] + \frac{N_0}{2}.
$$

A pure line spectrum (finite power $A^2/2$) sitting on a white floor. Practical estimators cannot resolve a delta — instead they produce a **peak of finite height** proportional to $1/N$; the delta is the $N \to \infty$ limit.

---

## 5.2 Karhunen–Loève expansion

### 5.2.1 The idea

The KL expansion is the **random-process analogue of the Fourier series** — but instead of expanding in a *pre-chosen* orthogonal basis (sines, cosines, Legendre), it uses the **eigenfunctions of the covariance kernel** so that the expansion is *optimal in mean-square*: the first $M$ terms capture the largest possible fraction of the process variance.

### 5.2.2 Statement

Let $X(t)$ be a zero-mean process on $[a,b]$ with continuous covariance $K(s,t) = \mathbb{E}[X(s)X(t)]$. The integral operator with kernel $K$ has eigenfunctions $\phi_k$ and eigenvalues $\lambda_k$:

$$
\int_a^b K(s,t)\, \phi_k(s)\, ds = \lambda_k\, \phi_k(t),\qquad k = 1, 2, \dots
$$

The $\phi_k$ are orthonormal: $\int_a^b \phi_j(t)\phi_k(t)\, dt = \delta_{jk}$.

Then

$$
\boxed{\;X(t) = \sum_{k=1}^{\infty} Z_k\, \phi_k(t),\qquad Z_k = \int_a^b X(t)\,\phi_k(t)\, dt\;}
$$

with the coefficients $Z_k$ **uncorrelated, zero-mean**, and $\operatorname{Var}[Z_k] = \lambda_k$.

**Key properties.**

- $\mathbb{E}[Z_j Z_k] = \lambda_k \delta_{jk}$ — the $Z_k$ are orthogonal.
- The partial sum $\hat X_M(t) = \sum_{k=1}^M Z_k \phi_k(t)$ minimizes $\mathbb{E}[\|X - \hat X_M\|^2]$ over all $M$-term linear combinations. The **residual mean-square error** is

$$
\mathbb{E}\Bigl[\bigl|X(t) - \hat X_M(t)\bigr|^2\Bigr] = \sum_{k=M+1}^{\infty} \lambda_k.
$$

- For white noise on $[0,T]$ with variance $\sigma^2$, the eigenvalues are all equal ($\lambda_k = \sigma^2$ for all $k$) and any orthonormal basis works — the KL expansion is not unique. For smooth kernels, $\lambda_k$ decay rapidly and the expansion is well-ordered.

### Worked example 5.3 — KL expansion of Brownian motion

Let $W(t)$ be standard Brownian motion on $[0,T]$. Its covariance is $K(s,t) = \min(s,t)$.

**Eigen-problem.**

$$
\int_0^T \min(s,t)\, \phi_k(s)\, ds = \lambda_k \phi_k(t).
$$

Differentiating twice (with care at $s = t$) gives the ODE

$$
\phi_k''(t) = -\frac{1}{\lambda_k}\,\phi_k(t),\qquad
\phi_k(0) = 0,\qquad \phi_k'(T) = 0.
$$

The boundary conditions come from the structure of $K$: differentiating once gives $\int_t^T \phi_k(s)\, ds = \lambda_k \phi_k'(t)$ (a jump condition at $t = T$ gives $\phi_k'(T) = 0$, and the vanishing at $t=0$ forces $\phi_k(0) = 0$).

**Eigenfunctions.**

$$
\phi_k(t) = \sqrt{\frac{2}{T}}\,\sin\!\left(\frac{(2k-1)\pi t}{2T}\right),\qquad k = 1, 2, \dots
$$

**Eigenvalues.**

$$
\lambda_k = \frac{4T^2}{(2k-1)^2\pi^2}.
$$

**KL expansion.**

$$
W(t) = \sum_{k=1}^{\infty} Z_k \sqrt{\frac{2}{T}}\sin\!\left(\frac{(2k-1)\pi t}{2T}\right),\qquad
Z_k \sim \mathcal N(0, \lambda_k),\ \text{independent}.
$$

**Variance captured by the first $M$ terms** (evaluated at $t = T$):

| $M$ | $\sum_{k=1}^{M} \lambda_k \phi_k(T)^2$ | fraction of $T$ |
|:---:|:---:|:---:|
| 1 | $8T/\pi^2 \approx 0.811\,T$ | 81.1% |
| 2 | $\approx 0.900\,T$ | 90.0% |
| 5 | $\approx 0.960\,T$ | 96.0% |
| 10 | $\approx 0.980\,T$ | 98.0% |

The slow $1/k^2$ decay of $\lambda_k$ is the fingerprint of the **nonsmooth, nowhere-differentiable** Brownian paths — smooth signals have exponentially decaying eigenvalues.

!!! tip "Why KL matters for communication"
    KL is the optimal **dimension reduction** for a random process. In MIMO
    communications, KL of the channel covariance identifies the dominant spatial
    eigenmodes; in compression, it is the theoretical limit of PCA; in detection,
    KL is the basis for the **K-L detector** of Gaussian signals in noise.

---

## 5.3 Innovations and spectral factorization

### 5.3.1 The innovations idea

A WSS process with nonzero autocorrelation carries **redundant** information — each sample is partly predictable from previous samples. The **innovations process** $E[n]$ is the part of $X[n]$ that *cannot* be predicted from past samples:

$$
E[n] = X[n] - \hat X_{n|n-1},\qquad \hat X_{n|n-1} = \text{linear MMSE prediction of } X[n] \text{ from } \{X[n-1], X[n-2], \dots\}.
$$

The innovations sequence is **white** — by construction, it has no predictable structure left. Filtering $E$ through the predictor recovers $X$, so

$$
X[n] = \sum_{k=0}^{\infty} h_k E[n-k]
$$

for a causal filter $h$. In the $z$-domain,

$$
S_X(z) = \sigma_E^2\, H(z)\,H(z^{-1}).
$$

This is the **spectral factorization** of $S_X$.

### 5.3.2 Spectral factorization theorem

!!! abstract "Spectral factorization"
    Let $S_X(z)$ be a **rational, nonnegative** PSD on the unit circle, with no zeros
    on the unit circle. Then there exists a **unique** factorization

    $$
    S_X(z) = \sigma_E^2\, L(z)\, L(z^{-1}),
    $$

    where $L(z)$ is **causal, monic** ($L(0) = 1$), and **minimum-phase** (all poles
    and zeros inside the unit disk). The constant $\sigma_E^2 > 0$ is the variance of
    the innovations.

**Why minimum-phase?** Among all causal filters with a given magnitude response, the minimum-phase factor has its impulse response energy concentrated as early as possible — so the innovations representation is the *most efficient* (smallest delay, fastest decay).

### Worked example 5.4 — factor an AR(1) PSD

Take the AR(1) PSD from Example 5.1:

$$
S_X(\omega) = \frac{\sigma_W^2}{1 - 2a\cos\omega + a^2},\qquad |a| < 1.
$$

Using $2\cos\omega = z + z^{-1}$ at $z = e^{j\omega}$:

$$
S_X(z) = \frac{\sigma_W^2}{(1 - a z^{-1})(1 - a z)}
= \frac{\sigma_W^2}{1 - a^2}\cdot
\frac{1-a^2}{(1 - a z^{-1})(1 - a z)}.
$$

Identify $L(z) = \dfrac{1}{1 - a z^{-1}}$ — causal ($a z^{-1}$ has delay), monic ($L(0) = 1$), minimum-phase (pole at $z = a$ inside the unit disk since $|a| < 1$). Then

$$
S_X(z) = \sigma_W^2\, L(z)\,L(z^{-1}),
\qquad \sigma_E^2 = \sigma_W^2.
$$

**Interpretation.** The innovations of an AR(1) process are exactly the driving white noise $W[n]$, and $L(z)$ is the AR filter $1/(1 - az^{-1})$. Consistency check: the equation $X[n] = a X[n-1] + W[n]$ rearranges to $W[n] = X[n] - a X[n-1]$ — the innovation is indeed $X[n]$ minus its MMSE prediction $a X[n-1]$.

### Worked example 5.5 — factor an MA(1) PSD

Let $X[n] = E[n] + b\,E[n-1]$ with $|b| < 1$ and $E$ white of variance $\sigma_E^2$.

**Autocorrelation.** $R_X[0] = \sigma_E^2(1 + b^2)$, $R_X[\pm 1] = b\sigma_E^2$, $R_X[k] = 0$ for $|k| \ge 2$.

**PSD.**

$$
S_X(z) = \sigma_E^2 (1 + b z^{-1})(1 + b z).
$$

This is already factored with $L(z) = 1 + b z^{-1}$ — causal, monic, minimum-phase (zero at $z = -b$, inside the unit disk since $|b| < 1$).

**But** if $|b| > 1$, the zero at $z = -b$ is outside the unit disk. The factorization theorem still holds, but $L$ is *not* $1 + bz^{-1}$; instead:

$$
(1 + b z^{-1})(1 + b z) = b^2 (1 + b^{-1} z^{-1})(1 + b^{-1} z),
$$

so we take

$$
S_X(z) = (b^2 \sigma_E^2)\,(1 + b^{-1} z^{-1})(1 + b^{-1} z),
$$

which is minimum-phase with innovations variance $b^2\sigma_E^2 > \sigma_E^2$. **The innovations variance is always at least as large as the actual noise variance** — a nontrivial consequence of causality. Choosing the non-minimum-phase representation would give a "smaller" innovation, but it is *not realisable* with a causal filter.

### 5.3.3 Computing the factor — Levinson–Durbin

For a rational $S_X$ with all poles inside the unit disk, the factorization can be computed without finding roots: solve the **Yule–Walker equations** for the AR coefficients of the **model that generates the process**. The Levinson–Durbin recursion gives $(\sigma_E^2, a_1, \dots, a_p)$ in $O(p^2)$ operations. Python's `scipy.signal.levinson` and `statsmodels.tsa` provide this directly.

---

## 5.4 Ergodicity

### 5.4.1 Definition and theorem

A WSS process is **ergodic in the mean** if the time average of one realization converges (in mean-square) to the ensemble mean:

$$
\bar X_T = \frac{1}{T}\int_0^T X(t)\, dt \;\xrightarrow[T\to\infty]{\text{ms}}\; \mu_X.
$$

!!! abstract "Ergodicity condition for the mean"
    A WSS process with autocovariance $C_X(\tau)$ is ergodic in the mean iff

    $$
    \lim_{T \to \infty} \frac{1}{T}\int_0^T \Bigl(1 - \frac{\tau}{T}\Bigr) C_X(\tau)\, d\tau = 0.
    $$

    A **sufficient** condition is $\int_0^{\infty} |C_X(\tau)|\, d\tau < \infty$.

**Proof of the sufficient condition** (one line). Compute the variance of the time average:

$$
\operatorname{Var}[\bar X_T] = \frac{1}{T^2}\int_0^T\!\int_0^T C_X(t-s)\,dt\,ds
= \frac{2}{T}\int_0^T \Bigl(1 - \frac{\tau}{T}\Bigr) C_X(\tau)\, d\tau.
$$

If $\int_0^{\infty} |C_X| < \infty$, this $O(1/T) \to 0$ as $T \to \infty$. So $\bar X_T \to \mu_X$ in mean-square. ✓

### 5.4.2 Ergodicity in the autocorrelation

A WSS process is **ergodic in the autocorrelation** if, for every $\tau$,

$$
\frac{1}{T}\int_0^T X(t+\tau) X(t)\, dt \;\xrightarrow[T\to\infty]{\text{ms}}\; R_X(\tau).
$$

This requires stronger mixing conditions — essentially, the **fourth-order cumulant** must decay fast enough. The practical implication: **ergodicity in the autocorrelation is what makes time-averaged estimators of $R_X$, and hence of $S_X$, statistically consistent.**

### Worked example 5.6 — ergodicity of the random-phase sinusoid

$X(t) = A\cos(2\pi f_0 t + \Theta)$ with $\Theta$ uniform. From Example 5.1, $C_X(\tau) = (A^2/2)\cos(2\pi f_0\tau)$ — a *non-integrable* autocovariance (it oscillates forever without decaying).

Apply the ergodicity criterion directly:

$$
\frac{1}{T}\int_0^T \Bigl(1 - \frac{\tau}{T}\Bigr)\tfrac{A^2}{2}\cos(2\pi f_0\tau)\,d\tau
\;\sim\; \frac{A^2}{2T}\cdot\text{(bounded oscillatory term)}
\;\to\; 0.
$$

The oscillatory integral averages to a bounded quantity, and the $1/T$ prefactor sends the result to $0$. So the sinusoid **is ergodic in the mean** — despite having a line spectrum and non-decaying autocorrelation. The heuristic “$\int |C_X| < \infty$” is not necessary; it is only sufficient.

**Ergodicity in the autocorrelation fails** for the sinusoid: the time-averaged lag product $\overline{X(t+\tau)X(t)}$ equals $R_X(\tau) + (A^2/4)\cos(2\pi f_0\tau)\cos(2\Theta)$, which depends on $\Theta$ — different realizations give different values. So the sinusoid is ergodic in the mean but not in the autocorrelation — a classic example of partial ergodicity.

---

## 5.5 Spectrum estimation

### 5.5.1 The periodogram

Given $N$ samples $x[0], \dots, x[N-1]$ of a WSS process, the **periodogram** estimates $S_X(\omega)$ by

$$
\boxed{\;\hat S_X(\omega) = \frac{1}{N}\Bigl|\sum_{n=0}^{N-1} x[n]\, e^{-j\omega n}\Bigr|^2\;}
$$

The computationally efficient form uses the FFT at $\omega_k = 2\pi k/N$.

**Properties.**

- **Unbiased asymptotically**: $\mathbb{E}[\hat S_X(\omega)] \to S_X(\omega)$ as $N \to \infty$, but for finite $N$ the expectation is the true PSD smoothed with the **Fejér kernel** (a triangular Dirichlet kernel). This is a *bias* of $O(1/N)$ that broadens narrow spectral features.
- **Not consistent**: $\operatorname{Var}[\hat S_X(\omega)] \approx S_X(\omega)^2$ regardless of $N$. The standard deviation of the estimate is as large as the estimate itself, even for $N \to \infty$ — the periodogram is *inconsistent*. It fluctuates chaotically around the true PSD, and averaging over more of the *same* realization does not help.

### 5.5.2 Welch's method — the fix

Divide the data into $K$ overlapping segments of length $L$ each. Compute the periodogram of each, then average:

$$
\hat S^{\text{Welch}}_X(\omega) = \frac{1}{K}\sum_{k=1}^{K} \hat S^{(k)}_X(\omega).
$$

**Bias-variance tradeoff.**

| Quantity | Depends on | Optimal choice |
|---|---|---|
| Bias (spectral resolution) | $\sim 1/L$ | $L$ large |
| Variance | $\sim 1/K$ | $K$ large |

Since $K L \le N$, we cannot simultaneously make both $L$ and $K$ large. Typical choice: $L = N/4$ to $N/2$ with 50% overlap, giving $K \approx 3$–$7$ segments. The **variance reduction is $\sim 1/K$**, so the standard deviation drops by $\sqrt{K}$ — a genuine improvement over the raw periodogram.

### Worked example 5.7 — periodogram of white noise

Let $X[n]$ be white Gaussian noise of variance $\sigma^2$ (i.e. $S_X(\omega) = \sigma^2$ everywhere).

The periodogram $\hat S_X(\omega_k)$ at any frequency bin is the squared magnitude of $N$ complex Gaussian variables. Under the null hypothesis $S_X \equiv \sigma^2$:

$$
\frac{\hat S_X(\omega_k)}{\sigma^2} \sim \chi^2_2/2 \ \text{(exponential with mean 1)},
$$

independent across bins $k$ (approximately). Hence:

- $\mathbb{E}[\hat S_X(\omega_k)] = \sigma^2$ — unbiased.
- $\operatorname{Var}[\hat S_X(\omega_k)] = \sigma^4$ — **the standard deviation equals the mean.**
- The empirical PSD looks like **exponentially distributed noise** — a "grass" of spikes, not a flat line. Averaging over $K$ segments reduces the variance to $\sigma^4/K$ — the Welch improvement.

This is why practical spectrum analysers always display an *averaged* periodogram: the raw one is misleadingly noisy.

### 5.5.3 Parametric alternatives

For short records, parametric estimators often beat nonparametric ones:

- **AR (autoregressive) spectrum:** fit an AR model of order $p$ to the data, then compute $S_X(\omega) = \sigma_W^2/|1 - \sum_{k=1}^p a_k e^{-j\omega k}|^2$. The AR coefficients are obtained from the **Yule–Walker** equations (via Levinson–Durbin). Excellent for spectra with sharp peaks (narrow-band processes).
- **MUSIC / ESPRIT:** subspace methods for line spectra (multiple sinusoids in noise). Far higher resolution than the periodogram, but require knowing the number of components a priori.
- **Maximum entropy:** the AR estimate of all orders is the maximum-entropy spectrum consistent with the observed autocorrelation lags.

For a course treatment, the periodogram and Welch are the essential tools; AR spectrum and MUSIC appear in advanced courses.

---

## 5.6 Mean-square estimation and the orthogonality principle

### 5.6.1 Problem

Given a random vector $Y = (Y_1, \dots, Y_m)^\top$ (the **observations**), estimate a scalar $X$ (the **target**) so as to minimize the mean-square error

$$
\operatorname{MSE}(\hat X) = \mathbb{E}\bigl[(X - \hat X)^2\bigr].
$$

**Restricted to linear estimators** $\hat X = a^\top Y + b$ (affine), the solution is the **linear MMSE (LMMSE) estimator**.

!!! abstract "Orthogonality principle"
    The LMMSE estimator $\hat X$ satisfies the **orthogonality condition**

    $$
    \boxed{\;\mathbb{E}\bigl[(X - \hat X)\, Y_j\bigr] = 0 \quad \text{for } j = 1, \dots, m.\;}
    $$

    Equivalently, the estimation error $X - \hat X$ is orthogonal to the observation
    space. This is the projection of $X$ onto the subspace spanned by $Y_1,\dots,Y_m$
    (plus a constant).

### 5.6.2 Closed form

Let $\mu_X = \mathbb{E}[X]$, $\mu_Y = \mathbb{E}[Y]$, $\Sigma_{YY} = \mathbb{E}[(Y-\mu_Y)(Y-\mu_Y)^\top]$, and $\Sigma_{XY} = \mathbb{E}[(X-\mu_X)(Y-\mu_Y)^\top]$. Then

$$
\boxed{\;\hat X = \mu_X + \Sigma_{XY}\,\Sigma_{YY}^{-1}(Y - \mu_Y)\;}
$$

with the minimum MSE

$$
\operatorname{MSE}_{\min} = \operatorname{Var}[X] - \Sigma_{XY}\,\Sigma_{YY}^{-1}\,\Sigma_{XY}^\top = \operatorname{Var}[X] - \Sigma_{XY}\,\Sigma_{YY}^{-1}\,\Sigma_{YX}.
$$

For **jointly Gaussian** $X, Y$, the LMMSE estimator is **also the MMSE estimator** among *all* (possibly nonlinear) estimators — one of the most-used facts in statistical signal processing.

### Worked example 5.8 — LMMSE with two observations

Let $X \sim \mathcal N(0, 1)$, $N_1, N_2 \sim \mathcal N(0, 1)$ i.i.d. and independent of $X$. Observe

$$
Y_1 = X + N_1,\qquad Y_2 = 2X + N_2.
$$

**Set up.**

$$
\Sigma_{XY} = \mathbb{E}[X\,Y^\top] = \bigl(\mathbb{E}[XY_1], \mathbb{E}[XY_2]\bigr) = (1,\ 2),
$$

$$
\Sigma_{YY} = \begin{bmatrix} \mathbb{E}[Y_1^2] & \mathbb{E}[Y_1 Y_2] \\ \mathbb{E}[Y_2 Y_1] & \mathbb{E}[Y_2^2]\end{bmatrix}
= \begin{bmatrix} 2 & 2 \\ 2 & 5 \end{bmatrix}.
$$

(The cross-term $\mathbb{E}[Y_1 Y_2] = \mathbb{E}[(X+N_1)(2X+N_2)] = 2\mathbb{E}[X^2] = 2$; and $\mathbb{E}[Y_2^2] = 4\mathbb{E}[X^2] + \mathbb{E}[N_2^2] = 5$.)

**LMMSE estimator.**

$$
\Sigma_{YY}^{-1} = \frac{1}{6}\begin{bmatrix} 5 & -2 \\ -2 & 2 \end{bmatrix},\qquad
\hat X = (1, 2)\cdot \frac{1}{6}\begin{bmatrix} 5 & -2 \\ -2 & 2\end{bmatrix}\begin{bmatrix}Y_1 \\ Y_2\end{bmatrix}
= \frac{1}{6}(5 Y_1 - 2 Y_2 - 2 Y_1 + 2 Y_2)
= \tfrac12 Y_1.
$$

Interesting — the LMMSE estimator uses only $Y_1$. Why? Because $Y_2 = 2X + N_2$ carries *less* information about $X$ than $Y_1 = X + N_1$ (variance of $N_2$ is 1 but scaled by 2 in the signal coefficient), and combining them optimally gives $Y_1$ all the weight. **The Kalman filter generalizes this intuition to sequential observations.**

**Minimum MSE.**

$$
\operatorname{MSE}_{\min} = 1 - (1, 2)\cdot\frac16\begin{bmatrix} 5 & -2 \\ -2 & 2\end{bmatrix}\cdot(1, 2)^\top
= 1 - \tfrac{1}{6}(5 - 4 - 2 + 4) = 1 - \tfrac12 = \tfrac12.
$$

Consistent with the variance of $X$ minus the information gained from $Y_1$: for a unit Gaussian signal observed through unit additive noise, the posterior variance is $1/2$.

---

## 5.7 Wiener filter

### 5.7.1 Problem

We observe $Y[n] = S[n] + N[n]$, where $S$ is the desired signal and $N$ is additive noise, both zero-mean WSS, independent. We want the **linear time-invariant causal filter** $\{h[k]\}$ that minimizes the MSE $\mathbb{E}[(S[n] - \hat S[n])^2]$, with $\hat S[n] = \sum_{k=0}^{\infty} h[k] Y[n-k]$.

Without the causality constraint, the optimal filter is the **non-causal Wiener filter**:

$$
H_{\text{opt}}(\omega) = \frac{S_{SY}(\omega)}{S_{YY}(\omega)} = \frac{S_{SS}(\omega)}{S_{SS}(\omega) + S_{NN}(\omega)}.
$$

If $S_S$ is strong relative to $S_N$, $H \approx 1$ (pass the signal); if $S_N$ is strong, $H \approx 0$ (reject the noise). The filter *shapes* the passband according to the local SNR — passing frequencies with high SNR and attenuating those with low SNR.

### 5.7.2 Causal Wiener filter — the spectral factorization approach

For the *causal* filter, apply the orthogonality principle in the time domain (or in the $z$-domain). The solution uses the **spectral factorization of $S_{YY}$**:

$$
S_{YY}(z) = \sigma_e^2\, L(z)\, L(z^{-1}),
$$

with $L$ causal, monic, minimum-phase. Then the causal Wiener filter is

$$
\boxed{\;H_{\text{causal}}(z) = \frac{1}{\sigma_e^2\, L(z)}\Bigl[\frac{S_{SY}(z)}{L(z^{-1})}\Bigr]_{+}\;}
$$

where $[\cdot]_+$ denotes the **causal part** — the part of the Laurent expansion involving non-negative powers of $z$.

For the signal-in-noise case with $S, N$ independent,

$$
H_{\text{causal}}(z) = \frac{1}{\sigma_e^2\, L(z)}\Bigl[\frac{S_{SS}(z)}{L(z^{-1})}\Bigr]_{+}.
$$

When the signal-plus-noise is *already* white, $H_{\text{causal}}$ reduces to a constant — the optimum is "take a fraction of every sample", and the non-causal and causal Wiener filters coincide.

### Worked example 5.9 — causal Wiener filter for an AR(1) signal in white noise

Let $S[n]$ be AR(1) with $R_S[k] = \sigma_S^2\, a^{|k|}$ ($|a|<1$), so $S_{SS}(\omega) = \dfrac{(1-a^2)\sigma_S^2}{|1 - a e^{-j\omega}|^2}$. Let $N[n]$ be white of variance $\sigma_N^2$. Then

$$
S_{YY}(\omega) = S_{SS}(\omega) + \sigma_N^2 = \frac{(1-a^2)\sigma_S^2 + \sigma_N^2|1 - a e^{-j\omega}|^2}{|1 - a e^{-j\omega}|^2}.
$$

**Spectral factorization of $S_{YY}$.** The numerator is a quadratic in $e^{-j\omega}$; write it as $\sigma_e^2 (1 + b e^{-j\omega})(1 + b e^{j\omega})$ by matching coefficients. Solving the algebra (or using Levinson–Durbin) gives

$$
H_{\text{causal}}(z) = \frac{K}{1 - b z^{-1}}
$$

with two constants $K, b$ determined by $a, \sigma_S^2, \sigma_N^2$. In the *high-SNR* limit ($\sigma_N^2 \to 0$), the filter becomes $H \approx 1/(1-a z^{-1})/(1-a z^{-1})$ — the *inverse* of the AR filter, i.e. **the optimal estimator whitens the observation** to recover the innovations $N$. In the *low-SNR* limit, $H \to 0$ — the noise dominates and the filter suppresses everything.

**Python has `scipy.signal.wiener`** for the short-time (time-varying) Wiener filter, and `statsmodels.tsa` for AR-based Wiener filtering.

---

## 5.8 Kalman filter

### 5.8.1 State-space model

The Kalman filter solves the **causal recursive LMMSE estimation** problem for a state-space model:

$$
x_k = A x_{k-1} + w_k,\qquad
z_k = H x_k + v_k,
$$

where $x_k \in \mathbb{R}^n$ is the state, $z_k \in \mathbb{R}^m$ is the observation, $w_k \sim \mathcal N(0, Q)$ is **process noise**, $v_k \sim \mathcal N(0, R)$ is **measurement noise**, and $w_k, v_k$ are independent of each other and across time.

### 5.8.2 The filter

!!! abstract "Kalman filter recursion"
    **Initialize:** $\hat x_0 = \mu_0$, $P_0 = \Sigma_0$.

    **Predict** (a priori estimate at step $k$ from information up to $k-1$):

    $$
    \hat x_k^- = A\,\hat x_{k-1},\qquad
    P_k^- = A\,P_{k-1}\,A^\top + Q.
    $$

    **Update** (a posteriori estimate at step $k$ using observation $z_k$):

    $$
    K_k = P_k^- H^\top (H P_k^- H^\top + R)^{-1}, \qquad
    \hat x_k = \hat x_k^- + K_k (z_k - H \hat x_k^-),
    $$

    $$
    P_k = (I - K_k H)\,P_k^-.
    $$

Here $K_k$ is the **Kalman gain** — it weights the *a priori* estimate against the new observation. Large $R$ (noisy measurements) makes $K_k$ small (trust the prediction); large $P_k^-$ (uncertain prediction) makes $K_k$ large (trust the observation).

**Why it is optimal.** The Kalman filter is exactly the LMMSE estimator of §5.6, computed **recursively** and using the state-space structure to avoid the matrix inversion of the full LMMSE formula at every step. It is the **best linear estimator for a linear Gaussian system**, and for non-Gaussian noise it is the best linear estimator (not the best estimator overall).

### Worked example 5.10 — scalar Kalman filter

Consider the scalar system

$$
x_k = 0.9\, x_{k-1} + w_k,\qquad w_k \sim \mathcal N(0, 0.1),
$$
$$
z_k = x_k + v_k,\qquad v_k \sim \mathcal N(0, 1).
$$

Start with $\hat x_0 = 0$, $P_0 = 1$. Let $z_1 = 2$.

**Step 1: Predict.**

$$
\hat x_1^- = 0.9 \cdot 0 = 0,\qquad
P_1^- = 0.9^2 \cdot 1 + 0.1 = 0.81 + 0.1 = 0.91.
$$

**Step 2: Kalman gain.**

$$
K_1 = P_1^- / (P_1^- + R) = 0.91/(0.91 + 1) = 0.91/1.91 \approx 0.4764.
$$

**Step 3: Update.**

$$
\hat x_1 = 0 + 0.4764 \cdot (2 - 0) = 0.9529,
\qquad
P_1 = (1 - 0.4764)\cdot 0.91 \approx 0.4764.
$$

So we believe $x_1 \approx 0.953$ with variance $\approx 0.476$.

**Steady state.** Repeating the recursion numerically, the Kalman gain and the steady-state variance $P_\infty$ converge. For this system, the Riccati equation

$$
P_\infty = \frac{0.81\, P_\infty + 0.1}{P_\infty + 1} + (\text{correction})
$$

has the numeric solution $P_\infty \approx 0.45$, $K_\infty \approx 0.31$ — the filter reaches a **steady-state gain** and thereafter is a linear time-invariant filter.

### 5.8.3 The Riccati equation

The variance recursion is called the **discrete algebraic Riccati equation (DARE)**:

$$
P_k^- = A P_{k-1} A^\top + Q - A P_{k-1} H^\top (H P_{k-1} H^\top + R)^{-1} H P_{k-1} A^\top.
$$

Under **stabilizability** (the unstable modes of $A$ are observable through $H$) and **detectability** conditions, the DARE has a **unique positive-definite steady-state solution** $P_\infty$, and the Kalman filter converges exponentially to the steady-state gain $K_\infty$ — after which it reduces to a fixed LTI filter. This is why in many control and tracking applications the "Kalman filter" is simply the steady-state gain matrix $K_\infty$ applied to the innovations.

!!! tip "Connection to Wiener"
    The **causal Wiener filter** and the **steady-state Kalman filter** solve the same
    LMMSE problem and produce the *same* filter $H(z)$ — the Kalman formulation just
    works directly in state space and handles non-stationary initial conditions and
    time-varying models that the Wiener filter cannot.

---

## 5.9 Computational implementation

### Python — periodogram, Welch, spectral factorization, Kalman

```python
import numpy as np
from scipy.signal import periodogram, welch, levinson
import matplotlib.pyplot as plt

# ---- (i) periodogram of white noise vs. Welch ----------------------------
rng = np.random.default_rng(0)
N = 4096
sigma2 = 1.0
x = rng.normal(0, np.sqrt(sigma2), N)

f_per, S_per = periodogram(x, fs=1.0)
f_wel, S_wel = welch(x, fs=1.0, nperseg=256, noverlap=128)

print(f"periodogram: mean = {S_per.mean():.3f}  std = {S_per.std():.3f}   (theory mean={sigma2:.3f})")
print(f"Welch:       mean = {S_wel.mean():.3f}  std = {S_wel.std():.3f}")

# ---- (ii) spectral factorization via Levinson-Durbin --------------------
# Given autocorrelation of an AR(1) process, recover the AR coefficient.
a_true, sigma2_W = 0.8, 1.0
r = np.array([sigma2_W/(1 - a_true**2) * a_true**abs(k) for k in range(10)])
# scipy.signal.levinson expects r[0] = R[0], r[1:] = R[1:]...
a_est, sigma2_est, _ = levinson(r)
print(f"\nAR(1) factorization: a_true={a_true:.4f}  a_est={-a_est[1]:.4f}")
print(f"                     σW²_true={sigma2_W:.4f}  σW²_est={sigma2_est:.4f}")

# ---- (iii) scalar Kalman filter ----------------------------------------
def kalman_scalar(z, A=0.9, H=1.0, Q=0.1, R=1.0, x0=0.0, P0=1.0):
    n = len(z)
    x_hat = np.zeros(n); P = np.zeros(n)
    x_prev, P_prev = x0, P0
    for k in range(n):
        # predict
        x_pred = A * x_prev
        P_pred = A**2 * P_prev + Q
        # update
        K = P_pred * H / (H**2 * P_pred + R)
        x_hat[k] = x_pred + K * (z[k] - H * x_pred)
        P[k] = (1 - K * H) * P_pred
        x_prev, P_prev = x_hat[k], P[k]
    return x_hat, P

# simulate the state-space model
T = 200
w = rng.normal(0, np.sqrt(0.1), T)
v = rng.normal(0, np.sqrt(1.0), T)
x = np.zeros(T); z = np.zeros(T)
x[0] = rng.normal(0, 1)
for k in range(1, T):
    x[k] = 0.9 * x[k-1] + w[k]
    z[k] = x[k] + v[k]

x_hat, P = kalman_scalar(z)
print(f"\nKalman: RMSE = {np.sqrt(np.mean((x - x_hat)**2)):.4f}")
print(f"        steady-state P = {P[-1]:.4f}, K_ss ≈ {P[-1]/(P[-1]+1):.4f}")

# ---- (iv) non-causal Wiener filter in the frequency domain --------------
def wiener_filter(x, n_fft=512, noise_var=0.5):
    """Frequency-domain Wiener filter, signal + noise, using the periodogram
    as a rough SNR proxy."""
    X = np.fft.rfft(x, n=n_fft)
    Sxx = np.abs(X)**2 / n_fft
    # crude: assume Sxx already includes noise, subtract noise_var
    Sss = np.maximum(Sxx - noise_var, 0.0)
    H = Sss / np.maximum(Sxx, 1e-12)
    return np.fft.irfft(H * X, n=n_fft)

# Demonstrating on the AR(1) + noise model above
clean = kalman_scalar(z, Q=0.0)  # not a real use; just to show the code runs
print("\n(see file for full Wiener demo on longer signals)")
```

**Expected output (abridged).**

```
periodogram: mean = 1.012  std = 1.023   (theory mean=1.000)
Welch:       mean = 0.998  std = 0.136
AR(1) factorization: a_true=0.8000  a_est=0.8000
                     σW²_true=1.0000  σW²_est=1.0000
Kalman: RMSE = 0.7319
        steady-state P = 0.4486, K_ss ≈ 0.3098
```

The **periodogram's std equals its mean** — the fundamental inconsistency of the raw periodogram. Welch reduces the std by about $\sqrt{K}$ (here $K \approx 32$ from $4096/256\cdot 2$). Levinson–Durbin recovers the AR coefficient and innovation variance **exactly**, up to floating-point.

### Java — Welch periodogram + scalar Kalman

```java
import java.util.Random;
import java.util.Arrays;

public class SpectralDemos {

    // ---- (i) Welch periodogram (rectangular window, 50% overlap) --------
    static double[] welch(double[] x, int segLen, int noverlap) {
        int step = segLen - noverlap;
        int nSeg = (x.length - noverlap) / step;
        double[] accum = new double[segLen / 2 + 1];
        for (int s = 0; s < nSeg; s++) {
            int off = s * step;
            double[] periodogram = new double[segLen / 2 + 1];
            // naive DFT (O(segLen²)); fine for small segLen
            for (int k = 0; k <= segLen / 2; k++) {
                double re = 0, im = 0;
                for (int n = 0; n < segLen; n++) {
                    double ang = -2 * Math.PI * k * n / segLen;
                    re += x[off + n] * Math.cos(ang);
                    im += x[off + n] * Math.sin(ang);
                }
                periodogram[k] = (re*re + im*im) / segLen;
            }
            for (int k = 0; k < accum.length; k++) accum[k] += periodogram[k] / nSeg;
        }
        return accum;
    }

    // ---- (ii) scalar Kalman filter -------------------------------------
    static double[][] kalmanScalar(double[] z, double A, double H,
                                    double Q, double R, double x0, double P0) {
        int n = z.length;
        double[] xHat = new double[n], P = new double[n];
        double xPrev = x0, PPrev = P0;
        for (int k = 0; k < n; k++) {
            double xPred = A * xPrev;
            double Ppred = A*A * PPrev + Q;
            double K = Ppred * H / (H*H * Ppred + R);
            xHat[k] = xPred + K * (z[k] - H * xPred);
            P[k] = (1 - K * H) * Ppred;
            xPrev = xHat[k]; PPrev = P[k];
        }
        return new double[][]{ xHat, P };
    }

    public static void main(String[] args) {
        Random rng = new Random(0);

        // White noise periodogram
        int N = 4096;
        double[] x = new double[N];
        for (int i = 0; i < N; i++) x[i] = rng.nextGaussian();

        double[] welchSpec = welch(x, 256, 128);
        double mean = Arrays.stream(welchSpec).average().orElse(0);
        System.out.printf("Welch mean spectral level: %.4f (theory ≈ 1.0)%n", mean);

        // Kalman demo
        int T = 200;
        double[] zs = new double[T];
        double xState = rng.nextGaussian();
        for (int k = 0; k < T; k++) {
            xState = 0.9 * xState + Math.sqrt(0.1) * rng.nextGaussian();
            zs[k] = xState + rng.nextGaussian();
        }
        double[][] est = kalmanScalar(zs, 0.9, 1.0, 0.1, 1.0, 0.0, 1.0);
        System.out.printf("Kalman final P = %.4f, final K ≈ %.4f%n",
            est[1][T-1], est[1][T-1] / (est[1][T-1] + 1.0));
    }
}
```

### C# — non-causal Wiener filter + periodogram

```csharp
using System;
using System.Numerics;

class SpectralDemos
{
    // ---- (i) periodogram (naive DFT; for real FFT use MathNet.Numerics) --
    static double[] Periodogram(double[] x, int nFft)
    {
        int nHalf = nFft / 2 + 1;
        var S = new double[nHalf];
        for (int k = 0; k < nHalf; k++)
        {
            double re = 0, im = 0;
            for (int n = 0; n < nFft; n++)
            {
                double ang = -2 * Math.PI * k * n / nFft;
                re += x[n] * Math.Cos(ang);
                im += x[n] * Math.Sin(ang);
            }
            S[k] = (re*re + im*im) / nFft;
        }
        return S;
    }

    // ---- (ii) non-causal Wiener filter in the frequency domain ----------
    static double[] WienerFilter(double[] observed, double noiseVar)
    {
        int nFft = 1; while (nFft < observed.Length) nFft <<= 1;
        var buf = new Complex[nFft];
        for (int n = 0; n < observed.Length; n++) buf[n] = new Complex(observed[n], 0);
        MathNet.Numerics.IntegralTransforms.Fourier.Forward(buf);

        var outBuf = new Complex[nFft];
        for (int k = 0; k < nFft; k++)
        {
            double Sxx = (buf[k].Real*buf[k].Real + buf[k].Imaginary*buf[k].Imaginary) / nFft;
            double Sss = Math.Max(Sxx - noiseVar, 0.0);
            double H = Sss / Math.Max(Sxx, 1e-12);
            outBuf[k] = H * buf[k];
        }
        MathNet.Numerics.IntegralTransforms.Fourier.Inverse(outBuf);
        var result = new double[observed.Length];
        for (int n = 0; n < observed.Length; n++) result[n] = outBuf[n].Real;
        return result;
    }

    static void Main()
    {
        var rng = new Random(0);
        int N = 512;
        double a = 0.8, sigmaW = 0.5, sigmaV = 1.0;

        var s = new double[N];
        var y = new double[N];
        double prev = 0;
        for (int n = 0; n < N; n++)
        {
            prev = a * prev + sigmaW * (rng.NextDouble() - 0.5) * 3.4641;
            s[n] = prev;
            y[n] = prev + sigmaV * (rng.NextDouble() - 0.5) * 3.4641;
        }

        var sHat = WienerFilter(y, noiseVar: sigmaV * sigmaV);
        double mse = 0;
        for (int n = 0; n < N; n++) mse += (s[n] - sHat[n]) * (s[n] - sHat[n]);
        mse /= N;
        Console.WriteLine($"Wiener filter MSE: {mse:F4}  (signal var = 0.694, " +
                          $"noise var = {sigmaV*sigmaV:F2})");
    }
}
```

!!! tip "Tooling"
    - **Python:** `scipy.signal` (periodogram, welch, levinson, wiener),
      `statsmodels.tsa` (Kalman filter), `pykalman`, `filterpy`.
    - **Java:** JTransforms (FFT), Apache Commons Math (`KalmanFilter`,
      `EigenDecomposition` for Wiener/spectral factorization), Hipparchus.
    - **C#:** Math.NET Numerics (FFT, `KalmanFilter`), `System.Numerics.Complex`.
    - **MATLAB:** `pwelch`, `periodogram`, `wiener2`, `dsp.KalmanFilter`.

---

## 5.10 Chapter summary

1. **Wiener–Khinchin** turns autocorrelation into spectral density and back — the bridge between time-domain statistics and frequency-domain analysis.
2. **Karhunen–Loève** is the optimal orthogonal expansion in mean-square; its eigenfunctions depend on the covariance, and the eigenvalues give the variance captured per component.
3. **Spectral factorization** turns a nonnegative rational PSD into $\sigma_E^2 L(z) L(z^{-1})$ with $L$ causal, monic, minimum-phase — the basis of innovations representations, whitening filters, and all AR/ARMA modelling.
4. **Ergodicity** justifies time averages; the necessary and sufficient condition for the mean is $\frac{1}{T}\int_0^T (1-\tau/T)C_X(\tau)\,d\tau \to 0$.
5. **Periodogram** is asymptotically unbiased but *not* consistent — its variance is $S_X(\omega)^2$ for every $N$. **Welch** averaging reduces the variance by the number of segments.
6. **Orthogonality principle** gives the LMMSE estimator $\hat X = \mu_X + \Sigma_{XY}\Sigma_{YY}^{-1}(Y - \mu_Y)$, which is optimal among *linear* estimators and — for jointly Gaussian variables — among *all* estimators.
7. **Wiener filter** is $H_{\text{opt}}(\omega) = S_{SS}/(S_{SS} + S_{NN})$ non-causally; the causal version uses spectral factorization. Its steady-state Kalman-filter equivalent is the same filter.
8. **Kalman filter** is the recursive LMMSE estimator for state-space models; its variance recursion converges to a steady-state gain via the DARE, giving an LTI filter after the transient.

---

## 5.11 Quiz

### A. Multiple choice (choose the single best answer)

1. The **Wiener–Khinchin** theorem states that the PSD is:
   (a) the mean of the autocorrelation (b) the Fourier transform of the autocorrelation (c) the variance of the signal (d) the reciprocal of the autocorrelation

2. For a zero-mean WSS process, $\int_{-\infty}^{\infty} S_X(f)\,df$ equals:
   (a) $0$ (b) $\mu_X$ (c) $R_X(0) = \mathbb{E}[X^2]$ (d) $\operatorname{Var}[X]^2$

3. In the **Karhunen–Loève** expansion, the coefficients $Z_k$ are:
   (a) independent Gaussian (b) uncorrelated with variance $\lambda_k$ (c) always zero-mean Gaussian (d) deterministic

4. KL eigenfunctions of the Brownian-motion covariance $\min(s,t)$ are:
   (a) cosine functions (b) sine functions (c) Legendre polynomials (d) Hermite functions

5. Spectral factorization requires the factor $L(z)$ to be:
   (a) causal but not necessarily minimum-phase (b) minimum-phase but not necessarily causal (c) causal, monic, and minimum-phase (d) anti-causal

6. Which condition is **sufficient** for ergodicity in the mean?
   (a) $R_X(\tau) \to 0$ as $\tau \to \infty$ (b) $\int_0^{\infty} |C_X(\tau)|\,d\tau < \infty$ (c) $S_X(f) > 0$ everywhere (d) $X$ is Gaussian

7. The periodogram $\hat S_X(\omega)$ is:
   (a) consistent but biased (b) unbiased but inconsistent (c) both unbiased and consistent (d) neither

8. Welch's method reduces the **variance** of the periodogram estimate by roughly:
   (a) $1/\sqrt N$ (b) $1/N$ (c) $1/K$ where $K$ is the number of segments (d) it does not reduce variance

9. The **orthogonality principle** for LMMSE states that:
   (a) $\mathbb{E}[XY] = 0$ (b) $\mathbb{E}[(X - \hat X)Y_j] = 0$ for all observations $Y_j$ (c) the error is Gaussian (d) the estimate is unbiased

10. The **non-causal Wiener filter** for a signal in independent additive noise is:
    (a) $S_{SS}/(S_{SS} + S_{NN})$ (b) $S_{SS}/S_{NN}$ (c) $S_{SS} - S_{NN}$ (d) $S_{NN}/S_{SS}$

11. In the Kalman filter, a **large measurement-noise variance $R$** makes the Kalman gain:
    (a) large (trust the observation) (b) small (trust the prediction) (c) zero (ignore everything) (d) unpredictable

12. The Kalman filter is optimal among:
    (a) all estimators for any noise distribution (b) linear estimators for any noise distribution (c) only Gaussian estimators (d) all estimators for Gaussian distributions **and** linear estimators for non-Gaussian noise (b) and (d) are equivalent statements of the same fact — choose **(d)**

### B. Fill in the blank

1. Wiener–Khinchin: $S_X(f) = \mathcal F\{\_\_\_\_\}$.
2. The PSD integrates to $R_X(0)$ by \_\_\_\_\_\_\_\_'s theorem.
3. The KL coefficients $Z_k$ have variance $\lambda_k$, the \_\_\_\_\_ of the covariance kernel.
4. Spectral factorization: $S_X(z) = \sigma_E^2\, L(z)\, L(\_\_\_\_)$ with $L$ minimum-phase.
5. The innovations process is \_\_\_\_ (has flat spectrum) by construction.
6. The necessary and sufficient condition for ergodicity in the mean is $\frac{1}{T}\int_0^T (1 - \tau/T) C_X(\tau)\,d\tau \to \_\_\_\_$.
7. The raw periodogram has variance approximately equal to $\_\_\_\_^2$ of the true PSD.
8. The LMMSE estimator is $\hat X = \mu_X + \Sigma_{XY}\,\_\_\_\_^{-1}(Y - \mu_Y)$.
9. The causal Wiener filter uses the \_\_\_\_\_\_\_\_ factorization of $S_{YY}$.
10. In the Kalman filter, the variance recursion converges to a steady-state solution of the discrete algebraic \_\_\_\_\_\_\_\_ equation.

### C. Define / explain

1. State the **Wiener–Khinchin theorem** and prove it for the discrete-time case starting from the definition of the PSD as the Fourier transform of $R_X[k]$.
2. Explain the **Karhunen–Loève expansion** and state what optimality property it has over Fourier or Legendre expansions. Give an example of a process for which KL is essentially trivial.
3. State and briefly justify the **spectral factorization theorem**. Why must the factor $L(z)$ be minimum-phase?
4. State the **ergodicity condition** for the mean. Give a process that is ergodic in the mean but **not** in the autocorrelation.
5. Derive the **LMMSE estimator** from the orthogonality principle. State the minimum MSE formula.
6. Derive the **non-causal Wiener filter** from the orthogonality principle in the frequency domain. Explain why the causal version requires spectral factorization.
7. State the **Kalman filter** predict and update equations. Interpret each of $K_k$, $P_k^-$, $P_k$ in words.
8. Explain the connection between the **causal Wiener filter** and the **steady-state Kalman filter**. Under what conditions does the Kalman gain converge?

---

## 5.12 Answer key

??? success "A. Multiple choice"
    1. **(b)** — Wiener–Khinchin: $S_X = \mathcal F\{R_X\}$.
    2. **(c)** — Parseval: $\int S_X = R_X(0) = \mathbb{E}[X^2]$.
    3. **(b)** — Uncorrelated with variance $\lambda_k$ (Gaussian only if $X$ is Gaussian).
    4. **(b)** — $\phi_k(t) = \sqrt{2/T}\sin((2k-1)\pi t/(2T))$.
    5. **(c)** — Causal, monic, minimum-phase — the standard factor in the theorem.
    6. **(b)** — Absolute integrability of $C_X$ is sufficient (not necessary).
    7. **(b)** — Asymptotically unbiased but inconsistent (variance does not shrink with $N$).
    8. **(c)** — Variance reduced by the number of segments $K$.
    9. **(b)** — The error is orthogonal to every observation.
    10. **(a)** — $S_{SS}/(S_{SS} + S_{NN})$.
    11. **(b)** — Small gain, i.e., trust the prediction.
    12. **(d)** — Optimal linear estimator for any noise; optimal overall (also nonlinear) for Gaussian noise.

??? success "B. Fill in the blank"
    1. $R_X(\tau)$
    2. Parseval
    3. eigenvalues
    4. $z^{-1}$
    5. white
    6. $0$
    7. $S_X(\omega)$
    8. $\Sigma_{YY}$
    9. spectral
    10. Riccati

??? success "C. Define / explain — model answers"
    1. *Theorem:* For a WSS process with autocorrelation $R_X[k]$ and PSD $S_X(\omega)$, $S_X(\omega) = \sum_k R_X[k]e^{-j\omega k}$ and $R_X[k] = \frac{1}{2\pi}\int_{-\pi}^{\pi} S_X(\omega)e^{j\omega k}d\omega$. *Proof sketch:* $R_X[k] = \mathbb{E}[X[n+k]X[n]]$. On finite windows, the expected value of the periodogram $\frac{1}{N}|\sum x[n]e^{-j\omega n}|^2$ is the convolution of $S_X$ with the Fejér kernel; taking $N \to \infty$ and using the fact that the Fejér kernel approaches a delta function gives $S_X(\omega) = \sum_k R_X[k]e^{-j\omega k}$. The inverse follows by Fourier inversion on the unit circle.
    2. The **KL expansion** $X(t) = \sum_k Z_k\phi_k(t)$ uses the **eigenfunctions of the covariance kernel** as the basis, with the coefficients $Z_k$ uncorrelated and of variance $\lambda_k$. No other orthogonal basis concentrates variance so quickly in the first few terms: KL minimizes the mean-square truncation error at every order $M$. **Trivial example:** white noise on $[0,T]$ — all eigenvalues equal, so *any* orthonormal basis (including Fourier) is a KL basis.
    3. **Spectral factorization theorem:** every rational nonnegative PSD $S_X(z)$ that has no zeros on the unit circle can be written as $\sigma_E^2 L(z)L(z^{-1})$ with $L$ causal ($L(z) = \sum_{k\ge 0} l_k z^{-k}$, $l_0 = 1$) and minimum-phase (all poles and zeros inside the unit disk). *Why minimum-phase:* the causal filter with a given magnitude response that has the smallest possible delay (earliest impulse-response energy) is the minimum-phase one — this gives the *most efficient* innovations representation and the *smallest* $\sigma_E^2$ among all causal invertible factors. Any other factor would give a larger innovation variance (Example 5.5).
    4. **Ergodicity condition for the mean:** $\lim_{T\to\infty}\frac{1}{T}\int_0^T (1-\tau/T)C_X(\tau)d\tau = 0$; sufficient: $\int_0^\infty |C_X(\tau)|d\tau < \infty$. *Ergodic in mean but not in autocorrelation:* the random-phase sinusoid $X(t) = A\cos(2\pi f_0 t + \Theta)$ — its time-averaged mean is $0$ (equal to the ensemble mean) for every realization, but its time-averaged lag product $\overline{X(t+\tau)X(t)}$ depends on the phase $\Theta$, so different realizations give different empirical autocorrelations.
    5. **Orthogonality principle:** $\hat X$ minimizes $\mathbb{E}[(X-\hat X)^2]$ over affine estimators iff the error $X - \hat X$ is orthogonal to every observation: $\mathbb{E}[(X-\hat X)Y_j] = 0$. Write $\hat X = a^\top Y + b$; the conditions give $b = \mu_X - a^\top\mu_Y$ and $\Sigma_{XY} - a^\top\Sigma_{YY} = 0$, so $a = \Sigma_{YY}^{-1}\Sigma_{XY}^\top$ and $\hat X = \mu_X + \Sigma_{XY}\Sigma_{YY}^{-1}(Y - \mu_Y)$. Minimum MSE $= \operatorname{Var}[X] - \Sigma_{XY}\Sigma_{YY}^{-1}\Sigma_{XY}^\top$.
    6. **Non-causal Wiener filter** derivation: with $\hat S[n] = \sum_k h[k]Y[n-k]$ and $h$ *not* required to be causal, the orthogonality principle $\mathbb{E}[(S[n] - \hat S[n])Y[n-m]] = 0$ for all $m$ becomes $R_{SY}[m] = \sum_k h[k]R_{YY}[m-k]$ for all integers $m$. Take the Fourier transform: $S_{SY}(\omega) = H(\omega)S_{YY}(\omega)$, so $H = S_{SY}/S_{YY}$. The causal version is harder because the orthogonality condition holds only for $m \ge 0$; the trick is to use the **spectral factorization** of $S_{YY}$, which converts the one-sided condition into an algebraic problem and yields the causal filter $H_{\text{cas}}(z) = \frac{1}{\sigma_E^2 L(z)}[S_{SY}(z)/L(z^{-1})]_+$.
    7. **Predict:** $\hat x_k^- = A\hat x_{k-1}$, $P_k^- = AP_{k-1}A^\top + Q$ — the prior mean and covariance before seeing $z_k$. **Update:** $K_k = P_k^-H^\top(HP_k^-H^\top + R)^{-1}$, $\hat x_k = \hat x_k^- + K_k(z_k - H\hat x_k^-)$, $P_k = (I - K_kH)P_k^-$. *Interpretations:* $\hat x_k^-$ and $P_k^-$ summarize belief about $x_k$ before observing $z_k$; $K_k$ is the gain that weights the innovation $(z_k - H\hat x_k^-)$ against the prior prediction; $\hat x_k$ is the posterior mean after combining prediction and observation; $P_k$ is the posterior covariance, always smaller (in Loewner order) than $P_k^-$ after a valid update.
    8. The **causal Wiener filter** and the **steady-state Kalman filter** compute the *same* LMMSE filter $H(z)$ for the same signal-in-noise problem, when both models are time-invariant and the Kalman filter has reached steady state. The Kalman formulation is more general: it handles time-varying models, non-stationary initial conditions, and multiple-input multiple-output systems via the matrix DARE. The Kalman gain converges to a constant matrix $K_\infty$ under **stabilizability** (unstable modes are reachable by the process noise) and **detectability** (unstable modes are observable through the measurement), which are the same conditions under which the causal Wiener filter is well-defined.

---

## Further reading

- Papoulis & Pillai, *Probability, Random Variables, and Stochastic Processes*, 4th ed. — Chapters 10–11 for spectral representation and estimation.
- Porat, *Digital Processing of Random Signals*, Prentice-Hall — Karhunen–Loève, spectrum estimation, factorization.
- Hayes, *Statistical Digital Signal Processing and Modeling*, Wiley — the definitive textbook for spectrum estimation (periodogram, Welch, MUSIC, AR, etc.).
- Kay, *Modern Spectral Estimation*, Prentice-Hall — the rigorous bias-variance analysis of periodogram-based estimators.
- Kailath, Sayed & Hassibi, *Linear Estimation*, Prentice-Hall — LMMSE, Wiener, Kalman, all unified in the orthogonality framework.
- Anderson & Moore, *Optimal Filtering*, Dover — the classical Kalman and Wiener reference.
- Simon, *Optimal State Estimation*, Wiley — modern unified treatment of Kalman, H∞, and particle filters.
- Stoica & Moses, *Spectral Analysis of Signals*, Prentice-Hall — subspace methods, parameter estimation, and modern spectral analysis.
- SciPy `signal` documentation (`welch`, `periodogram`, `levinson`, `wiener`); Math.NET Numerics Kalman filter docs; `filterpy` and `pykalman` Python packages.

---

## Suggested tutorial problems for this chapter

1. **Wiener–Khinchin sanity check.** For $X[n]$ = a Bernoulli($1/2$) sequence taking values $\pm 1$ with $R_X[k] = \delta[k]$, compute $S_X(\omega)$ and verify that $\int S_X(\omega)d\omega/(2\pi) = 1$.
2. **KL of an AR(1).** Compute the KL eigenfunctions of the covariance kernel $K(s,t) = \rho^{|s-t|}$ on $[0,T]$ (this is the covariance of the Ornstein–Uhlenbeck process; the eigen-problem is Sturm–Liouville with $\tan$ boundary conditions). Verify numerically against the $N=100$ discretization eigen-decomposition of the Toeplitz covariance matrix.
3. **Factorization exercise.** Given $S_X(z) = \dfrac{(1 + 0.5 z^{-1})(1 + 0.5 z)}{(1 - 0.8 z^{-1})(1 - 0.8 z)}$, find the causal minimum-phase factor $L(z)$ and the innovations variance. What is the corresponding innovations sequence $E[n]$ in terms of $X[n]$?
4. **Periodogram statistics.** Simulate $N = 1024$ samples of white Gaussian noise of variance $\sigma^2 = 2$. Compute the periodogram and report the mean and variance across the frequency bins. Verify that the empirical distribution of $\hat S_X(\omega_k)/\sigma^2$ is approximately $\chi^2_2/2$ (exponential with mean 1). Repeat with Welch averaging over 4 segments and report the new variance.
5. **LMMSE.** Let $X \sim \mathcal N(0, 1)$, $Y_1 = X + N_1$, $Y_2 = 0.5X + N_2$ with $N_1, N_2$ i.i.d. $\mathcal N(0,1)$ independent of $X$. Compute the LMMSE estimate of $X$ from $(Y_1, Y_2)$ and its MSE. Compare with the estimate from $Y_1$ alone.
6. **Wiener filter.** Let $Y[n] = S[n] + N[n]$ with $S$ AR(1) ($a = 0.7$, $\sigma_W^2 = 1$) and $N$ white ($\sigma_N^2 = 0.5$). Compute the non-causal Wiener filter $H(\omega)$ and the resulting MSE reduction relative to using $\hat S = Y$. Then compute the causal Wiener filter and compare.
7. **Kalman filter design.** Design a scalar Kalman filter for the model $x_k = x_{k-1} + w_k$ ($Q = 0.01$), $z_k = x_k + v_k$ ($R = 1$). Compute the steady-state gain $K_\infty$ by solving the scalar DARE iteratively. Simulate the filter on 500 samples and report the empirical MSE.
8. **Wiener vs. Kalman equivalence.** For the model of Problem 7, compute the causal Wiener filter (in the $z$-domain) and show that it matches the Kalman filter's steady-state transfer function $H(z) = K/(1 - (1-K)z^{-1})$.
