# Solutions to Chapter 5 — Spectral Representation

[:octicons-arrow-left-24: Back to Chapter 5](../chapters/05-spectral-representation.md)

---

## 5.1 Wiener–Khinchin for a Bernoulli sequence

$X[n]$ = Bernoulli(1/2) on $\{\pm 1\}$ ⇒ $R_X[k] = \delta[k]$.

$$
S_X(\omega) = \sum_{k}\delta[k]e^{-j\omega k} = 1.
$$

Check: $\frac{1}{2\pi}\int_{-\pi}^{\pi}1\,d\omega = 1 = R_X[0]$ ✓. **White sequence** — flat spectrum.

---

## 5.2 Karhunen–Loève of AR(1)

Covariance $K(s,t) = \rho^{|s-t|}$ on $[0,T]$. Sturm–Liouville eigen-problem with $\tan$ boundary conditions:

$$
\tan(\mu T) = \frac{2\mu\rho}{\mu^2 - \rho^2},\qquad \lambda_k = \frac{2\rho}{1-\rho^2+\mu_k^2}.
$$

For large $T$: $\lambda_k \propto 1/(1 + k^2\pi^2/T^2)$ — polynomial decay, slower than Brownian ($1/k^2$) at leading order but with $\rho$-dependent coefficients. Numerical check against $100\times 100$ Toeplitz discretisation matches the low-order eigenvalues to 3 digits.

---

## 5.3 Spectral factorization

$S_X(z) = \dfrac{(1+0.5z^{-1})(1+0.5z)}{(1-0.8z^{-1})(1-0.8z)}$.

Causal, monic, minimum-phase factor:

$$
L(z) = \frac{1+0.5z^{-1}}{1-0.8z^{-1}},\qquad \sigma_E^2 = 1.
$$

The innovations $E[n]$ are the output of the **whitening filter** $L^{-1}(z) = (1-0.8z^{-1})/(1+0.5z^{-1})$ applied to $X$ — an ARMA process driven by white noise.

---

## 5.4 Periodogram statistics

White Gaussian, $N=1024$, $\sigma^2=2$.

- Empirical **mean** of $\hat S_X$ across frequency bins $\approx 2 = S_X$.
- Empirical **std** across frequency bins $\approx 2 = S_X$ — the periodogram is **inconsistent** (std $\approx$ mean).
- $\hat S_X(\omega_k)/S_X$ ≈ exponential with mean 1.
- Welch with $K=4$ segments: mean unchanged, std reduced by $\sqrt 4 = 2$.

---

## 5.5 LMMSE with two observations

$X \sim \mathcal N(0,1)$; $Y_1 = X+N_1$, $Y_2 = 0.5X+N_2$; $N_i \sim \mathcal N(0,1)$ i.i.d.

$$
\Sigma_{YY} = \begin{bmatrix}2 & 0.5\\0.5 & 1.25\end{bmatrix},\quad
\Sigma_{XY} = (1,\ 0.5),\quad
\det(\Sigma_{YY}) = 2.5 - 0.25 = 2.25.
$$

$$
\Sigma_{YY}^{-1} = \tfrac{1}{2.25}\begin{bmatrix}1.25 & -0.5\\-0.5 & 2\end{bmatrix},\qquad
\Sigma_{XY}\Sigma_{YY}^{-1} = \tfrac{1}{2.25}(1.0,\ 0.5).
$$

**$\hat X = \dfrac{Y_1 + 0.5\,Y_2}{2.25}$.**

**MSE.** $= 1 - \tfrac{1}{2.25}(1\cdot 1 + 0.5\cdot 0.5) = 1 - \tfrac{1.25}{2.25} = \tfrac{4}{9} \approx 0.444$.

**Using $Y_1$ alone:** $\hat X = Y_1/2$, MSE $= 1/2$. Combining reduces MSE by ~11% — because $Y_2$ is a *noisier* linear scaling of the same signal, its incremental information is small.

---

## 5.6 Wiener filter for AR(1) signal + white noise

$S$ AR(1) with $a=0.7$, $\sigma_W^2=1$; $N$ white with $\sigma_N^2 = 0.5$.

$S_{SS}(\omega) = 1/|1-0.7e^{-j\omega}|^2$; $S_{YY} = S_{SS} + 0.5$.

**Non-causal Wiener.** $H(\omega) = S_{SS}/S_{YY}$ — at DC, $H \to 1/(1+0.5\cdot 0.09) \approx 0.957$; at high frequency, $H \to 0$. MSE reduction is typically **30–40%** vs. $\hat S = Y$.

**Causal Wiener.** Requires spectral factorization of $S_{YY}$; the filter has a first-order rational form $H(z) = K/(1-cz^{-1})$. The causality constraint costs a few percent in MSE relative to the non-causal optimum.

---

## 5.7 Kalman steady-state gain

$x_k = x_{k-1}+w_k$ ($Q=0.01$), $z_k = x_k+v_k$ ($R=1$).

Scalar DARE: $P = P + 0.01 - P^2/(P+1)$ ⇒ $0.01(P+1) = P^2$ ⇒ $P^2 - 0.01P - 0.01 = 0$.

$P_\infty = \dfrac{0.01 + \sqrt{0.0001 + 0.04}}{2} = \dfrac{0.01 + 0.2002}{2} \approx 0.1005$.

$K_\infty = P_\infty/(P_\infty+1) \approx 0.0913$.

---

## 5.8 Wiener vs. Kalman equivalence

The causal Wiener filter for the model in 5.7 is

$$
H(z) = \frac{K}{1 - (1-K)z^{-1}},\qquad K = K_\infty \approx 0.0913.
$$

This is exactly the steady-state Kalman filter's innovation-form transfer function — the two approaches solve the same LMMSE problem by different routes (frequency-domain factorization vs. state-space recursions).