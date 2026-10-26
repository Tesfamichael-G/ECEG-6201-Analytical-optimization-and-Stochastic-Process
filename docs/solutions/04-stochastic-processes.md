# Solutions to Chapter 4 — Stochastic Processes

[:octicons-arrow-left-24: Back to Chapter 4](../chapters/04-stochastic-processes.md)

---

## 4.1 WSS check

- $X_1(t) = A\cos(\omega_0 t + \Theta)$, $\Theta$ uniform: mean $0$, $R_X(\tau) = (A^2/2)\cos(\omega_0\tau)$ → **WSS**.
- $X_2(t) = A\cos(\omega_0 t) + B$, $A,B,\Theta$ independent: $\mathbb E[X_2] = \mathbb E[B]$, $R_{X_2}$ depends only on $\tau$ → **WSS**.
- $X_3(t) = W(t)\cos(\omega_0 t)$, $W$ zero-mean WSS: $R_{X_3}(t+\tau,t) = R_W(\tau)\cos(\omega_0(t+\tau))\cos(\omega_0 t)$ — depends on $t$ → **not WSS**, but **cyclostationary** with period $T_0/2 = \pi/\omega_0$.

---

## 4.2 Random-walk moments

$S_n = \sum_{k=1}^n X_k$, $X_k \in \{\pm 1\}$, $\Pr\{X_k=+1\}=p$.

$\mathbb E[S_n] = n(2p-1)$; $\operatorname{Var}[S_n] = 4p(1-p)n$.

For $p \ne 1/2$: $S_n/n \to 2p-1$ a.s. (SLLN). So $|S_n| \to \infty$ a.s. — the walk drifts to $\pm\infty$; it is **transient**.

---

## 4.3 Poisson superposition

PGFs: $G_{N_1}(z) = e^{\lambda_1 t(z-1)}$, $G_{N_2}(z) = e^{\lambda_2 t(z-1)}$.

For $N = N_1+N_2$ independent: $G_N = G_{N_1}G_{N_2} = e^{(\lambda_1+\lambda_2)t(z-1)}$ — the PGF of Poisson($(\lambda_1+\lambda_2)t$). ✓

---

## 4.4 Shot noise with rectangular pulse

$h(s) = 1$ for $0 \le s \le T_p$, $0$ otherwise.

**Mean.** $\mathbb E[X] = \lambda\int_0^{T_p} ds = \lambda T_p$.
**Variance.** $\operatorname{Var}[X] = \lambda\int_0^{T_p} ds = \lambda T_p$.
**Autocorrelation.** $R_X(\tau) = \lambda(T_p - |\tau|)^+ + (\lambda T_p)^2$ — a triangular term plus a constant. Correlation time: $T_p$.

---

## 4.5 Cyclostationary PAM

$p(t) = \sqrt{2/T}\cos(2\pi t/T)\cdot\mathbf 1_{|t|\le T/2}$.

$\bar R_X(\tau) = \dfrac{\sigma_a^2}{T}\int p(s+\tau)\,p(s)\,ds$ — a raised-cosine-shaped autocorrelation. Its Fourier transform is the **raised-cosine spectrum** with bandwidth $W = 1/T$; this is the pulse shape that satisfies the **Nyquist zero-ISI criterion** and is standard in digital communications.

---

## 4.6 Sampling of a band-limited Gaussian process

$W = 10$ Hz, sampled at $20$ Hz (Nyquist): sinc reconstruction error $\sim 10^{-14}$ — machine precision.

Sampled at $15$ Hz (below Nyquist): spectrum above $7.5$ Hz folds into lower frequencies — aliasing. Reconstruction error depends on the power above $7.5$ Hz; for a typical smooth spectrum, 5–15% of total power is aliased.

---

## 4.7 Matched filter for a triangular pulse

$s(t) = 1-|t|/T$ on $[-T, T]$. Energy $E_s = 2T/3$.

Matched filter $h(t) = s(T-t)$; maximum SNR $= 2E_s/N_0 = 4T/(3N_0)$.

**Comparison with a rectangular pulse of the same energy:** the matched-filter SNR depends only on $E_s$, so the two are identical — shape does not affect the maximum achievable SNR. Shape affects robustness to timing errors and spectral containment, not the peak SNR.

---

## 4.8 Entropy chain

**(a) Discrete:** Fair coin $H = 1$ bit; biased $p=0.1$ → $H = -0.1\log_2(0.1) - 0.9\log_2(0.9) \approx 0.469$ bits; degenerate (p=0 or 1) → $H = 0$.

**(b) Differential:** $\mathcal N(0,1)$ → $h = \tfrac12\log(2\pi e) \approx 1.419$ nats. $\mathcal N(0,4)$ → $h = \tfrac12\log(8\pi e) \approx 2.305$ nats. Uniform on $[-1,1]$ (variance $= 1/3$) → $h = \log 2 \approx 0.693$ nats.

**Gaussian maximizes differential entropy for a given variance** — the uniform's entropy at variance $1/3$ is smaller than the Gaussian's entropy at variance $1/3$ (which would be $\approx 0.919$ nats).