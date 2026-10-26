# Solutions to Chapter 7 — Queuing Theory

[:octicons-arrow-left-24: Back to Chapter 7](../chapters/07-queuing-theory.md)

---

## 7.1 Little's formula

18 patients, average stay 4 days. $\lambda = 18/4 = 4.5$ arrivals/day.

New stay 3 days at same $\lambda$: $L = 4.5\cdot 3 = \mathbf{13.5}$ patients.

---

## 7.2 M/M/1 by hand

$\lambda = 8$ pkt/s, $\mu = 10$ pkt/s, $\rho = 0.8$.

| Quantity | Value |
|---|---|
| $L$ | $0.8/0.2 = 4$ |
| $L_q$ | $0.64/0.2 = 3.2$ |
| $W$ | $1/(10-8) = 0.5$ s |
| $W_q$ | $0.8/2 = 0.4$ s |
| $\pi_0$ | $0.2$ |
| $\Pr\{N \ge 5\}$ | $\rho^5 = 0.328$ |
| $\Pr\{W > 0.5\text{s}\}$ | $e^{-(\mu-\lambda)\cdot 0.5} = e^{-1} \approx 0.368$ |

Checks: $L = \lambda W = 8\cdot 0.5 = 4$ ✓; $L_q = \lambda W_q = 8\cdot 0.4 = 3.2$ ✓.

---

## 7.3 Unstable queue

$\lambda = 12$, $\mu = 10$: $\rho = 1.2 > 1$.

The expression $L = \rho/(1-\rho)$ is negative — **not** the answer: the stationary distribution does not exist. The queue length grows as $(\lambda-\mu)t = 2t$. After $T=10^4$: expected queue length $\approx 2\cdot 10^4$. Simulation shows linear growth with slope $\lambda-\mu = 2$.

---

## 7.4 Erlang-C sensitivity

$\lambda = 50$, $\mu = 20$ calls/min ⇒ $a = 2.5$ Erlangs. Find smallest $s$ with $C(s,a) \le 0.10$.

| $s$ | $\rho$ | $C(s,a)$ | Wait? |
|:---:|:---:|:---:|:---:|
| 3 | 0.833 | 0.578 | 57.8% — too high |
| 4 | 0.625 | 0.164 | 16.4% — too high |
| **5** | **0.500** | **0.036** | **3.6%** ✓ |

At $s=5$: $L_q = C\cdot \rho/(1-\rho) = 0.036\cdot 1 = 0.036$; $W_q = L_q/\lambda = 0.036/50 \approx 0.00072$ min.

At $s=6$: $C(6, 2.5) \approx 0.005$, $W_q \approx 0.0001$ min — 7× lower wait, but 20% more servers. The marginal cost of the extra server is large; **$s=5$ is the practical answer**.

---

## 7.5 M/M/1/K buffer

$\lambda = 5$, $\mu = 10$, $\rho = 0.5$. $\pi_K = \dfrac{(1-\rho)\rho^K}{1-\rho^{K+1}} = \dfrac{0.5^{K+1}}{1-0.5^{K+1}}$.

| $K$ | Blocking $\pi_K$ |
|:---:|:---:|
| 1 | 0.333 |
| 2 | 0.143 |
| 5 | 0.0154 |
| 10 | 0.00049 |

Converges to 0 as $K \to \infty$ — blocking vanishes with an infinite buffer, as expected.

---

## 7.6 PK formula comparison

$\lambda = 6$, $1/\mu = 0.1$ ⇒ $\mu = 10$, $\rho = 0.6$.

$$
L_q = \frac{\rho^2 + \lambda^2\sigma_s^2}{2(1-\rho)} = \frac{0.36 + 36\,\sigma_s^2}{0.8}.
$$

| Service | $\sigma_s^2$ | $L_q$ |
|---|:---:|:---:|
| Deterministic | 0 | $0.36/0.8 = 0.45$ |
| Erlang-2 | $1/(2\mu^2) = 0.005$ | $(0.36+0.18)/0.8 = 0.675$ |
| Exponential | $1/\mu^2 = 0.01$ | $(0.36+0.36)/0.8 = 0.9$ |

**Order:** $L_q^{M/D/1} < L_q^{M/E_2/1} < L_q^{M/M/1}$ ✓. Reducing service variance nearly halves $L_q$ — a cheap win over adding servers.

---

## 7.7 Jackson network

Two queues in series. Queue 1 has external arrivals $r_1 = 4$, $\mu_1 = 10$. Queue 2 has $\mu_2 = 8$, no external arrivals.

**Traffic equations.** $\lambda_1 = 4$; $\lambda_2 = 0 + \lambda_1 = 4$.

**Stability.** $\rho_1 = 0.4 < 1$ ✓; $\rho_2 = 0.5 < 1$ ✓.

**Per-queue means.**

$$
L_1 = \frac{0.4}{0.6} \approx 0.667,\qquad L_2 = \frac{0.5}{0.5} = 1.
$$

**Network.** $L_{\text{net}} = 0.667 + 1 = 1.667$.

**Network arrival rate** $= r_1 = 4$. **Network mean time:** $W_{\text{net}} = L_{\text{net}}/4 \approx 0.417$ time units.

---

## 7.8 Simulation project

$\lambda = 7$, $\mu = 10$. Theory:

- $L = 7/3 \approx 2.333$, $W = 1/3 \approx 0.333$
- $L_q = 1.633$, $W_q = 0.233$

Simulation with $T=10^6$ typically gives $L$ within $\pm 0.05$ of $2.333$ with batch-means 95% CI.

**M/D/1 comparison.** Deterministic service with mean $0.1$:

$$
L_q^{M/D/1} = \frac{\rho^2}{2(1-\rho)} = \frac{0.49}{0.6} \approx 0.817,\qquad L^{M/D/1} \approx 0.817 + 0.7 = 1.517.
$$

That is **35% lower** than M/M/1's $L = 2.333$ — deterministic service is fundamentally more efficient than exponential at the same utilisation.

---

## 7.9 Real-data analysis

On public traces (WIDE / CAIDA), interarrival distributions are **heavier-tailed than exponential**. Log-normal or Weibull fits give a much better description of the tail.

Simulation with the empirical distribution therefore produces **longer** queues than M/M/1 predicts at the same mean arrival rate — the exponential assumption underestimates burstiness. This is why practical network engineering uses:

- **M/G/1** with empirically-measured service variance, or
- **self-similar / heavy-tailed** traffic models (e.g. fractional Brownian motion), or
- **traffic shaping** to make the arrival process *closer* to Poisson (smoothing the burstiness).

M/M/1 is the analytically tractable reference — it is usually an **optimistic** predictor for real traffic.