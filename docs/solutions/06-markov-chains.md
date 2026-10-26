# Solutions to Chapter 6 — Markov Chains

[:octicons-arrow-left-24: Back to Chapter 6](../chapters/06-markov-chains.md)

---

## 6.1 Two-state chain $P^n$

$P = \begin{bmatrix}0.9&0.1\\0.5&0.5\end{bmatrix}$.

Eigenvalues: $1$ and $0.4$. Eigenvectors give $P^n = \mathbf 1\pi + 0.4^n\begin{bmatrix}2/3 & -2/3\\-1/3 & 1/3\end{bmatrix}$ with $\pi = (0.8333, 0.1667)$.

Convergence rate: $0.4^n$. Mixing time to $\|P^n - \mathbf 1\pi\|_\infty < 10^{-3}$: $0.4^n < 10^{-3}$ ⇒ $n \ge 7.5$ ⇒ **$n = 8$**.

---

## 6.2 Gambler's ruin with $p=0.45$

$u_k = \dfrac{1 - (1.2222)^k}{1 - (1.2222)^4}$.

| $k$ | $u_k$ |
|:---:|:---:|
| 1 | 0.205 |
| 2 | 0.395 |
| 3 | 0.606 |

Expected absorption times, solving $E[T_k] = 1 + pE[T_{k+1}] + (1-p)E[T_{k-1}]$ with $E[T_0]=E[T_4]=0$:

$E[T_1] \approx 2.77$, $E[T_2] \approx 3.44$, $E[T_3] \approx 2.77$.

---

## 6.3 Ehrenfest urn ($N=20$)

The chain has period 2 — parity of the number of balls in urn 1 alternates every step. The raw empirical histogram (one lag per step) alternates between two distributions concentrated on even and odd states. **Averaging over even and odd steps separately** gives $\binom{20}{k}/2^{20}$ — the true stationary distribution.

---

## 6.4 Classification table

| Matrix | Classes | Recurrent? | Period | $\pi$ | $P^n$ converges? |
|---|---|---|---|---|---|
| $P_1 = I$ | $\{1\}, \{2\}$ | both absorbing, recurrent | 1 | not unique | yes |
| $P_2 = \begin{bmatrix}0&1\\1&0\end{bmatrix}$ | $\{1,2\}$ | irreducible, recurrent | 2 | $(0.5, 0.5)$ | no (Cesàro average does) |
| $P_3 = \begin{bmatrix}0.5&0.5\\0.5&0.5\end{bmatrix}$ | $\{1,2\}$ | irreducible, recurrent | 1 | $(0.5, 0.5)$ | yes ($P_3^n = P_3$ for $n \ge 1$) |

---

## 6.5 Null recurrence of symmetric random walk on $\mathbb Z$

$f_0 = 1$ (returns a.s.). But the first-return probabilities are

$$
\Pr\{\text{first return at } n\} = \frac{1}{n-1}\binom{n}{n/2}2^{-n}
$$

for even $n = 2m$, which is $\sim \dfrac{1}{n^{3/2}\sqrt{\pi/2}}$. Then

$$
E[T_0] = \sum_n n \cdot \Pr\{\text{first return at } n\} \sim \sum_n n^{-1/2} \to \infty.
$$

So $m_0 = \infty$ — **null recurrent**, not positive recurrent.

---

## 6.6 Birth–death CTMCs

**M/M/1** ($\lambda_k = \lambda$, $\mu_k = \mu$): $\pi_k = (1-\rho)\rho^k$ for $\rho < 1$.

**M/M/∞** ($\lambda_k = \lambda$, $\mu_k = k\mu$): $\pi_k = e^{-a}a^k/k!$ with $a = \lambda/\mu$ — **Poisson**.

The M/M/∞ distribution follows from $\pi_k = \pi_0\prod_{j=0}^{k-1}\lambda/((j+1)\mu) = \pi_0 a^k/k!$.

---

## 6.7 Detailed balance check

Three-state directed cycle with $P_{12}=P_{23}=P_{31}=1$, all other entries zero.

$\pi = (1/3,1/3,1/3)$ satisfies $\pi P = \pi$ ✓ (each state has exactly one incoming and one outgoing edge). But $\pi_1 P_{12} = 1/3 \ne \pi_2 P_{21} = 0$ — **detailed balance fails**. This chain is not reversible: forward time goes $1\to2\to3\to1$, backward gives $1\to3\to2\to1$, a different chain.

---

## 6.8 M/M/1 simulation

$T = 10^5$, $\lambda=0.8$, $\mu=1.0$: empirical $\bar L \approx 4$ (theory $\rho/(1-\rho) = 4$). Empirical distribution of $K$ matches $(1-\rho)\rho^k = 0.2\cdot 0.8^k$.

With $\lambda = 1.1$: $\rho = 1.1 > 1$. The queue grows as $(\lambda-\mu)t = 0.1t$; empirical histogram never stabilizes — this **is** what instability looks like: no stationary distribution exists.

---

## 6.9 PageRank preview

The PageRank chain has transition

$$
P_{ij} = \frac{1-d}{N} + d\,\frac{a_{ij}}{L_i},\qquad d = 0.85,
$$

where $a_{ij}\in\{0,1\}$ indicates a link from $i$ to $j$ and $L_i$ is the out-degree of $i$.

Solving $\pi P = \pi$ gives the PageRank vector. For a 4-page graph with a mix of incoming links, $\pi$ values are on the order of $1/N$ to $3/N$, weighted by link count.

**Why $d < 1$ matters:** without damping ($d=1$), a page with no outgoing links (a "sink") or disconnected subgraph makes the chain reducible or periodic — $\pi$ may not exist or may not be unique. The damping term guarantees irreducibility and aperiodicity, so $\pi$ exists and is unique.

---

## 6.10 (extra) Two-state chain eigen-decomposition for mixing time

$P = \begin{bmatrix}0.9&0.1\\0.5&0.5\end{bmatrix}$: subdominant eigenvalue $0.4$. The $\ell_\infty$ distance from stationarity decays as $0.4^n$. To reach $10^{-3}$: $n \ge 8$ (see 6.1).