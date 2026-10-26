# Chapter 6 — Markov Chains

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 tutorial · **Home study:** 5 h
    
    **Coverage (syllabus):** concepts of the Markov process · discrete-time
    Markov chains (DTMC) · continuous-time Markov chains (CTMC) · classification
    of states · stationary distributions · ergodic properties and limiting
    probabilities.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. State the **Markov property** and decide whether a given process satisfies it.
2. Build a **transition matrix** for a DTMC; verify that its rows sum to one.
3. Apply the **Chapman–Kolmogorov equation** to compute $n$-step transition probabilities.
4. Classify states as **recurrent/transient**, **periodic/aperiodic**, and identify **communicating classes**.
5. Compute **stationary distributions** for small DTMCs, both by solving $\pi P = \pi$ and by detailed balance when applicable.
6. State the **ergodic theorem** for finite irreducible aperiodic chains and use it to justify long-run averages.
7. Write the **generator matrix** $Q$ for a CTMC, derive the **Kolmogorov forward equations**, and compute stationary distributions by solving $\pi Q = 0$.
8. Model simple engineering systems (buffers, servers, two-state links) as birth–death chains and derive their stationary distributions.
9. Implement DTMC / CTMC computations in Python, Java, or C#.

---

## 6.1 The Markov property

!!! abstract "Markov property"
    A stochastic process $\{X_n\}_{n \ge 0}$ taking values in a countable set
    $\mathcal S$ (the **state space**) satisfies the **Markov property** if, for
    every $n$ and every choice of states $i_0, \dots, i_n, j$,

    $$
    \Pr\{X_{n+1} = j \mid X_0 = i_0, \dots, X_n = i_n\}
    = \Pr\{X_{n+1} = j \mid X_n = i_n\}.
    $$

    In words: the **future depends on the past only through the present** —
    “the process has no memory of how it got here”.

A process satisfying this property is called a **Markov chain**. When $n \in \mathbb Z_{\ge 0}$ it is a **discrete-time Markov chain (DTMC)**; when $t \in \mathbb R_{\ge 0}$ it is a **continuous-time Markov chain (CTMC)**.

The quantity

$$
P_{ij} = \Pr\{X_{n+1} = j \mid X_n = i\}
$$

is the **transition probability** from state $i$ to state $j$. The chain is **time-homogeneous** if $P_{ij}$ does not depend on $n$; in that case the entire chain is determined by the **transition matrix** $P = (P_{ij})$.

**Basic properties of $P$.**

- $P_{ij} \ge 0$ for all $i, j$.
- $\sum_j P_{ij} = 1$ for each $i$ (each row sums to one) — such a matrix is called **stochastic**.
- The initial distribution $\mu_0$ is a row vector with $\mu_0(i) = \Pr\{X_0 = i\}$, $\sum_i \mu_0(i) = 1$.

### Worked example 6.1 — the weather chain

States: $1$ = sunny, $2$ = rainy. Transition matrix

$$
P = \begin{bmatrix} 0.9 & 0.1 \\ 0.5 & 0.5 \end{bmatrix}.
$$

From sunny: 90% sunny tomorrow, 10% rainy. From rainy: 50% sunny, 50% rainy.

**Two-step transition.** Today sunny, what is the probability of being sunny in 2 days?

$$
\Pr\{X_2 = \text{sunny} \mid X_0 = \text{sunny}\} = (P^2)_{11}.
$$

Compute $P^2$:

$$
P^2 = \begin{bmatrix} 0.9\cdot 0.9 + 0.1\cdot 0.5 & 0.9\cdot 0.1 + 0.1\cdot 0.5 \\ 0.5\cdot 0.9 + 0.5\cdot 0.5 & 0.5\cdot 0.1 + 0.5\cdot 0.5 \end{bmatrix}
= \begin{bmatrix} 0.86 & 0.14 \\ 0.70 & 0.30 \end{bmatrix}.
$$

So the answer is $0.86$. In general, $P^n$ gives $n$-step transition probabilities: $(P^n)_{ij} = \Pr\{X_n = j \mid X_0 = i\}$.

**Long-run behaviour.** Compute $P^\infty = \lim_{n\to\infty} P^n$ by repeatedly squaring $P$:

$$
P^2 = \begin{bmatrix}0.86 & 0.14\\0.70&0.30\end{bmatrix}, \quad
P^4 = \begin{bmatrix}0.842 & 0.158\\0.79&0.21\end{bmatrix}, \quad
P^8 \approx \begin{bmatrix}0.8337 & 0.1663 \\ 0.8315 & 0.1685\end{bmatrix}.
$$

The rows converge to the same vector, roughly $(0.8333,\ 0.1667)$. This row is the **stationary distribution** of the chain — computed in Example 6.4.

!!! tip "Markov property check"
    A process satisfies the Markov property if its future is *conditionally
    independent* of the past given the present. Note that this does **not**
    mean it is memoryless in the everyday sense — the *state* can encode memory
    (e.g. a queue's state encodes how full it is). It means the state variable
    is *rich enough* that knowing it is sufficient to predict the future.

---

## 6.2 The Chapman–Kolmogorov equations

The **Chapman–Kolmogorov (CK) equations** relate the $(m+n)$-step transition probabilities to the $m$-step and $n$-step probabilities:

$$
\boxed{\;(P^{m+n})_{ij} = \sum_{k \in \mathcal S} (P^m)_{ik}\,(P^n)_{kj}\;}
$$

This is exactly matrix multiplication: $P^{m+n} = P^m P^n$. The proof is a one-line application of the **law of total probability** and the Markov property:

$$
\Pr\{X_{m+n} = j \mid X_0 = i\} = \sum_k \Pr\{X_{m+n} = j \mid X_m = k\}\Pr\{X_m = k \mid X_0 = i\}.
$$

The power $P^n$ therefore encodes all multi-step transition probabilities, and the marginal distribution at step $n$ is

$$
\mu_n = \mu_0 P^n,
$$

where $\mu_n(j) = \Pr\{X_n = j\}$ and $\mu_0$ is the initial distribution as a row vector.

### Worked example 6.2 — gambler's ruin (absorbing chain)

A gambler starts with $k \in \{1, 2, 3\}$ birr. Each round, with probability $p$ he wins one birr; with probability $1 - p$ he loses one. He stops when he reaches $\{0, 4\}$ (ruin or goal).

State space $\mathcal S = \{0, 1, 2, 3, 4\}$ with $0$ and $4$ **absorbing** ($P_{00} = P_{44} = 1$). Transition matrix:

$$
P = \begin{bmatrix}
1 & 0   & 0   & 0   & 0 \\
1-p & 0 & p & 0 & 0 \\
0 & 1-p & 0 & p & 0 \\
0 & 0 & 1-p & 0 & p \\
0 & 0   & 0   & 0   & 1
\end{bmatrix}.
$$

**Absorption probabilities.** Let $u_k = \Pr\{\text{reach } 4 \text{ before } 0 \mid X_0 = k\}$. Condition on the first step (first-step analysis):

$$
u_k = p\, u_{k+1} + (1-p)\, u_{k-1},\quad 1 \le k \le 3,\qquad u_0 = 0,\ u_4 = 1.
$$

This is a linear second-order difference equation. Rearranged:

$$
p\, u_{k+1} - u_k + (1-p)\, u_{k-1} = 0.
$$

The characteristic equation is $p r^2 - r + (1-p) = 0$, whose roots are $r = 1$ and $r = (1-p)/p$. For $p \ne 1/2$,

$$
u_k = A + B\left(\frac{1-p}{p}\right)^{k}.
$$

Boundary conditions give

$$
u_k = \frac{1 - \left(\frac{1-p}{p}\right)^{k}}{1 - \left(\frac{1-p}{p}\right)^{4}}.
$$

For $p = 1/2$ (fair game), the double root at $r = 1$ gives $u_k = k/4$.

**Sanity check.** For a fair game and $k = 2$: $u_2 = 2/4 = 0.5$ — starting from the middle, ruin and success are equally likely. ✓

**Numerical evaluation** for $p = 0.4$ (unfair against the gambler):

| $k$ | $u_k$ | Probability of ruin |
|:---:|:---:|:---:|
| 1 | 0.1525 | 0.8475 |
| 2 | 0.3295 | 0.6705 |
| 3 | 0.5368 | 0.4632 |

Even with a small disadvantage ($p = 0.4$ vs. $0.5$), the gambler's chance of reaching the goal drops dramatically.

### Worked example 6.3 — Ehrenfest urn (periodic chain)

Two urns contain a total of $N$ balls. At each step, pick a ball uniformly at random and move it to the other urn. Let $X_n$ be the number of balls in urn 1.

**Transition probabilities.** From state $k$ (urn 1 has $k$ balls, urn 2 has $N - k$):

- Move to $k - 1$: pick a ball from urn 1, probability $k/N$.
- Move to $k + 1$: pick a ball from urn 2, probability $(N - k)/N$.
- All other transitions have probability zero.

So $P_{k, k-1} = k/N$, $P_{k, k+1} = (N-k)/N$.

**Periodicity.** Every transition changes $k$ by exactly $\pm 1$, so the parity of $k$ alternates each step. Hence the chain is **periodic with period 2**.

**Stationary distribution.** Try the binomial $\pi(k) = \binom{N}{k} 2^{-N}$. Check detailed balance:

$$
\pi(k) P_{k, k+1} = \binom{N}{k}2^{-N} \cdot \frac{N - k}{N}.
$$

$$
\pi(k+1) P_{k+1, k} = \binom{N}{k+1}2^{-N} \cdot \frac{k+1}{N}.
$$

Using $\binom{N}{k+1}\cdot(k+1) = \binom{N}{k}\cdot(N-k)$:

$$
\pi(k+1) P_{k+1, k} = \binom{N}{k}2^{-N}\cdot\frac{N-k}{N} = \pi(k) P_{k, k+1}. \;\checkmark
$$

So $\pi(k) = \binom{N}{k}/2^N$ is the stationary distribution. It is also the **limiting time-average**: the long-run fraction of time spent at state $k$ is $\binom{N}{k}/2^N$, which is $\mathrm{Binomial}(N, 1/2)$.

!!! tip "Periodicity and the limit"
    For the Ehrenfest chain, $P^n$ does **not** converge as $n \to \infty$ — even and odd $n$
    give different distributions (from parity). But the **Cesàro average**
    $\frac{1}{n}\sum_{k=0}^{n-1} P^k$ **does** converge to the stationary distribution.
    This is why the ergodic theorem for periodic chains refers to *time averages*,
    not to $\lim P^n$ — a distinction Example 6.4 develops.

---

## 6.3 Classification of states

### 6.3.1 Communication and irreducibility

State $j$ is **accessible** from $i$ (written $i \to j$) if $P^n_{ij} > 0$ for some $n \ge 0$. States $i$ and $j$ **communicate** (written $i \leftrightarrow j$) if $i \to j$ and $j \to i$. Communication is an equivalence relation; its equivalence classes are the **communicating classes**.

A chain is **irreducible** if it has a single communicating class (all states communicate).

### 6.3.2 Recurrence and transience

Let $f_i = \Pr\{\text{return to } i \text{ at some } n \ge 1 \mid X_0 = i\}$ be the **return probability**. Then:

- **Recurrent** state: $f_i = 1$ (returns almost surely).
- **Transient** state: $f_i < 1$ (positive probability of never returning).

For a recurrent state, the **mean return time** is $m_i = \mathbb E[T_i \mid X_0 = i]$ where $T_i$ is the first return time. If $m_i < \infty$, the state is **positive recurrent**; if $m_i = \infty$, it is **null recurrent**.

!!! note "Finite-state fact"
    In a **finite** state space, every recurrent state is positive recurrent, and
    there are no null-recurrent states.

### 6.3.3 Periodicity

The **period** of state $i$ is $d(i) = \gcd\{n \ge 1 : P^n_{ii} > 0\}$. If $d(i) = 1$, the state is **aperiodic**. In an irreducible chain all states have the same period, so we can speak of "the period of the chain".

### Worked example 6.4 — a finite irreducible chain with no absorbing states

A chain on $\{1,2,3\}$ with

$$
P = \begin{bmatrix} 0.5 & 0.5 & 0 \\ 0.3 & 0.4 & 0.3 \\ 0 & 0.5 & 0.5\end{bmatrix}.
$$

**Irreducibility.** All entries reachable in $\le 2$ steps; all states communicate. Irreducible. ✓

**Periodicity.** $P_{11} = 0.5 > 0$, so $d(1) = 1$. Aperiodic.

**Stationary distribution.** Solve $\pi P = \pi$, $\sum \pi_i = 1$:

$$
\pi_1 = 0.5\pi_1 + 0.3\pi_2,\quad
\pi_2 = 0.5\pi_1 + 0.4\pi_2 + 0.5\pi_3,\quad
\pi_3 = 0.3\pi_2 + 0.5\pi_3.
$$

From the first: $0.5\pi_1 = 0.3\pi_2$, so $\pi_2 = \tfrac{5}{3}\pi_1$. From the third: $0.5\pi_3 = 0.3\pi_2$, so $\pi_3 = 0.6\pi_2 = \pi_1$. Normalizing: $\pi_1(1 + 5/3 + 1) = 1$, so $\pi_1 = 3/11$ and

$$
\pi = \left(\tfrac{3}{11},\ \tfrac{5}{11},\ \tfrac{3}{11}\right) \approx (0.2727,\ 0.4545,\ 0.2727).
$$

**Long-run behaviour.** Because the chain is irreducible, aperiodic, and finite, $P^n \to \mathbf 1\pi$ as $n \to \infty$ — every row of $P^\infty$ equals $\pi$.

### Worked example 6.5 — periodicity in the two-state chain

A two-state chain with

$$
P = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}
$$

has $P^n = P$ if $n$ is odd, $P^n = I$ if $n$ is even. Clearly $P^n$ does not converge. This chain has period $2$. The **Cesàro average** converges:

$$
\frac{1}{n}\sum_{k=0}^{n-1} P^k \to \begin{bmatrix}0.5 & 0.5\\0.5 & 0.5\end{bmatrix}
$$

which is $\mathbf 1\pi$ with $\pi = (0.5, 0.5)$.

### Summary table of classification

| Property | Definition | Consequence |
|---|---|---|
| Irreducible | One communicating class | Convergence to a unique stationary distribution |
| Positive recurrent | Returns a.s., finite mean return time | Stationary exists with $\pi_i = 1/m_i$ |
| Aperiodic | Period 1 | $\lim P^n = \mathbf 1\pi$ |
| Ergodic | Irreducible + positive recurrent + aperiodic | Long-run time average = ensemble average |

---

## 6.4 Stationary distributions and ergodicity

### 6.4.1 Stationary distribution

A probability vector $\pi$ is a **stationary distribution** of the chain if

$$
\pi P = \pi,\qquad \pi_i \ge 0,\qquad \sum_i \pi_i = 1.
$$

This is the left eigenvector of $P$ at eigenvalue 1. If $\mu_0 = \pi$, then $\mu_n = \pi$ for all $n$ — the chain is stationary in the sense of §4.1.

### 6.4.2 Ergodicity theorem

!!! abstract "Ergodic theorem for DTMC"
    Let $\{X_n\}$ be an irreducible, positive recurrent, aperiodic chain on a
    countable state space with stationary distribution $\pi$. Then:

    1. **Convergence to stationarity:** $P^n \to \mathbf 1\pi$ as $n \to \infty$,
       i.e. $\Pr\{X_n = j \mid X_0 = i\} \to \pi_j$ for all $i, j$, independent of the initial state.
    2. **Time average = ensemble average:** for any bounded function $f$ on the state space,

    $$
    \frac{1}{N}\sum_{n=1}^{N} f(X_n) \xrightarrow{\text{a.s.}} \sum_j \pi_j f(j) \quad \text{as } N \to \infty.
    $$

This is the **strong law of large numbers for Markov chains** — the analogue of the ergodic theorem of Chapter 5 in discrete state space, and the theoretical basis for simulating Markov chains (MCMC).

For **periodic** chains, statement 1 fails but statement 2 still holds with $\lim_{N} \frac{1}{N}\sum f(X_n) \to \sum \pi_j f(j)$ — the time average converges, even though the distribution does not.

### 6.4.3 Detailed balance (reversibility)

A probability vector $\pi$ satisfies **detailed balance** (equivalently, the chain is **time-reversible** with respect to $\pi$) if

$$
\boxed{\;\pi_i P_{ij} = \pi_j P_{ji} \quad \text{for all } i, j.\;}
$$

Detailed balance implies $\pi P = \pi$ (sum over $i$: $\sum_i \pi_i P_{ij} = \pi_j \sum_i P_{ji} = \pi_j$). It is **stronger** than stationarity — many chains have a stationary distribution but are *not* reversible. When it holds, it provides the easiest way to compute $\pi$: just find any positive vector satisfying detailed balance, then normalize.

### Worked example 6.6 — detailed balance for a birth–death chain

Consider a chain on $\{0, 1, 2, \dots\}$ with $P_{k, k+1} = p_k$, $P_{k, k-1} = q_k$ ($p_k + q_k = 1$, $p_k, q_k > 0$). This is a **birth–death chain** — the discrete-time analogue of the M/M/1 queue of Chapter 7.

**Detailed balance.** For a birth–death chain, transitions can only happen between adjacent states, so detailed balance is a *local* condition between $k$ and $k+1$:

$$
\pi_k\, p_k = \pi_{k+1}\, q_{k+1} \;\Longrightarrow\; \pi_{k+1} = \pi_k\, \frac{p_k}{q_{k+1}}.
$$

By induction from $\pi_0$:

$$
\pi_k = \pi_0 \prod_{j=0}^{k-1} \frac{p_j}{q_{j+1}}.
$$

Normalizing, $\pi_0 = \Bigl(1 + \sum_{k=1}^\infty \prod_{j=0}^{k-1} \frac{p_j}{q_{j+1}}\Bigr)^{-1}$ — the chain is positive recurrent iff this sum converges.

**Special case: reflecting random walk on $\{0, 1, \dots, N\}$** with $p_k = p$, $q_k = 1-p$ for $1 \le k \le N-1$:

$$
\pi_k = \pi_0\, \left(\frac{p}{1-p}\right)^{k}.
$$

This is a truncated geometric distribution. If $p < 1/2$ (drift toward 0), $\pi$ concentrates on small $k$; if $p > 1/2$ (drift toward $N$), $\pi$ concentrates on large $k$.

---

## 6.5 Continuous-time Markov chains

### 6.5.1 Definition and generator

A **CTMC** $\{X(t)\}_{t \ge 0}$ on a countable state space satisfies the continuous-time analogue of the Markov property: for any $0 \le t_1 < \dots < t_n < t$,

$$
\Pr\{X(t) = j \mid X(t_n) = i_n, \dots, X(t_1) = i_1\} = \Pr\{X(t) = j \mid X(t_n) = i_n\}.
$$

The chain is **time-homogeneous** if the right side depends only on $t - t_n$. The transition probabilities

$$
P_{ij}(t) = \Pr\{X(t+s) = j \mid X(s) = i\}
$$

form the **transition matrix** $P(t) = (P_{ij}(t))$. It satisfies:

- $P(0) = I$ (identity),
- $P(s + t) = P(s)P(t)$ (semigroup property — the continuous analogue of Chapman–Kolmogorov),
- $P(t)$ is stochastic for each $t$.

The **generator (infinitesimal generator) matrix** is

$$
\boxed{\;Q = \lim_{h \to 0^+} \frac{P(h) - I}{h} = P'(0)\;}
$$

with entries

$$
Q_{ij} = \begin{cases} q_{ij} \ge 0 & \text{(rate of transition } i \to j, \; i \ne j), \\ -\sum_{j \ne i} q_{ij} & \text{(i.e. } Q_{ii} = -q_i \text{, where } q_i \text{ is the total exit rate from } i\text{)}. \end{cases}
$$

Each row of $Q$ sums to zero: $\sum_j Q_{ij} = 0$.

**Relationship to $P(t)$.** The **Kolmogorov forward equation** is

$$
\frac{d}{dt} P(t) = P(t)\, Q \qquad (\text{forward}),
$$

and the **Kolmogorov backward equation** is

$$
\frac{d}{dt} P(t) = Q\, P(t) \qquad (\text{backward}).
$$

Both have the same solution: $P(t) = e^{Qt}$ (matrix exponential). In practice you rarely compute $e^{Qt}$; you solve for the stationary distribution or for specific time-dependent quantities.

### 6.5.2 Holding times and jumps

In a CTMC, once in state $i$, the process stays there for a random time $\tau_i$ that is **exponentially distributed** with rate $q_i = -Q_{ii}$:

$$
\Pr\{\tau_i > t\} = e^{-q_i t}.
$$

At the end of the holding time, the process **jumps** to state $j \ne i$ with probability $q_{ij}/q_i$, independent of the holding time and of the past. This *exponential-holding-time + jump-chain* description is the **uniformization / embedded chain** representation, and it is the basis for simulation algorithms.

### 6.5.3 Stationary distribution

A probability vector $\pi$ is stationary if $\pi P(t) = \pi$ for all $t \ge 0$, which is equivalent to

$$
\boxed{\;\pi Q = 0 \quad \text{with } \sum_j \pi_j = 1\;}
$$

Note $\pi Q = 0$ is **not** a linear system in $\pi$ — its rank is $n - 1$ (because the rows of $Q$ sum to zero), so you supplement it with the normalization $\sum \pi_j = 1$.

### Worked example 6.7 — two-state CTMC

Model a communication link as up/down:

- Up state (0): fails at rate $\alpha$ (goes to Down).
- Down state (1): repairs at rate $\beta$ (goes to Up).

Generator:

$$
Q = \begin{bmatrix} -\alpha & \alpha \\ \beta & -\beta \end{bmatrix}.
$$

**Stationary distribution.** Solve $\pi Q = 0$:

$$
-\alpha \pi_0 + \beta \pi_1 = 0 \;\Longrightarrow\; \pi_1 = \frac{\alpha}{\beta}\,\pi_0.
$$

With normalization:

$$
\pi_0 = \frac{\beta}{\alpha + \beta},\qquad \pi_1 = \frac{\alpha}{\alpha + \beta}.
$$

**Interpretation.** The **availability** of the link is $\pi_0 = \beta/(\alpha+\beta)$ — the ratio of repair rate to total rate. If $\beta \gg \alpha$ (fast repair), the link is up almost all the time.

| $\alpha$ (MTTF⁻¹) | $\beta$ (MTTR⁻¹) | $\pi_0$ |
|:---:|:---:|:---:|
| 0.01 | 1 | 0.990 |
| 0.1 | 1 | 0.909 |
| 0.1 | 0.1 | 0.500 |

**Embedded chain.** The jump chain has transitions $0 \to 1$ with probability 1 and $1 \to 0$ with probability 1 — it alternates deterministically. All randomness in the CTMC is in the *timing*, not the *sequence* of states.

### Worked example 6.8 — birth–death CTMC and M/M/1 preview

The **M/M/1 queue** is a CTMC on $\{0, 1, 2, \dots\}$ with arrival rate $\lambda$ and service rate $\mu$:

- From state $k$: arrivals at rate $\lambda$ (jump to $k+1$), services at rate $\mu$ (jump to $k-1$ if $k \ge 1$).
- From state 0: only arrivals.

Generator rows:

$$
Q_{k, k+1} = \lambda,\quad Q_{k, k-1} = \mu \ (k \ge 1),\quad Q_{kk} = -(\lambda + \mu).
$$

**Detailed balance** (which holds for any birth–death chain):

$$
\pi_k \lambda = \pi_{k+1} \mu \;\Longrightarrow\; \pi_k = \pi_0\, \rho^k,\quad \rho = \frac{\lambda}{\mu}.
$$

Normalizing: if $\rho < 1$,

$$
\boxed{\;\pi_k = (1-\rho)\,\rho^k,\quad k = 0, 1, 2, \dots\;}
$$

the stationary distribution of an M/M/1 queue is **geometric** with parameter $1-\rho$. If $\rho \ge 1$ the queue is **not positive recurrent** — the stationary distribution does not exist (queue grows unboundedly).

**Mean queue length.**

$$
\mathbb E[K] = \sum_{k=0}^\infty k (1-\rho)\rho^k = \frac{\rho}{1-\rho}.
$$

This is the fundamental **M/M/1 waiting-time formula** and is the starting point for all queueing analysis of Chapter 7.

---

## 6.6 Comparison: DTMC vs. CTMC

| Feature | DTMC | CTMC |
|---|---|---|
| Time parameter | $n \in \mathbb Z_{\ge 0}$ | $t \in \mathbb R_{\ge 0}$ |
| Transition rule | $P_{ij}$ (probability) | $Q_{ij}$ (rate) |
| Rows sum | $P$: rows sum to 1 | $Q$: rows sum to 0 |
| Self-transition | $P_{ii} > 0$ allowed | $Q_{ii} = -q_i$ (exit rate) |
| Holding time | One unit (deterministic) | $\text{Exp}(q_i)$ (random) |
| Evolution | $\mu_n = \mu_0 P^n$ | $\mu(t) = \mu_0 e^{Qt}$ |
| Stationary condition | $\pi P = \pi$ | $\pi Q = 0$ |
| Ergodic theorem | $P^n \to \mathbf 1\pi$ if aperiodic | $P(t) \to \mathbf 1\pi$ always (if ergodic) |

Note: **CTMCs have no periodicity** — the exponential holding times smooth out any parity effects, so aperiodicity is automatic.

---

## 6.7 Computational implementation

### Python — DTMC, classification, CTMC

```python
import numpy as np
from numpy.linalg import eig

# ---- (i) two-state weather chain ----------------------------------------
P = np.array([[0.9, 0.1], [0.5, 0.5]])

# n-step distribution from a fixed start
mu0 = np.array([1.0, 0.0])  # start sunny
for n in [1, 2, 5, 10]:
    mu_n = mu0 @ np.linalg.matrix_power(P, n)
    print(f"n={n:2d}: μ_n = {mu_n}")

# stationary by iterative application
mu = mu0.copy()
for _ in range(200):
    mu = mu @ P
print(f"stationary (iterated): {mu}")

# stationary by eigenvector of Pᵀ
eigvals, eigvecs = eig(P.T)
idx = np.argmin(np.abs(eigvals - 1))
pi = eigvecs[:, idx].real
pi = pi / pi.sum()
print(f"stationary (eigen):    {pi}")

# ---- (ii) solve π P = π by linear system -------------------------------
n = P.shape[0]
A = P.T - np.eye(n)
A[-1, :] = np.ones(n)                     # replace last row with normalization
b = np.zeros(n); b[-1] = 1
pi_ls = np.linalg.solve(A, b)
print(f"stationary (linear):   {pi_ls}")

# ---- (iii) classification: irreducibility check --------------------------
def communicating_classes(P, tol=1e-12):
    """Return list of communicating classes by reaching relation."""
    n = P.shape[0]
    # reachable[i] = set of j reachable from i (using BFS on adjacency)
    reachable = [set() for _ in range(n)]
    for i in range(n):
        stack = [i]; seen = {i}
        while stack:
            u = stack.pop()
            for v in range(n):
                if P[u, v] > tol and v not in seen:
                    seen.add(v); stack.append(v)
        reachable[i] = seen
    # i ~ j iff j in reachable[i] and i in reachable[j]
    classes = []
    assigned = set()
    for i in range(n):
        if i in assigned:
            continue
        cls = {j for j in range(n)
               if j in reachable[i] and i in reachable[j]}
        classes.append(sorted(cls))
        assigned |= cls
    return classes

P_irred = np.array([[0.5, 0.5, 0.0],
                    [0.3, 0.4, 0.3],
                    [0.0, 0.5, 0.5]])
print(f"\nclasses: {communicating_classes(P_irred)}")

P_reduc = np.array([[1.0, 0.0, 0.0],
                    [0.5, 0.0, 0.5],
                    [0.0, 0.0, 1.0]])
print(f"classes: {communicating_classes(P_reduc)}")   # absorbing 0 and 2; {1}

# ---- (iv) CTMC: two-state up/down --------------------------------------
alpha, beta = 0.1, 1.0
Q = np.array([[-alpha, alpha], [beta, -beta]])
print(f"\nstationary CTMC by π Q = 0:")

# replace last row of Q with ones
A = Q.copy()
A[-1, :] = 1.0
b = np.array([0.0, 1.0])
pi_ct = np.linalg.solve(A, b)
print(f"  π = {pi_ct}  (theory: β/(α+β), α/(α+β))")

# time-dependent P(t) via matrix exponential (scipy)
from scipy.linalg import expm
t = 1.0
Pt = expm(Q * t)
print(f"P(1) =\n{Pt}")

# ---- (v) M/M/1 stationary distribution ---------------------------------
lam, mu = 0.7, 1.0
rho = lam / mu
pi_MM1 = np.array([(1 - rho) * rho**k for k in range(20)])
print(f"\nM/M/1: ρ = {rho}, E[K] = {rho/(1-rho):.4f}")
print(f"  π[0:5] = {pi_MM1[:5]}")
```

### Java — iterative stationary distribution + classification

```java
import java.util.*;

public class MarkovChains {

    // Power iteration to compute stationary distribution: π_{k+1} = π_k P
    static double[] stationaryIter(double[][] P, int maxIter, double tol) {
        int n = P.length;
        double[] pi = new double[n];
        pi[0] = 1.0;                                // start at state 0
        for (int k = 0; k < maxIter; k++) {
            double[] next = new double[n];
            for (int j = 0; j < n; j++)
                for (int i = 0; i < n; i++)
                    next[j] += pi[i] * P[i][j];
            double diff = 0;
            for (int j = 0; j < n; j++)
                diff = Math.max(diff, Math.abs(next[j] - pi[j]));
            pi = next;
            if (diff < tol) break;
        }
        return pi;
    }

    // Irreducibility: single communicating class iff BFS from any i reaches all
    static boolean irreducible(double[][] P, double tol) {
        int n = P.length;
        for (int start = 0; start < n; start++) {
            boolean[] seen = new boolean[n];
            ArrayDeque<Integer> stack = new ArrayDeque<>();
            stack.push(start); seen[start] = true;
            while (!stack.isEmpty()) {
                int u = stack.pop();
                for (int v = 0; v < n; v++)
                    if (P[u][v] > tol && !seen[v]) {
                        seen[v] = true; stack.push(v);
                    }
            }
            for (int j = 0; j < n; j++) if (!seen[j]) return false;
        }
        return true;
    }

    public static void main(String[] args) {
        double[][] P = { {0.9, 0.1}, {0.5, 0.5} };
        double[] pi = stationaryIter(P, 1000, 1e-12);
        System.out.printf("weather stationary: [%.6f, %.6f]%n", pi[0], pi[1]);
        System.out.println("irreducible: " + irreducible(P, 1e-12));

        double[][] P3 = { {0.5, 0.5, 0.0},
                          {0.3, 0.4, 0.3},
                          {0.0, 0.5, 0.5} };
        double[] pi3 = stationaryIter(P3, 1000, 1e-12);
        System.out.printf("3-state stationary: [%.6f, %.6f, %.6f]%n",
                          pi3[0], pi3[1], pi3[2]);
        System.out.println("irreducible: " + irreducible(P3, 1e-12));

        double[][] P_absorbing = { {1.0, 0.0, 0.0},
                                    {0.5, 0.0, 0.5},
                                    {0.0, 0.0, 1.0} };
        System.out.println("absorbing irreducible: "
                            + irreducible(P_absorbing, 1e-12));
    }
}
```

### C# — DTMC simulation + CTMC stationary by πQ = 0

```csharp
using System;
using System.Linq;

class MarkovChains
{
    // Simulate DTMC: draw n steps, return empirical frequencies
    static double[] Simulate(double[,] P, int[] initialState, int nSteps, int seed)
    {
        int n = P.GetLength(0);
        var rng = new Random(seed);
        var counts = new int[n];
        int state = initialState[0];
        for (int step = 0; step < nSteps; step++)
        {
            double u = rng.NextDouble(), cum = 0.0;
            for (int j = 0; j < n; j++)
            {
                cum += P[state, j];
                if (u <= cum) { state = j; break; }
            }
            counts[state]++;
        }
        double total = counts.Sum();
        return counts.Select(c => c / total).ToArray();
    }

    // Solve πQ = 0 with normalization by Gaussian elimination
    static double[] StationaryCTMC(double[,] Q)
    {
        int n = Q.GetLength(0);
        var A = (double[,])Q.Clone();
        // replace last row with ones
        for (int j = 0; j < n; j++) A[n - 1, j] = 1.0;
        var b = new double[n]; b[n - 1] = 1.0;

        // Gauss-Jordan
        for (int col = 0; col < n; col++)
        {
            int piv = col;
            for (int r = col + 1; r < n; r++)
                if (Math.Abs(A[r, col]) > Math.Abs(A[piv, col])) piv = r;
            for (int j = col; j < n; j++)
                (A[col, j], A[piv, j]) = (A[piv, j], A[col, j]);
            (b[col], b[piv]) = (b[piv], b[col]);

            double d = A[col, col];
            for (int j = col; j < n; j++) A[col, j] /= d;
            b[col] /= d;

            for (int r = 0; r < n; r++)
            {
                if (r == col) continue;
                double factor = A[r, col];
                for (int j = col; j < n; j++) A[r, j] -= factor * A[col, j];
                b[r] -= factor * b[col];
            }
        }
        return b;
    }

    static void Main()
    {
        // Simulate weather chain
        var P = new double[,] { {0.9, 0.1}, {0.5, 0.5} };
        var emp = Simulate(P, new[]{0}, 1_000_000, 42);
        Console.WriteLine($"simulated freq: [{emp[0]:F4}, {emp[1]:F4}]  (theory: [0.8333, 0.1667])");

        // Two-state CTMC
        double alpha = 0.1, beta = 1.0;
        var Q = new double[,] { {-alpha, alpha}, {beta, -beta} };
        var pi = StationaryCTMC(Q);
        Console.WriteLine($"CTMC stationary: [{pi[0]:F4}, {pi[1]:F4}]  " +
                          $"(theory: [{beta/(alpha+beta):F4}, {alpha/(alpha+beta):F4}])");
    }
}
```

!!! tip "Tooling"
    - **Python:** `numpy.linalg` for linear systems and matrix powers, `scipy.linalg.expm`
      for CTMC transition matrices, `numpy.random` for simulation. `markovchain` and
      `PyDTMC` provide higher-level APIs.
    - **Java:** no special library needed for small problems; Apache Commons Math
      provides `RealMatrix` + `DecompositionSolver` for larger ones.
    - **C#:** Math.NET Numerics provides `Matrix<double>.Solve` and `Matrix.Exp`.
    - **MATLAB:** `dtmc`, `mc` (econometrics toolbox), `ctmc` for building chains; `simulate`.

---

## 6.8 Chapter summary

1. **Markov property:** the future depends on the past only through the present. This is what makes Markov chains tractable.
2. **Transition matrix $P$** encodes the entire DTMC. Rows sum to 1. Multi-step probabilities come from $P^n$ via the Chapman–Kolmogorov equation.
3. **Classification** — irreducible, aperiodic, positive recurrent — determines the long-run behaviour: an irreducible, positive-recurrent, aperiodic finite chain converges to a unique stationary distribution.
4. **Stationary distribution $\pi$** solves $\pi P = \pi$. **Detailed balance** $\pi_i P_{ij} = \pi_j P_{ji}$ is a stronger condition (reversibility) that makes $\pi$ easy to find, especially for birth–death chains.
5. **Ergodic theorem:** for ergodic chains, time averages equal ensemble averages with respect to $\pi$ — the basis of MCMC and simulation.
6. **Continuous-time chains** use the **generator $Q$** instead of $P$; the stationary distribution solves $\pi Q = 0$. Holding times are exponentially distributed, and the chain's randomness is in the timing plus the jump chain.
7. **Birth–death CTMCs** give the M/M/1 queue: geometric stationary distribution $\pi_k = (1-\rho)\rho^k$, mean queue length $\rho/(1-\rho)$, the gateway to Chapter 7.

---

## 6.9 Quiz

### A. Multiple choice (choose the single best answer)

1. The **Markov property** says that:
   (a) the process has no memory (b) the future depends on the past only through the present (c) the state space is finite (d) the transition matrix is symmetric

2. For a DTMC with transition matrix $P$, the $n$-step transition probabilities are:
   (a) $nP$ (b) $P^n$ (c) $P/n$ (d) the stationary distribution of $P$

3. A probability vector $\pi$ is **stationary** if:
   (a) $\pi P = \pi$ (b) $P \pi = \pi$ (c) $\pi^\top P = \mathbf 1$ (d) $P = \pi\pi^\top$

4. Which condition is **not** required for a chain to be ergodic?
   (a) irreducible (b) positive recurrent (c) aperiodic (d) symmetric transition matrix

5. In an **irreducible, positive-recurrent, aperiodic** finite chain, $\lim_{n\to\infty} P^n$ equals:
   (a) $I$ (b) $0$ (c) $\mathbf 1\pi$ (d) $P$

6. **Detailed balance** holds when:
   (a) $\pi_i P_{ij} = \pi_j P_{ji}$ for all $i,j$ (b) $P$ is symmetric (c) the chain is periodic (d) $P^n$ converges

7. For a two-state CTMC with generator $Q = \begin{bmatrix} -\alpha & \alpha \\ \beta & -\beta \end{bmatrix}$, the stationary distribution is:
   (a) $(\alpha/(\alpha+\beta),\ \beta/(\alpha+\beta))$ (b) $(\beta/(\alpha+\beta),\ \alpha/(\alpha+\beta))$ (c) $(1/2, 1/2)$ (d) $(1, 0)$

8. In a CTMC, the **holding time** in state $i$ is:
   (a) deterministic (b) uniform (c) exponentially distributed with rate $-Q_{ii}$ (d) geometric

9. The Kolmogorov **forward** equation is:
   (a) $\dot P(t) = Q P(t)$ (b) $\dot P(t) = P(t) Q$ (c) $\dot P(t) = P(t) + Q$ (d) $P(t) = Q t$

10. For an M/M/1 queue with arrival rate $\lambda$ and service rate $\mu$, the stationary distribution of queue length is:
    (a) uniform (b) geometric with parameter $1-\rho$ where $\rho = \lambda/\mu$ (c) Poisson with mean $\rho$ (d) normal

11. If $\rho \ge 1$ in an M/M/1 queue, then:
    (a) the queue is positive recurrent (b) the stationary distribution does not exist and the queue grows unboundedly (c) the stationary distribution is uniform (d) the queue empties instantly

12. Which statement about CTMCs is **true**?
    (a) they always have a period (b) $Q$ is a stochastic matrix (c) exponential holding times make aperiodicity automatic (d) stationary distribution is $\pi = Q$

### B. Fill in the blank

1. The Markov property: $\Pr\{X_{n+1}=j \mid X_0,\dots,X_n\} = \Pr\{X_{n+1}=j \mid \_\_\_\_\}$.
2. Each row of a DTMC transition matrix $P$ sums to \_\_\_\_.
3. Each row of a CTMC generator $Q$ sums to \_\_\_\_.
4. The Chapman–Kolmogorov equation says $P^{m+n} = \_\_\_\_$.
5. A state $i$ is **recurrent** if the probability of returning to $i$ eventually is \_\_\_\_.
6. The stationary distribution of a DTMC satisfies $\pi P = \_\_\_\_$ and $\sum \pi_i = \_\_\_\_$.
7. The **period** of an aperiodic state is \_\_\_\_.
8. Detailed balance: $\pi_i P_{ij} = \_\_\_\_$.
9. The CTMC stationary distribution solves $\pi \_\_\_\_ = 0$.
10. In an M/M/1 queue, the mean number in the system is \_\_\_\_.

### C. Define / explain

1. State the **Markov property** for a discrete-time process. Give one example of a process that satisfies it and one that does *not*, and explain why.
2. State the **Chapman–Kolmogorov equation** and prove it using the law of total probability and the Markov property.
3. Define **recurrence**, **transience**, and **positive recurrence**. Give an example of each (or explain why an example of the third is hard to construct for finite state spaces).
4. State and prove the **ergodic theorem** for finite irreducible aperiodic DTMCs. What changes if the chain is periodic?
5. State **detailed balance**. Explain why it is a *sufficient* condition for stationarity, but not a necessary one, and give an example of a chain that has a stationary distribution but does not satisfy detailed balance.
6. Write down the **generator matrix** of a birth–death CTMC with arrival rates $\lambda_k$ and departure rates $\mu_k$. Derive the stationary distribution via detailed balance, and state the positive-recurrence condition.
7. Derive the **M/M/1 stationary distribution** and the mean queue length. What happens when $\rho = \lambda/\mu \ge 1$?
8. Compare **DTMCs and CTMCs** in three respects: rows of the transition/generator matrix, holding-time distribution, and existence of periodicity.

---

## 6.10 Answer key

??? success "A. Multiple choice"
    1. **(b)** — future depends on past only through present.
    2. **(b)** — $P^n$ gives $n$-step probabilities.
    3. **(a)** — $\pi P = \pi$ and $\sum \pi_i = 1$.
    4. **(d)** — symmetric $P$ is not required for ergodicity (and symmetric $P$ would make $\pi$ uniform, but that is a consequence, not a requirement).
    5. **(c)** — $\mathbf 1\pi$, every row equal to $\pi$.
    6. **(a)** — $\pi_i P_{ij} = \pi_j P_{ji}$.
    7. **(b)** — $(\beta/(\alpha+\beta),\ \alpha/(\alpha+\beta))$.
    8. **(c)** — exponential with rate $-Q_{ii}$.
    9. **(b)** — $\dot P = PQ$ (forward); $\dot P = QP$ is backward.
    10. **(b)** — geometric with parameter $1-\rho$.
    11. **(b)** — not positive recurrent; unbounded growth.
    12. **(c)** — exponential holding times ⇒ aperiodicity automatic.

??? success "B. Fill in the blank"
    1. $X_n$
    2. $1$
    3. $0$
    4. $P^m P^n$
    5. $1$ (returns almost surely)
    6. $\pi$; $1$
    7. $1$
    8. $\pi_j P_{ji}$
    9. $Q$
    10. $\rho/(1-\rho)$

??? success "C. Define / explain — model answers"
    1. *Markov property:* for every $n$ and states $i_0,\dots,i_n,j$, $\Pr\{X_{n+1}=j \mid X_0=i_0,\dots,X_n=i_n\} = \Pr\{X_{n+1}=j \mid X_n=i_n\}$ — the future depends on the past only through the present. *Example that satisfies it:* a queue-length process at packet-arrival epochs (state = current queue length, transition = arrival/service). *Example that does not:* the position of a moving object with inertia (velocity is not in the state, so knowing current position is insufficient to predict the next position — you also need the recent past). Adding velocity to the state restores the Markov property. Another non-example: an autoregressive process $X_{n+1} = \sum_{k=1}^{p} \phi_k X_{n-k+1} + W_n$ is not Markov in $X_n$ alone, but the vector $(X_n, X_{n-1}, \dots, X_{n-p+1})$ is.
    2. *CK equation:* $P^{m+n} = P^m P^n$, i.e. $\Pr\{X_{m+n}=j \mid X_0=i\} = \sum_k \Pr\{X_{m+n}=j \mid X_m=k\}\Pr\{X_m=k \mid X_0=i\}$. *Proof:* by the law of total probability, $\Pr\{X_{m+n}=j \mid X_0=i\} = \sum_k \Pr\{X_{m+n}=j \mid X_m=k, X_0=i\}\Pr\{X_m=k \mid X_0=i\}$. By the Markov property (conditioned on $X_m=k$, the future does not depend on $X_0$), $\Pr\{X_{m+n}=j \mid X_m=k, X_0=i\} = \Pr\{X_{m+n}=j \mid X_m=k\} = (P^n)_{kj}$. Hence the equation.
    3. *Recurrent:* $f_i = 1$; the chain returns to $i$ almost surely. *Transient:* $f_i < 1$; positive probability of never returning. *Positive recurrent:* recurrent **and** finite mean return time $m_i < \infty$. *Examples:* (a) **Recurrent:** simple symmetric random walk on $\mathbb{Z}$; two-state alternating chain; any finite irreducible chain. (b) **Transient:** simple symmetric random walk on $\mathbb{Z}^3$; biased random walk on $\mathbb{Z}$ with drift; gambler's ruin with state 0 absorbing (states $1,2,3$ are transient). (c) **Null recurrent (positive-recurrence failure):** simple symmetric random walk on $\mathbb{Z}$ is recurrent but **null recurrent** — $m_i = \infty$. In **finite** state spaces every recurrent state is positive recurrent, so null recurrence requires an infinite state space; this is why example (c) is impossible for a finite chain.
    4. *Ergodic theorem:* for a finite irreducible aperiodic DTMC with stationary $\pi$, $P^n \to \mathbf 1\pi$ and $\frac{1}{N}\sum_{n=1}^N f(X_n) \to \sum \pi_j f(j)$ a.s. for bounded $f$. *Proof sketch:* Irreducibility gives a unique stationary distribution. Aperiodicity ensures that $P$ has only one eigenvalue of modulus 1 (the Perron eigenvalue 1), and all other eigenvalues have magnitude $< 1$. Jordan decomposition of $P = \mathbf 1\pi + R$ with $\rho(R) < 1$ gives $P^n = \mathbf 1\pi + R^n \to \mathbf 1\pi$. The time-average convergence follows from the Birkhoff ergodic theorem or by the same eigenvalue argument on $P$ applied to $f$. *If periodic:* $P^n$ does not converge, but the Cesàro average $\frac{1}{N}\sum P^k \to \mathbf 1\pi$ still does — periodicity affects pointwise-in-$n$ limits, not Cesàro limits.
    5. *Detailed balance:* $\pi_i P_{ij} = \pi_j P_{ji}$ for all $i, j$ — the chain looks the same forward and backward in time when started from $\pi$. *Sufficiency:* summing over $i$ gives $\sum_i \pi_i P_{ij} = \pi_j \sum_i P_{ji} = \pi_j$ (since rows sum to 1), i.e. $\pi P = \pi$. *Non-necessity example:* Consider the three-state directed cycle $P_{12}=P_{23}=P_{31}=1$, all other entries 0. This chain has stationary $\pi = (1/3, 1/3, 1/3)$ but $\pi_1 P_{12} = 1/3 \ne \pi_2 P_{21} = 0$, so detailed balance fails. The chain is not time-reversible — every step goes $1 \to 2 \to 3 \to 1$, and reversing time gives $1 \to 3 \to 2 \to 1$, a different chain.
    6. *Birth–death generator:* $Q_{k,k+1} = \lambda_k$, $Q_{k,k-1} = \mu_k$ (for $k \ge 1$), $Q_{kk} = -(\lambda_k + \mu_k)$. *Detailed balance:* $\pi_k \lambda_k = \pi_{k+1} \mu_{k+1}$, so $\pi_k = \pi_0 \prod_{j=0}^{k-1} \lambda_j/\mu_{j+1}$. *Normalization:* $\pi_0 = \left(1 + \sum_{k=1}^\infty \prod_{j=0}^{k-1} \lambda_j/\mu_{j+1}\right)^{-1}$. *Positive recurrence condition:* the infinite sum converges. If it diverges, the chain is transient or null recurrent.
    7. *M/M/1 derivation:* state $k = $ number in system; births at rate $\lambda$, deaths at rate $\mu$ (only from $k \ge 1$). Detailed balance gives $\pi_k \lambda = \pi_{k+1} \mu$, so $\pi_k = \pi_0 \rho^k$ with $\rho = \lambda/\mu$. Normalizing $\sum_k \pi_k = 1$ requires $\rho < 1$, giving $\pi_0 = 1 - \rho$ and $\pi_k = (1-\rho)\rho^k$. *Mean queue length:* $E[K] = \rho/(1-\rho)$. *If $\rho \ge 1$:* the normalization sum diverges — no stationary distribution exists. The queue grows unboundedly at average rate $\lambda - \mu > 0$, and its expected length is infinite.
    8. *Comparison:* (i) *Row sums:* $P$ rows sum to 1 (probability); $Q$ rows sum to 0 (rates). (ii) *Holding times:* deterministic unit time in DTMC; $\text{Exp}(-Q_{ii})$ in CTMC. (iii) *Periodicity:* DTMCs can be periodic (e.g., two-state alternating chain has period 2); CTMCs are always aperiodic because exponential holding times smear out any parity structure.

---

## Further reading

- Norris, *Markov Chains*, Cambridge University Press, 1997 — the standard reference on discrete- and continuous-time chains.
- Karlin & Taylor, *A First Course in Stochastic Processes*, 2nd ed., Academic Press — the classical textbook with excellent coverage of classification and limiting behaviour.
- Ross, *Introduction to Probability Models*, 11th ed., Academic Press — Chapters 4–6 cover DTMC, exponential distributions, and CTMC.
- Kemeny & Snell, *Finite Markov Chains*, Springer — the definitive reference for finite-state chains and matrix methods.
- Levin, Peres & Wilmer, *Markov Chains and Mixing Times*, AMS, 2009 — modern treatment with applications to MCMC and mixing bounds.
- Stewart, *Introduction to the Numerical Solution of Markov Chains*, Princeton University Press, 1994 — the reference for computational methods (power iteration, GTH algorithm, block methods).
- Häggström, *Finite Markov Chains and Algorithmic Applications*, Cambridge University Press, 2002 — concise, with MCMC applications.
- Robert & Casella, *Monte Carlo Statistical Methods*, Springer — Markov chain theory in the service of simulation.
- Python `markovchain`, `PyDTMC`, `networkx` (for graph/transition diagrams), `scipy.linalg.expm`.

---

## Suggested tutorial problems for this chapter

1. **Two-state chain.** For the weather chain of Example 6.1, compute $P^n$ in closed form using the eigen-decomposition of $P$. Verify that $P^n \to \mathbf 1\pi$ and that the convergence rate is governed by the subdominant eigenvalue $0.4$. What is the mixing time to reach $\|P^n - \mathbf 1\pi\|_\infty < 10^{-3}$?
2. **Gambler's ruin.** For the chain of Example 6.2 with $p = 0.45$, compute the absorption probabilities and the expected absorption time for $k = 1, 2, 3$. Verify against a simulation of 10,000 independent runs.
3. **Ehrenfest urn.** Simulate the Ehrenfest chain for $N = 20$ and confirm that the empirical frequencies of $k$ converge to the binomial $\binom{20}{k}/2^{20}$, **averaged over even and odd steps separately**. Why is the parity effect visible in the raw empirical histogram?
4. **Classification.** For the transition matrices
   $P_1 = \begin{bmatrix}1&0\\0&1\end{bmatrix}$,
   $P_2 = \begin{bmatrix}0&1\\1&0\end{bmatrix}$,
   $P_3 = \begin{bmatrix}0.5&0.5\\0.5&0.5\end{bmatrix}$,
   determine for each: communicating classes, recurrence/transience, period, stationary distribution, and whether $P^n$ converges.
5. **Positive recurrence from scratch.** Show that a simple symmetric random walk on $\mathbb{Z}$ is null recurrent (not positive recurrent) — i.e., $f_0 = 1$ but $E[T_0] = \infty$.
6. **Birth–death CTMC.** For the birth–death chain with $\lambda_k = \lambda$ for all $k$ and $\mu_k = \mu$ for all $k \ge 1$ (M/M/1), derive $\pi$. Then take $\lambda_k = \lambda/(k+1)$ and $\mu_k = \mu$ (an M/M/∞ queue — unlimited servers) and derive the stationary distribution. Verify it is Poisson.
7. **Detailed balance check.** Show that the two-state chain of Example 6.4 satisfies detailed balance with respect to its stationary distribution. Then construct a 3-state chain with stationary $\pi = (1/3, 1/3, 1/3)$ that does *not* satisfy detailed balance (a directed cycle).
8. **M/M/1 simulation.** Simulate an M/M/1 queue with $\lambda = 0.8$, $\mu = 1.0$ for $T = 10^5$ time units. Plot the queue-length trajectory and the empirical distribution. Compare with $(1-\rho)\rho^k$. Repeat with $\lambda = 1.1$ and observe the divergence.
9. **PageRank preview.** The Google PageRank algorithm models web navigation as a DTMC on the graph of pages, with transition $P_{ij} = (1-d)/N + d\, a_{ij}/L_i$ for damping $d \in (0,1)$. Show that the PageRank vector is the stationary distribution of this chain. For a 4-page graph, compute it numerically.
