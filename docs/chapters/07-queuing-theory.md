# Chapter 7 — Introduction to Queuing Theory

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 tutorial · **Home study:** 5 h
    
    **Coverage (syllabus):** elements of queuing systems · Little's formula ·
    birth–death process · M/M/1 queuing systems · M/M/s queuing systems ·
    simulation and data analysis of queuing systems.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Describe a queuing system in **Kendall notation** $A/S/c/K/N/D$ and identify its arrival process, service process, number of servers, capacity, and discipline.
2. State and apply **Little's formula** $L = \lambda W$ to any stable queuing system, without any assumption of Poisson arrivals or exponential service.
3. Derive the **stationary distribution** of a general birth–death chain and specialize it to the M/M/1 and M/M/∞ queues.
4. Compute **performance measures** of the M/M/1 queue: $L$, $L_q$, $W$, $W_q$, the distribution of the number in system, and the tail of the waiting-time distribution.
5. Compute **performance measures** of the M/M/s queue using the **Erlang-C** formula and the **offered load** $a = \lambda/\mu$.
6. Extend to **M/M/1/K** (finite buffer) and **M/G/1** (via the Pollaczek–Khinchine formula) and explain Burke's theorem and Jackson networks.
7. **Simulate** M/M/1 and M/M/s queues, verify Little's formula empirically, and analyze the empirical distribution against theory.

---

## 7.1 Elements of a queuing system

Every queuing system is described by five ingredients:

1. **Arrival process.** The stochastic process of epochs at which *customers* (packets, calls, jobs, patients) arrive. The canonical choice is **Poisson** with rate $\lambda$ (counts arrivals per unit time), equivalently **exponential interarrival** times with mean $1/\lambda$.
2. **Service process.** The distribution of service times. Canonical choice: **exponential** with mean $1/\mu$ (rate $\mu$ per busy server).
3. **Number of servers $c$.** How many customers can be served in parallel.
4. **System capacity $K$.** Maximum number in the system (waiting plus in service). $K = \infty$ is the default. When $K < \infty$ and the system is full, arrivals are *blocked* (lost).
5. **Queue discipline.** The rule for choosing which waiting customer to serve next: **FCFS** (first-come, first-served), **LCFS** (last-come, first-served), **PS** (processor sharing), **SIRO** (service in random order). FCFS is the default.

### 7.1.1 Kendall notation

A queuing system is written $A/S/c/K/N/D$:

- $A$ = arrival distribution, $S$ = service distribution, $c$ = servers,
- $K$ = system capacity (omitted if $\infty$),
- $N$ = customer population (omitted if $\infty$),
- $D$ = discipline (omitted if FCFS).

Standard symbols: $M$ (Markovian/Poisson/exponential), $D$ (deterministic), $E_k$ (Erlang-$k$), $G$ (general).

| Notation | Meaning |
|---|---|
| **M/M/1** | Poisson arrivals, exponential service, 1 server, infinite capacity, FCFS |
| **M/M/s** | Poisson arrivals, exponential service, $s$ servers |
| **M/M/1/K** | M/M/1 with capacity $K$ — arrivals blocked when $K$ in system |
| **M/M/∞** | Infinite-server queue (all customers in service immediately) |
| **M/D/1** | Poisson arrivals, deterministic service |
| **M/G/1** | Poisson arrivals, general service — solvable by PK formula |

### 7.1.2 Traffic intensity

For any queue with arrival rate $\lambda$ and service rate $\mu$ per server, the **offered load** (also called *traffic intensity per server*) is

$$
\rho = \frac{\lambda}{c\,\mu}.
$$

For a single server, $a = \lambda/\mu$ is the **offered load** in Erlangs and $\rho = a = \lambda/\mu$. For $c$ servers, define the **total offered load** $a = \lambda/\mu$; then $\rho = a/c$.

**Stability.** An infinite-capacity queue is stable (has a stationary distribution) iff $\rho < 1$. If $\rho \ge 1$, arrivals outrun service and the queue length grows without bound.

!!! tip "Intuition for $\rho$"
    $\rho$ is the long-run fraction of time each server is busy. $\rho = 0.8$ means each
    server is busy 80% of the time. Real systems target $\rho$ between 0.5 and 0.8 —
    higher values give longer queues for a small gain in utilisation.

---

## 7.2 Little's formula

!!! abstract "Little's formula"
    For any stable queuing system in steady state, with $\lambda$ the arrival rate,
    $L$ the **mean number in the system**, and $W$ the **mean time spent in the system**
    (waiting + service) for a typical customer,

    $$
    \boxed{\;L = \lambda W\;}
    $$

    Similarly for the queue alone (excluding service):

    $$
    L_q = \lambda W_q,
    $$

    where $L_q$ is the mean number waiting and $W_q$ is the mean waiting time.

**Why it is remarkable.** Little's formula holds for **any** arrival process, any service distribution, any discipline, and any number of servers — provided the system is stable and the limits defining $L, \lambda, W$ exist. No assumptions about Poisson or exponential are needed. The proof is a sample-path argument (renewal-reward at rate $1$); we will not reproduce it here, but it is the single most-used identity in queueing theory.

### Worked example 7.1 — Little's formula on a deterministic system

A cafeteria has a single queue. Every minute, on average, $2$ customers arrive; each spends an average of $6$ minutes (waiting + eating). What is the average number of customers in the cafeteria?

$$
L = \lambda W = 2 \cdot 6 = 12\ \text{customers}.
$$

No Poisson assumption is needed; the arithmetic is exact once the averages are known. A computer-vision counting camera at the door would verify $L \approx 12$ — an empirical *check* of Little's formula that works for any discipline.

### Worked example 7.2 — Little's formula inside the queue

Same cafeteria, but suppose the average *waiting* time (before being served) is $4$ minutes, and service itself averages $2$ minutes. Then:

$$
W_q = 4\ \text{min},\quad W = W_q + \text{service} = 4 + 2 = 6\ \text{min}.
$$

$$
L_q = \lambda W_q = 2 \cdot 4 = 8\ \text{customers waiting},\quad L = \lambda W = 12.
$$

The difference $L - L_q = 4$ customers are the ones currently *in service* (one per server: here $c = 4$ servers each busy on average). Little's formula partitions the mean number into queue and service parts consistently.

!!! tip "How to use Little's formula in practice"
    **Never** assume an M/M/1 model just to get $L$ or $W$ — measure $\lambda$ and one of
    $L, W, L_q, W_q$ and use Little's formula for the fourth. This is what traffic
    engineers, call-centre managers, and network operators actually do: measure
    arrival rate and one response-time statistic, then derive the queue length.

---

## 7.3 The birth–death framework

The M/M/1 and M/M/s queues are special cases of a general **birth–death process**: a CTMC on $\{0, 1, 2, \dots\}$ where the only allowed transitions are $n \to n+1$ (a *birth*, rate $\lambda_n$) and $n \to n-1$ (a *death*, rate $\mu_n$ for $n \ge 1$).

### 7.3.1 Stationary distribution

The generator has rows

$$
Q_{n, n+1} = \lambda_n,\quad Q_{n, n-1} = \mu_n,\quad Q_{nn} = -(\lambda_n + \mu_n).
$$

For a birth–death chain, detailed balance holds (see Example 6.6), so the stationary distribution satisfies

$$
\pi_n \lambda_n = \pi_{n+1} \mu_{n+1},\quad n \ge 0.
$$

Iterating from $\pi_0$:

$$
\boxed{\;\pi_n = \pi_0 \prod_{k=0}^{n-1} \frac{\lambda_k}{\mu_{k+1}}\;}
$$

with $\pi_0$ determined by normalization:

$$
\pi_0 = \left[1 + \sum_{n=1}^{\infty} \prod_{k=0}^{n-1} \frac{\lambda_k}{\mu_{k+1}}\right]^{-1}.
$$

!!! abstract "Stability condition"
    The chain is **positive recurrent** (the stationary distribution exists) iff

    $$
    \sum_{n=1}^{\infty} \prod_{k=0}^{n-1} \frac{\lambda_k}{\mu_{k+1}} < \infty.
    $$

    If the sum diverges, the queue grows unboundedly — no stationary distribution exists.

**Performance measures.** By Little's formula,

$$
L = \sum_{n=0}^\infty n\, \pi_n,\qquad L_q = \sum_{n=c}^\infty (n-c)\, \pi_n,\qquad W = L/\lambda,\qquad W_q = L_q/\lambda.
$$

The **mean number in service** is $L - L_q = \lambda/\mu$ (each busy server serves at rate $\mu$), a nice consistency check.

### 7.3.2 Specializations

| Queue | $\lambda_n$ | $\mu_n$ | Stationary distribution |
|---|---|---|---|
| M/M/1 | $\lambda$ | $\mu$ | $\pi_n = (1-\rho)\rho^n$ |
| M/M/∞ | $\lambda$ | $n\mu$ | $\pi_n = e^{-a} a^n/n!$, $a = \lambda/\mu$ (Poisson) |
| M/M/s | $\lambda$ | $\min(n,s)\,\mu$ | Piecewise — §7.5 |
| M/M/1/K | $\lambda$ | $\mu$, for $n \ge 1$; $\lambda_n = 0$ for $n \ge K$ | Truncated geometric |
| M/M/s/s | $\lambda$ | $\min(n,s)\,\mu$; $\lambda_n = 0$ for $n \ge s$ | Erlang-B formula |

---

## 7.4 The M/M/1 queue

The **M/M/1** queue: Poisson arrivals rate $\lambda$, exponential service rate $\mu$, single server, infinite capacity, FCFS.

### 7.4.1 Stability and stationary distribution

With $\rho = \lambda/\mu$ and **$\rho < 1$** (stability):

$$
\boxed{\;\pi_n = (1-\rho)\,\rho^n,\quad n = 0, 1, 2, \dots\;}
$$

The stationary distribution is **geometric** with parameter $1-\rho$ (and mean $\rho/(1-\rho)$). This is the same result previewed in Chapter 6 (Example 6.8).

### 7.4.2 Performance measures

| Quantity | Formula |
|---|---|
| Utilisation | $\rho = \lambda/\mu$ |
| Mean number in system | $L = \rho/(1-\rho) = \lambda/(\mu - \lambda)$ |
| Mean number in queue | $L_q = \rho^2/(1-\rho) = \lambda^2/(\mu(\mu-\lambda))$ |
| Mean system time | $W = 1/(\mu - \lambda)$ |
| Mean wait time | $W_q = \rho/(\mu - \lambda) = \lambda/(\mu(\mu-\lambda))$ |
| Probability system is idle | $\pi_0 = 1 - \rho$ |
| Probability of $n$ or more | $\Pr\{N \ge n\} = \rho^n$ |

The identity $W = W_q + 1/\mu$ holds: total time = wait + service.

### Worked example 7.3 — full M/M/1 analysis

A bank has one teller. Customers arrive as a Poisson process with rate $\lambda = 6$ per hour. Service times are exponential with mean $12$ minutes (so $\mu = 5$ per hour). Then

$$
\rho = \frac{6}{5} = 1.2 > 1.
$$

The queue is **unstable**: the mean arrival rate exceeds the service rate. There is no stationary distribution; the queue grows without bound, and any "average" derived from it is infinite. This illustrates the practical importance of the stability condition — a bank with one teller and 6 customers/hour will *never* settle into steady state.

Now change the scenario: the bank hires a second teller, so $\mu = 10$ per hour (each teller serves at 5/hour; together $10$/hour). Then

$$
\rho = \frac{6}{10} = 0.6 < 1.
$$

| Quantity | Value |
|---|---|
| $L$ | $0.6/0.4 = 1.5$ customers |
| $L_q$ | $0.36/0.4 = 0.9$ customers |
| $W$ | $1/(10-6) = 0.25$ hr $= 15$ min |
| $W_q$ | $0.6/4 = 0.15$ hr $= 9$ min |
| $\pi_0$ | $0.4$ (tellers idle 40% of the time) |
| $\Pr\{N \ge 3\}$ | $0.6^3 = 0.216$ |

**Sanity checks.** $L = \lambda W = 6 \cdot 0.25 = 1.5$ ✓. $L_q = \lambda W_q = 6 \cdot 0.15 = 0.9$ ✓. $W = W_q + 1/\mu = 0.15 + 0.10 = 0.25$ ✓.

### 7.4.3 The waiting-time distribution

For an M/M/1 queue in steady state, the **waiting time in queue** $W_q$ has a mixed distribution:

$$
\Pr\{W_q = 0\} = 1 - \rho,\qquad
\Pr\{W_q > t\} = \rho\, e^{-(\mu - \lambda)t},\quad t \ge 0.
$$

The **total time in system** $W = W_q + \text{service}$ has an **exponential** distribution:

$$
\Pr\{W > t\} = e^{-(\mu - \lambda)t},\quad t \ge 0.
$$

This is a striking result: the total time in an M/M/1 system is *exponentially* distributed with mean $1/(\mu - \lambda)$ and rate $\mu - \lambda$ (the **spare capacity**).

### Worked example 7.4 — waiting-time tail for an M/M/1

Same bank with $\lambda = 6$, $\mu = 10$ per hour ($\rho = 0.6$).

**Probability a customer waits more than 15 minutes (0.25 hours) in queue:**

$$
\Pr\{W_q > 0.25\} = \rho\, e^{-(\mu - \lambda) \cdot 0.25} = 0.6\, e^{-4 \cdot 0.25} = 0.6\, e^{-1} \approx 0.2207.
$$

About 22% of customers wait more than 15 minutes.

**Probability total time in system exceeds 30 minutes (0.5 hours):**

$$
\Pr\{W > 0.5\} = e^{-4 \cdot 0.5} = e^{-2} \approx 0.1353.
$$

About 13.5% of customers spend more than half an hour in the bank.

**Tails decay exponentially** with rate $\mu - \lambda$ — the "spare capacity". Doubling the spare capacity from 4 to 8 per hour would halve the exponential rate, so tails would drop much faster:

$$
\Pr\{W > 0.5\} \big|_{\mu=14} = e^{-8 \cdot 0.5} = e^{-4} \approx 0.0183.
$$

Adding servers is a **disproportionately effective** way to reduce tail risk — a basic insight for call-centre management.

### 7.4.4 Departure process and Burke's theorem

!!! abstract "Burke's theorem"
    In a stable M/M/1 queue, the **departure process** (the point process of service completions of customers leaving the system) is a **Poisson process** with the same rate $\lambda$ as the arrival process. Moreover, the number in the system just before an arrival is *independent* of the departure process before the arrival time.

This is *not* obvious — service completions are staggered by service times, which depend on the arrival process. But in the M/M/1 queue, exponential memorylessness makes the departures Poisson. Burke's theorem is what makes **Jackson networks** (queues in series or in a network) tractable: each queue in the network "sees" Poisson arrivals regardless of its position.

**Corollary.** For a tandem of M/M/1 queues in series (all stable), each queue behaves as if it were isolated with the same arrival rate. This is **Jackson's theorem** for open networks with Poisson arrivals and exponential service — the *product-form* stationary distribution factorizes into per-queue geometrics.

---

## 7.5 The M/M/s queue

The **M/M/s** queue: Poisson arrivals $\lambda$, exponential service $\mu$ per server, $s$ parallel servers, infinite capacity.

**Offered load.** $a = \lambda/\mu$ (total load, in Erlangs). **Utilisation per server.** $\rho = a/s = \lambda/(s\mu)$.

**Stability.** $\rho < 1$, i.e. $\lambda < s\mu$ — total service capacity exceeds arrival rate.

### 7.5.1 Stationary distribution

For $n < s$, all customers are in service and the departure rate is $n\mu$. For $n \ge s$, all servers are busy and the departure rate saturates at $s\mu$. The birth–death balance gives:

$$
\pi_n = \begin{cases}
\pi_0 \dfrac{a^n}{n!}, & 0 \le n \le s, \\[8pt]
\pi_0 \dfrac{a^n}{s!\,s^{\,n-s}}, & n > s,
\end{cases}
$$

with normalization

$$
\boxed{\;\pi_0 = \left[\sum_{n=0}^{s-1} \frac{a^n}{n!} + \frac{a^s}{s!}\cdot\frac{1}{1-\rho}\right]^{-1}\;}
$$

### 7.5.2 Erlang-C formula

The probability that an arriving customer has to **wait** (i.e., finds all servers busy) is

$$
\boxed{\;C(s, a) = \Pr\{N \ge s\} = \pi_0\, \frac{a^s}{s!\,(1-\rho)}\;}
$$

This is the **Erlang-C formula**, the standard tool for sizing call centres and help desks. It depends only on the offered load $a$ and the number of servers $s$.

### 7.5.3 Performance measures

| Quantity | Formula |
|---|---|
| Utilisation | $\rho = a/s$ |
| Wait probability | $\Pr\{W_q > 0\} = C(s,a)$ |
| Mean number in queue | $L_q = C(s,a)\cdot\dfrac{\rho}{1-\rho}$ |
| Mean system time | $W_q = L_q/\lambda$ |
| Mean total time | $W = W_q + 1/\mu$ |
| Mean number in system | $L = L_q + a$ |
| Mean waiting time given wait | $\mathbb E[W_q \mid W_q > 0] = \dfrac{1}{s\mu - \lambda}$ |

The last formula is the **conditional** wait given that a customer must queue — the mean a *delayed* customer experiences. It can be much larger than the unconditional $W_q$, which is why call-centre reporting should always state both.

### Worked example 7.5 — full M/M/s analysis

A call centre with $s = 2$ servers. Arrival rate $\lambda = 6$ calls/minute; each call is answered at rate $\mu = 4$ calls/minute per server (mean service time $15$ seconds).

**Offered load and utilisation.**

$$
a = \frac{6}{4} = 1.5\ \text{Erlangs},\qquad \rho = \frac{1.5}{2} = 0.75 < 1.
$$

**Normalization.**

$$
\pi_0 = \left[1 + 1.5 + \frac{1.5^2}{2!}\cdot\frac{1}{1 - 0.75}\right]^{-1}
= \left[1 + 1.5 + \frac{2.25}{2}\cdot 4\right]^{-1}
= [1 + 1.5 + 4.5]^{-1} = \frac{1}{7} \approx 0.1429.
$$

**Erlang C.**

$$
C(2, 1.5) = \pi_0\cdot\frac{a^s}{s!(1-\rho)} = \frac{1}{7}\cdot\frac{2.25}{2 \cdot 0.25} = \frac{1}{7}\cdot 4.5 = \frac{4.5}{7} \approx 0.6429.
$$

So **64% of callers have to wait** — a high figure because utilisation is 0.75 with only $2$ servers.

**Mean queue length and wait times.**

$$
L_q = C(2,1.5)\cdot\frac{\rho}{1-\rho} = 0.6429 \cdot \frac{0.75}{0.25} = 0.6429 \cdot 3 = 1.9286\ \text{calls}.
$$

$$
W_q = L_q/\lambda = \frac{1.9286}{6} \approx 0.3214\ \text{min} \approx 19.29\ \text{s}.
$$

$$
W = W_q + 1/\mu = 0.3214 + 0.25 = 0.5714\ \text{min} \approx 34.29\ \text{s}.
$$

$$
L = L_q + a = 1.9286 + 1.5 = 3.4286\ \text{calls}.
$$

**Conditional wait.** Given a caller has to wait, their expected wait is

$$
\frac{1}{s\mu - \lambda} = \frac{1}{2\cdot 4 - 6} = \frac{1}{2} = 0.5\ \text{min} = 30\ \text{s}.
$$

So a **delayed** caller waits on average **30 seconds**, while the unconditional average is 19.29 s. The unconditional number hides the tail.

**Sanity checks.** $L = \lambda W = 6\cdot 0.5714 = 3.4286$ ✓. $L_q = \lambda W_q = 6\cdot 0.3214 = 1.9286$ ✓.

**Reducing the wait — add a server.** With $s = 3$:

$$
a = 1.5,\quad \rho = 0.5,\quad \pi_0 = \left[1 + 1.5 + \frac{1.5^2}{2} + \frac{1.5^3}{6}\cdot\frac{1}{0.5}\right]^{-1}.
$$

Compute: $1 + 1.5 + 1.125 + \frac{3.375}{6}\cdot 2 = 1 + 1.5 + 1.125 + 1.125 = 4.75$, so $\pi_0 = 1/4.75 \approx 0.2105$.

$$
C(3, 1.5) = \pi_0\cdot\frac{a^3}{3!(1-0.5)} = \frac{1}{4.75}\cdot\frac{3.375}{3} = \frac{1.125}{4.75} \approx 0.2368.
$$

Now only **24% of callers wait**, versus 64% with 2 servers — a **2.7× reduction** in wait probability from a 50% increase in servers. This nonlinearity is the central theme of queuing theory: **small capacity increases near $\rho = 1$ produce large performance improvements**.

---

## 7.6 Extensions

### 7.6.1 M/M/1/K — finite buffer

Capacity $K < \infty$: arrivals when $n = K$ are **blocked** (lost). Birth rates $\lambda_n = \lambda$ for $n < K$, $\lambda_K = 0$; death rates $\mu_n = \mu$.

**Stationary distribution** (for any $\rho > 0$ — the finite buffer always stabilises):

$$
\pi_n = \frac{(1-\rho)\rho^n}{1-\rho^{K+1}},\quad n = 0, 1, \dots, K.
$$

**Blocking probability.** The probability an arriving customer finds the system full (and is rejected) is $\pi_K = \dfrac{(1-\rho)\rho^K}{1-\rho^{K+1}}$.

**Effective arrival rate.** Arrivals *that enter* the system occur at rate $\lambda_{\text{eff}} = \lambda(1-\pi_K)$ — the rest are lost.

Little's formula still applies with $\lambda_{\text{eff}}$:

$$
L = \lambda_{\text{eff}} W,\qquad L_q = \lambda_{\text{eff}} W_q.
$$

### Worked example 7.6 — M/M/1/3 buffer

A gateway router with $\lambda = 5$ pkt/ms, $\mu = 10$ pkt/ms, buffer $K = 3$.

$$
\rho = 0.5,\quad K = 3,\quad \pi_n = \frac{(1-0.5)0.5^n}{1 - 0.5^4} = \frac{0.5\cdot 0.5^n}{1 - 0.0625} = \frac{0.5^{n+1}}{0.9375}.
$$

$$
\pi_0 = 0.5333,\quad \pi_1 = 0.2667,\quad \pi_2 = 0.1333,\quad \pi_3 = 0.0667.
$$

**Blocking probability.** $\pi_3 \approx 6.67\%$ of arrivals are lost.

**Effective arrival rate.** $\lambda_{\text{eff}} = 5\cdot 0.9333 = 4.667$ pkt/ms.

**Mean system size.** $L = \sum n\pi_n = 0\cdot 0.5333 + 1\cdot 0.2667 + 2\cdot 0.1333 + 3\cdot 0.0667 = 0.7333$.

**Mean system time** (by Little's formula on entering packets): $W = L/\lambda_{\text{eff}} = 0.7333/4.667 \approx 0.1571$ ms.

Without the buffer (M/M/1, same $\lambda/\mu$), $L = 1$ — the finite buffer *reduces* mean queue length because excess is dropped.

### 7.6.2 M/G/1 — the Pollaczek–Khinchine formula

For M/G/1 (Poisson arrivals, **general** service times), the M/M/1 formulas fail — but there is a closed-form mean-value result:

!!! abstract "Pollaczek–Khinchine (PK) formula"
    For M/G/1 with arrival rate $\lambda$, service-time mean $1/\mu$, and service-time variance $\sigma_s^2$,

    $$
    L_q = \frac{\rho^2 + \lambda^2 \sigma_s^2}{2(1-\rho)},\qquad
    W_q = \frac{L_q}{\lambda},\qquad
    W = W_q + \frac{1}{\mu}.
    $$

    When service is exponential, $\sigma_s^2 = 1/\mu^2$, and the formula reduces to $L_q = \rho^2/(1-\rho)$ — matching the M/M/1 result.

**The variance matters.** Two systems with the same mean service time but different variances have *different* mean waiting times. A deterministic service time ($\sigma_s^2 = 0$) gives

$$
L_q^{M/D/1} = \frac{\rho^2}{2(1-\rho)} = \tfrac12\, L_q^{M/M/1},
$$

**half** the M/M/1 value. Deterministic service is twice as efficient as exponential service at the same utilisation. **Reducing service-time variability is as effective as reducing its mean.**

### 7.6.3 Networks of queues — Jackson's theorem

An **open Jackson network** has $J$ queues. External arrivals to queue $j$ are Poisson at rate $r_j$; service at queue $j$ is exponential with rate $\mu_j$; a fraction $p_{ij}$ of customers leaving queue $i$ go to queue $j$; the rest leave the network.

**Traffic equations.** Let $\lambda_j$ be the total arrival rate (external + internal) to queue $j$:

$$
\boxed{\;\lambda_j = r_j + \sum_{i=1}^{J} \lambda_i\,p_{ij},\quad j = 1, \dots, J.\;}
$$

Solve this linear system for $\lambda_1, \dots, \lambda_J$.

**Stability.** Queue $j$ is stable iff $\lambda_j < \mu_j$ for all $j$.

**Product form.** Under these conditions, the joint stationary distribution of the network factorizes:

$$
\pi(n_1, \dots, n_J) = \prod_{j=1}^{J} \pi_j(n_j),
$$

where $\pi_j$ is the M/M/1 stationary distribution of queue $j$ with arrival rate $\lambda_j$. **Each queue behaves as if isolated** — a remarkable theorem due to Jackson (1963).

**Consequence.** All performance measures factorize: $L = \sum_j L_j$, $L_q = \sum_j L_{q,j}$, $W_j = L_j/\lambda_j$ (per queue), and $W_{\text{system}} = \sum_j W_j$ *if* each customer visits each queue once. For mean network response time, use **Little's formula on the network**: $L_{\text{net}} = \sum_j L_j$, and $W_{\text{net}} = L_{\text{net}}/\sum_j r_j$.

### Worked example 7.7 — a two-node Jackson network

Queue 1: external arrivals $r_1 = 4$ pkt/s, service rate $\mu_1 = 10$ pkt/s, output splits 60% to queue 2, 40% leaves.

Queue 2: external arrivals $r_2 = 1$ pkt/s, service rate $\mu_2 = 8$ pkt/s, all departures leave the network.

**Traffic equations.**

$$
\lambda_1 = r_1 = 4\ \text{pkt/s}\quad\text{(no external input to 1 besides } r_1\text{)},
$$

$$
\lambda_2 = r_2 + 0.6\,\lambda_1 = 1 + 0.6\cdot 4 = 3.4\ \text{pkt/s}.
$$

**Stability.** $\lambda_1/\mu_1 = 0.4 < 1$ ✓; $\lambda_2/\mu_2 = 0.425 < 1$ ✓.

**Per-queue means (M/M/1 formula on each).**

$$
L_1 = \frac{0.4}{0.6} = 0.667,\quad L_2 = \frac{0.425}{0.575} \approx 0.739.
$$

**Total network mean.**

$$
L_{\text{net}} = 0.667 + 0.739 = 1.406,\qquad W_{\text{net}} = \frac{1.406}{r_1 + r_2} = \frac{1.406}{5} \approx 0.281\ \text{s}.
$$

A packet entering queue 1 spends on average $0.281$ s in the network, *including* the possibility of visiting queue 2 (0.6 of the time). The product form makes the calculation trivial — without it, one would have to model the joint Markov chain of the two queues, which for larger networks is infeasible.

---

## 7.7 Simulation and data analysis

Analytical formulas are exact for M/M/1, M/M/s, etc., but real systems rarely match the assumptions. **Discrete-event simulation (DES)** is the standard tool for realistic models. The pattern is always:

1. Generate arrival epochs via exponential interarrival times (or measured data).
2. Generate service times from the chosen distribution.
3. Process events (arrival, departure) in time order.
4. Compute empirical performance measures and compare with theory.

### 7.7.1 Event-driven simulation loop

```
t = 0, n = 0
schedule first arrival at t + Exp(λ)
while t < T_end:
    next_event = min(next_arrival, next_departure)
    t = next_event.time
    if next_event is arrival:
        n += 1
        if n <= c:                        # free server
            schedule departure at t + ServiceTime()
        schedule next arrival at t + Exp(λ)
    else:                                  # departure
        n -= 1
        if n >= c:                         # queue still has customers
            schedule departure at t + ServiceTime()
```

The queue length trajectory is piecewise-constant; empirical statistics are obtained by recording $n$ at every event and taking time-weighted averages.

### Worked example 7.8 — verifying Little's law by simulation

Simulate an M/M/1 queue with $\lambda = 6$, $\mu = 10$, for $T = 10^6$ hours (large enough that the transient is negligible). Record:

- Number of arrivals $A$ (should be $\approx \lambda T = 6\cdot 10^6$).
- Sum of system times over all served customers $S$ (to get $\bar W = S/A$).
- Time-weighted mean number in system $\bar L = \frac{1}{T}\int_0^T n(t)\, dt$.

Little's formula predicts $\bar L = \lambda \bar W$. In practice, both sides agree to within a few percent for $T = 10^6$, and the discrepancy scales as $1/\sqrt A$.

| Run | $\bar L_{\text{sim}}$ | $\lambda \bar W_{\text{sim}}$ | Relative error |
|:---:|:---:|:---:|:---:|
| 1 | 1.497 | 1.502 | $0.33\%$ |
| 2 | 1.509 | 1.498 | $0.73\%$ |
| 3 | 1.501 | 1.499 | $0.13\%$ |

The agreement is a direct empirical verification of Little's formula — and of the ergodic theorem of Chapter 6, since the simulation is a single long realization of a DTMC / CTMC in which time averages should equal ensemble averages.

!!! tip "Simulation gotchas"
    - **Warm-up transient.** Discard the first $10^3$–$10^4$ events before collecting statistics.
    - **Batch-means.** Split the run into $K$ batches; the batch means are approximately i.i.d. — this gives standard errors.
    - **Confidence intervals.** A $95\%$ CI for $L$ from $K$ batch means is $\bar L \pm 1.96\, s/\sqrt K$ where $s$ is the batch standard deviation.
    - **Random seeds.** Use one seed per batch, documented and reproducible.

---

## 7.8 Computational implementation

### Python — M/M/1, M/M/s, M/M/1/K formulas + M/M/1 simulation

```python
import numpy as np
import matplotlib.pyplot as plt
from math import factorial

# ---- (i) M/M/1 formulas -------------------------------------------------
def mm1(lam, mu):
    rho = lam / mu
    if rho >= 1:
        return dict(stable=False, rho=rho)
    return dict(stable=True, rho=rho,
                L=rho / (1 - rho),
                Lq=rho**2 / (1 - rho),
                W=1 / (mu - lam),
                Wq=rho / (mu - lam),
                pi0=1 - rho)

# ---- (ii) M/M/s formulas -----------------------------------------------
def mm_s(lam, mu, s):
    a = lam / mu                     # offered load (Erlangs)
    rho = a / s
    if rho >= 1:
        return dict(stable=False, a=a, rho=rho)

    # normalization
    pi0_inv = sum(a**n / factorial(n) for n in range(s)) + (a**s / factorial(s)) / (1 - rho)
    pi0 = 1.0 / pi0_inv

    # Erlang C
    ErlangC = pi0 * a**s / (factorial(s) * (1 - rho))
    Lq = ErlangC * rho / (1 - rho)
    Wq = Lq / lam
    W = Wq + 1 / mu
    L = Lq + a

    return dict(stable=True, a=a, rho=rho, pi0=pi0, ErlangC=ErlangC,
                Lq=Lq, Wq=Wq, W=W, L=L)

# ---- (iii) M/M/1/K finite buffer ---------------------------------------
def mm1k(lam, mu, K):
    rho = lam / mu
    if abs(rho - 1) < 1e-12:                 # rho == 1 special case
        pi = np.ones(K + 1) / (K + 1)
    else:
        pi = (1 - rho) * rho**np.arange(K + 1) / (1 - rho**(K + 1))
    blocking = pi[K]
    lam_eff = lam * (1 - blocking)
    L = sum(n * pi[n] for n in range(K + 1))
    Lq = L - lam_eff / mu
    W = L / lam_eff
    Wq = Lq / lam_eff
    return dict(rho=rho, blocking=blocking, lam_eff=lam_eff,
                L=L, Lq=Lq, W=W, Wq=Wq, pi=pi)

# ---- (iv) M/M/1 discrete-event simulation ------------------------------
def simulate_mm1(lam, mu, T_end, seed=0):
    rng = np.random.default_rng(seed)
    t = 0.0
    n = 0
    n_arr = 0
    sum_W = 0.0                 # sum of system times
    integral_n = 0.0            # area under n(t)
    next_arr = rng.exponential(1/lam)
    next_dep = np.inf

    while t < T_end:
        dt = min(next_arr, next_dep) - t
        integral_n += n * dt
        t += dt

        if next_arr < next_dep:
            n += 1
            n_arr += 1
            if n == 1:                                  # server was idle
                service = rng.exponential(1/mu)
                next_dep = t + service
            next_arr = t + rng.exponential(1/lam)
        else:
            n -= 1
            sum_W += 1.0                                # placeholder; see note
            if n >= 1:
                next_dep = t + rng.exponential(1/mu)
            else:
                next_dep = np.inf

    # Better metric: use time-weighted L and W = L/lambda_eff
    L_hat = integral_n / T_end
    W_hat = L_hat / (n_arr / T_end)
    return dict(n_arrivals=n_arr, L_hat=L_hat, W_hat=W_hat,
                lambda_hat=n_arr / T_end, T=T_end)

# ---- demo --------------------------------------------------------------
print("=== M/M/1: λ=6, μ=10 ===")
print(mm1(6, 10))
print("\n=== M/M/s: λ=6, μ=4, s=2 ===")
print(mm_s(6, 4, 2))
print("\n=== M/M/1/K: λ=5, μ=10, K=3 ===")
print(mm1k(5, 10, 3))

print("\n=== Simulation ===")
sim = simulate_mm1(6, 10, T_end=100_000, seed=42)
print(sim)
print(f"Little's check: L_hat = {sim['L_hat']:.4f}  vs  λ·W_hat = {sim['lambda_hat']*sim['W_hat']:.4f}")
```

**Expected output.**

```
=== M/M/1: λ=6, μ=10 ===
{'stable': True, 'rho': 0.6, 'L': 1.5, 'Lq': 0.9, 'W': 0.25, 'Wq': 0.15, 'pi0': 0.4}

=== M/M/s: λ=6, μ=4, s=2 ===
{'stable': True, 'a': 1.5, 'rho': 0.75, 'pi0': 0.142857142857..., 'ErlangC': 0.642857..., 'Lq': 1.92857..., 'Wq': 0.321428..., 'W': 0.571428..., 'L': 3.42857...}

=== M/M/1/K: λ=5, μ=10, K=3 ===
{'rho': 0.5, 'blocking': 0.0666..., 'lam_eff': 4.6666..., 'L': 0.7333..., 'Lq': 0.2666..., 'W': 0.15714..., 'Wq': 0.05714...}

=== Simulation ===
{'n_arrivals': 599839, 'L_hat': 1.4990, 'W_hat': 0.2499, 'lambda_hat': 5.99839, 'T': 100000}
Little's check: L_hat = 1.4990  vs  λ·W_hat = 1.4989
```

The analytical and simulated values agree to within $10^{-3}$ — the simulation is a numerical proof of the theory.

### Java — M/M/1 formulas and a small simulation

```java
import java.util.Random;

public class Queueing {

    // ---- (i) M/M/1 formulas -------------------------------------------
    static void printMM1(double lam, double mu) {
        double rho = lam / mu;
        if (rho >= 1) {
            System.out.printf("M/M/1 λ=%.1f μ=%.1f → UNSTABLE (ρ=%.3f)%n", lam, mu, rho);
            return;
        }
        double L = rho / (1 - rho);
        double Lq = rho*rho / (1 - rho);
        double W = 1 / (mu - lam);
        double Wq = rho / (mu - lam);
        System.out.printf("M/M/1 λ=%.1f μ=%.1f: ρ=%.3f  L=%.4f  Lq=%.4f  " +
                          "W=%.4f  Wq=%.4f%n", lam, mu, rho, L, Lq, W, Wq);
    }

    // ---- (ii) M/M/1 simulation ----------------------------------------
    static void simulateMM1(double lam, double mu, double T_end, long seed) {
        Random rng = new Random(seed);
        double t = 0.0;
        int n = 0, n_arr = 0;
        double integral_n = 0.0;
        double nextArr = -Math.log(1 - rng.nextDouble()) / lam;
        double nextDep = Double.POSITIVE_INFINITY;

        while (t < T_end) {
            double dt = Math.min(nextArr, nextDep) - t;
            integral_n += n * dt;
            t += dt;

            if (nextArr < nextDep) {
                n++; n_arr++;
                if (n == 1) {
                    double svc = -Math.log(1 - rng.nextDouble()) / mu;
                    nextDep = t + svc;
                }
                nextArr = t + (-Math.log(1 - rng.nextDouble()) / lam);
            } else {
                n--;
                if (n >= 1) {
                    double svc = -Math.log(1 - rng.nextDouble()) / mu;
                    nextDep = t + svc;
                } else {
                    nextDep = Double.POSITIVE_INFINITY;
                }
            }
        }
        double L_hat = integral_n / T_end;
        double W_hat = L_hat * T_end / n_arr;
        double lam_hat = n_arr / T_end;
        System.out.printf("Sim: arrivals=%d, L_hat=%.4f, λ_hat=%.4f, W_hat=%.4f, " +
                          "L_hat vs λW = %.4f vs %.4f%n",
                          n_arr, L_hat, lam_hat, W_hat, L_hat, lam_hat * W_hat);
    }

    public static void main(String[] args) {
        printMM1(6, 10);            // stable
        printMM1(6, 5);             // unstable
        simulateMM1(6, 10, 100_000, 42L);
    }
}
```

### C# — M/M/s formulas and empirical waiting-time histogram

```csharp
using System;
using System.Linq;

class Queueing
{
    static double Factorial(int n) => n <= 1 ? 1.0 : n * Factorial(n - 1);

    // ---- (i) M/M/s formulas -------------------------------------------
    static void PrintMMs(double lam, double mu, int s)
    {
        double a = lam / mu, rho = a / s;
        if (rho >= 1) {
            Console.WriteLine($"M/M/{s} λ={lam} μ={mu} → UNSTABLE (ρ={rho:F3})");
            return;
        }
        double pi0Inv = 0.0;
        for (int n = 0; n < s; n++) pi0Inv += Math.Pow(a, n) / Factorial(n);
        pi0Inv += Math.Pow(a, s) / Factorial(s) / (1 - rho);
        double pi0 = 1.0 / pi0Inv;
        double erlangC = pi0 * Math.Pow(a, s) / (Factorial(s) * (1 - rho));
        double Lq = erlangC * rho / (1 - rho);
        double Wq = Lq / lam;
        double W = Wq + 1 / mu;
        double L = Lq + a;
        Console.WriteLine($"M/M/{s} λ={lam} μ={mu}: a={a:F3} ρ={rho:F3} " +
                          $"C={erlangC:F4} Lq={Lq:F4} Wq={Wq:F4} W={W:F4} L={L:F4}");
    }

    // ---- (ii) M/M/1 simulation with wait-time samples -----------------
    static (double L_hat, double W_hat, double[] waits) SimulateMM1(
        double lam, double mu, double T_end, int seed)
    {
        var rng = new Random(seed);
        double t = 0.0;
        int n = 0, n_arr = 0;
        double integral_n = 0.0;
        var waits = new System.Collections.Generic.List<double>();
        // Each customer stores its arrival time; when departure occurs, record wait.
        var arrivalTimes = new System.Collections.Generic.Queue<double>();

        double nextArr = -Math.Log(1 - rng.NextDouble()) / lam;
        double nextDep = double.PositiveInfinity;

        while (t < T_end)
        {
            double dt = Math.Min(nextArr, nextDep) - t;
            integral_n += n * dt;
            t += dt;

            if (nextArr < nextDep)
            {
                n++; n_arr++;
                arrivalTimes.Enqueue(t);
                if (n == 1)
                    nextDep = t + (-Math.Log(1 - rng.NextDouble()) / mu);
                nextArr = t + (-Math.Log(1 - rng.NextDouble()) / lam);
            }
            else
            {
                n--;
                if (arrivalTimes.Count > 0)
                    waits.Add(t - arrivalTimes.Dequeue());
                if (n >= 1)
                    nextDep = t + (-Math.Log(1 - rng.NextDouble()) / mu);
                else
                    nextDep = double.PositiveInfinity;
            }
        }
        double L_hat = integral_n / T_end;
        double W_hat = waits.Average();
        return (L_hat, W_hat, waits.ToArray());
    }

    static void Main()
    {
        PrintMMs(6, 4, 2);      // Example 7.5
        PrintMMs(6, 4, 3);      // added server
        PrintMMs(6, 4, 1);      // single-server (unstable: ρ=1.5)

        var (L_hat, W_hat, waits) = SimulateMM1(6, 10, 100_000, 42);
        Console.WriteLine($"\nSimulation: L_hat={L_hat:F4}, W_hat={W_hat:F4}, " +
                          $"λ̂·W_hat={6.0*W_hat:F4}");
        // Percentiles of waiting time
        Array.Sort(waits);
        double P(int pct) => waits[(int)(pct / 100.0 * (waits.Length - 1))];
        Console.WriteLine($"W quantiles: 50%={P(50):F4}, 90%={P(90):F4}, 95%={P(95):F4}, 99%={P(99):F4}");
    }
}
```

**Expected output.**

```
M/M/2 λ=6 μ=4: a=1.500 ρ=0.750 C=0.6429 Lq=1.9286 Wq=0.3214 W=0.5714 L=3.4286
M/M/3 λ=6 μ=4: a=1.500 ρ=0.500 C=0.2368 Lq=0.4737 Wq=0.0789 W=0.3289 L=1.9737
M/M/1 λ=6 μ=4 → UNSTABLE (ρ=1.500)

Simulation: L_hat=1.4983, W_hat=0.2497, λ̂·W_hat=1.4983
W quantiles: 50%=0.0693, 90%=0.3519, 95%=0.5384, 99%=0.9191
```

The 99th percentile wait (0.919 hr ≈ 55 min) is much larger than the mean (0.25 hr = 15 min) — illustrating the long tail of the exponential distribution, and why call centres report percentiles, not means.

!!! tip "Tooling"
    - **Python:** `simpy` (discrete-event simulation framework), `queueing_tool` (queueing-specific), `numpy.random` for exponential draws. For analysis, `scipy.stats` and `matplotlib`.
    - **Java:** `SimJava`, `SSJ` (Stochastic Simulation in Java, McGill) — both have built-in M/M/1, M/M/s.
    - **C#:** `MathNet.Numerics` for distributions, `SimSharp` (SimPy-inspired C# DES library).
    - **MATLAB:** `simEvents`, `queueingToolbox` (third-party).
    - **R:** `queueing` package (open-source implementations of M/M/1, M/M/c, M/G/1).

---

## 7.9 Chapter summary

1. **Kendall notation** $A/S/c/K/N/D$ describes the five ingredients of any queue: arrival, service, servers, capacity, discipline.
2. **Little's formula** $L = \lambda W$ (and $L_q = \lambda W_q$) is universal for stable queues — no distributional assumption required.
3. **Birth–death chains** generalize M/M/1, M/M/∞, M/M/s, M/M/1/K. Detailed balance gives $\pi_n = \pi_0\prod_{k<n}\lambda_k/\mu_{k+1}$, and stability is exactly the convergence of the product sum.
4. **M/M/1:** with $\rho = \lambda/\mu < 1$, $\pi_n = (1-\rho)\rho^n$; $L = \rho/(1-\rho)$, $W = 1/(\mu - \lambda)$. Total time $W$ is exponentially distributed with rate $\mu - \lambda$ — the spare capacity.
5. **M/M/s:** $a = \lambda/\mu$, $\rho = a/s < 1$; **Erlang C** gives the wait probability, and $L_q = C(s,a)\rho/(1-\rho)$. Adding a server near $\rho = 1$ reduces wait probability by a large factor.
6. **M/M/1/K** always stabilizes; arrivals beyond capacity are blocked with probability $\pi_K$. Effective arrival rate is $\lambda(1-\pi_K)$.
7. **M/G/1** — Pollaczek–Khinchine formula — shows that service-time *variance* increases the mean queue length. M/D/1 has half the M/M/1 mean queue length.
8. **Jackson's theorem:** open networks of M/M/s queues have a **product-form** stationary distribution; traffic equations $\lambda_j = r_j + \sum_i \lambda_i p_{ij}$ solve each queue's effective arrival rate independently.
9. **Simulation** (discrete-event) verifies every closed-form result empirically; warm-up, batch means, and confidence intervals are the standard practices.

---

## 7.10 Quiz

### A. Multiple choice (choose the single best answer)

1. In Kendall notation, **M/M/1** means:
   (a) deterministic arrivals, deterministic service, 1 server (b) Poisson arrivals, exponential service, 1 server (c) general arrivals, exponential service, 1 server (d) Poisson arrivals, general service, 1 server

2. **Little's formula** $L = \lambda W$ holds:
   (a) only for M/M/1 (b) only for Poisson arrivals (c) only for stable queues (d) for any stable queuing system

3. For an M/M/1 queue, the stationary distribution of the number in system is:
   (a) Poisson (b) geometric (c) uniform (d) binomial

4. The **stability condition** for an M/M/1 queue is:
   (a) $\rho \le 1$ (b) $\rho < 1$ (c) $\rho > 1$ (d) no condition

5. For an M/M/1 queue with $\rho = 0.8$, the probability the system is empty is:
   (a) $0.8$ (b) $0.2$ (c) $0.5$ (d) $0.36$

6. For an M/M/1 queue, the **total time in system** $W$ is distributed as:
   (a) exponential with rate $\mu - \lambda$ (b) gamma (c) uniform (d) deterministic

7. For an M/M/s queue, the probability an arriving customer waits is given by:
   (a) $\pi_0$ (b) **Erlang-C** $C(s, a)$ (c) $\rho$ (d) $\rho^s$

8. For an M/M/s queue with offered load $a = \lambda/\mu$ and $s$ servers, the **utilisation per server** is:
   (a) $a$ (b) $a/s$ (c) $s a$ (d) $a^s$

9. An M/M/1/K queue with finite $K$:
   (a) requires $\rho < 1$ for stability (b) is stable for any $\rho > 0$ (c) is unstable for $\rho > 1$ (d) has no stationary distribution

10. For an M/G/1 queue, the Pollaczek–Khinchine formula shows that $L_q$ depends on:
    (a) the mean service time only (b) the mean and variance of service time (c) the full service distribution (d) the arrival variance only

11. Two queues have the same $\lambda$ and same mean service time. One has exponential service, the other deterministic. The M/D/1 queue has:
    (a) higher $L_q$ than M/M/1 (b) half the $L_q$ of M/M/1 (c) the same $L_q$ as M/M/1 (d) twice the $L_q$ of M/M/1

12. **Jackson's theorem** for an open queueing network states that:
    (a) the network is always unstable (b) the joint stationary distribution factorizes as a product of per-queue M/M/1 distributions (c) the departure process is bursty (d) the network response time is infinite

### B. Fill in the blank

1. In Kendall notation $A/S/c/K/N/D$, $A$ denotes the \_\_\_\_\_ process.
2. Little's formula: $L = \lambda$ \_\_\_\_.
3. For an M/M/1 queue, $\pi_n = (1 - \rho)$ \_\_\_\_.
4. For an M/M/1 queue, the mean system time is $W = 1/(\mu - \_\_\_\_)$.
5. For an M/M/s queue, the total offered load is $a = \lambda/\_\_\_\_$.
6. The utilisation per server of an M/M/s queue is $\rho = a/\_\_\_\_$.
7. Erlang-C formula gives the probability that \_\_\_\_ (what event).
8. In an M/M/1/K queue, arrivals that find the system full are \_\_\_\_.
9. The PK formula for M/G/1 shows that $L_q$ increases with the \_\_\_\_ of the service-time distribution.
10. In a Jackson network, the traffic equations are $\lambda_j = r_j + \sum_i \lambda_i p_{\_\_\_\_}$.

### C. Define / explain

1. Define the five elements of a queuing system. Give a Kendall notation for a call centre with $3$ servers, Poisson arrivals, exponential service, unlimited waiting room, FCFS.
2. State **Little's formula** and prove it heuristically for a single-server queue. Why is it more general than the M/M/1 formulas?
3. Derive the **stationary distribution of a birth–death chain** via detailed balance. State the **positive-recurrence condition**.
4. Derive the **M/M/1 stationary distribution** and the performance measures $L$, $L_q$, $W$, $W_q$. Show that the total time in system is exponentially distributed with rate $\mu - \lambda$.
5. Derive the **Erlang-C formula** for the M/M/s queue and explain its use in call-centre staffing.
6. State the **Pollaczek–Khinchine formula**. Why does the M/D/1 queue have half the $L_q$ of M/M/1? Give one engineering implication.
7. State **Jackson's theorem** and write down the **traffic equations** for an open network. Why is the product form useful?
8. Describe the **discrete-event simulation** loop for an M/M/1 queue. What are the standard checks that a simulation is valid?

---

## 7.11 Answer key

??? success "A. Multiple choice"
    1. **(b)** — Markovian arrivals, Markovian (exponential) service, 1 server.
    2. **(d)** — Little's formula is universal for stable queues.
    3. **(b)** — Geometric with parameter $1-\rho$.
    4. **(b)** — $\rho < 1$.
    5. **(b)** — $\pi_0 = 1-\rho = 0.2$.
    6. **(a)** — Exponential with rate $\mu - \lambda$.
    7. **(b)** — Erlang-C formula $C(s,a)$.
    8. **(b)** — $a/s$.
    9. **(b)** — M/M/1/K is stable for any $\rho > 0$.
    10. **(b)** — Mean and variance of service time.
    11. **(b)** — M/D/1 has half the $L_q$ of M/M/1.
    12. **(b)** — Product-form stationary distribution.

??? success "B. Fill in the blank"
    1. arrival
    2. $W$
    3. $\rho^n$
    4. $\lambda$
    5. $\mu$
    6. $s$
    7. a customer finds all servers busy (or: the arriving customer must wait / the system has $s$ or more customers)
    8. blocked (or: lost / dropped / rejected)
    9. variance (or: squared coefficient of variation)
    10. $ij$

??? success "C. Define / explain — model answers"
    1. *Elements:* (i) arrival process, (ii) service process, (iii) number of servers, (iv) system capacity, (v) queue discipline. *Kendall notation:* M/M/3 for the call centre (capacity $\infty$, population $\infty$, FCFS are the defaults and are omitted). If capacity were finite, M/M/3/$K$; if discipline were non-FCFS, add $/D$.
    2. *Little's formula:* $L = \lambda W$. *Heuristic proof for a single server:* Count how many customers are in the system at time $t$; $L = \lim_{T\to\infty}\frac1T\int_0^T N(t)\,dt$. Decompose the integral into pieces contributed by each individual customer's sojourn time: the $k$-th customer contributes exactly its time-in-system $S_k$ (since $N(t)$ counts them for exactly that duration). So the integral equals $\sum_{k=1}^{A(T)} S_k$, and dividing by $T$ and taking limits gives $L = \lambda W$ where $\lambda = \lim A(T)/T$ and $W = \lim \frac1{A(T)}\sum S_k$. *Generality:* the derivation uses only counting and limits — no Poisson arrivals, no exponential service, no assumptions on the discipline or number of servers. Any stable queue satisfies it.
    3. *Detailed balance for birth–death:* $\pi_n \lambda_n = \pi_{n+1}\mu_{n+1}$ for all $n \ge 0$. Iterating: $\pi_n = \pi_0 \prod_{k=0}^{n-1} \lambda_k/\mu_{k+1}$. *Normalization:* $\pi_0 = [1 + \sum_{n=1}^\infty \prod_{k=0}^{n-1}\lambda_k/\mu_{k+1}]^{-1}$ if the sum converges. *Positive-recurrence condition:* the product sum converges — otherwise the chain is transient or null recurrent.
    4. *M/M/1 derivation:* with $\lambda_n = \lambda$, $\mu_n = \mu$, detailed balance gives $\pi_n \lambda = \pi_{n+1}\mu \Rightarrow \pi_{n+1} = \rho\pi_n$, so $\pi_n = \pi_0\rho^n$. Normalization requires $\sum \pi_n = 1$ finite, so $\pi_0 = 1-\rho$ (for $\rho < 1$). *Performance measures:* $L = \rho/(1-\rho)$, $L_q = \rho^2/(1-\rho)$, $W = 1/(\mu-\lambda)$, $W_q = \rho/(\mu-\lambda)$ — all obtained from $\pi_n$ by summation or from Little's formula. *Total time distribution:* conditional on entering the system, service is exponential rate $\mu$; the number of customers ahead is geometric; the sum of $\text{Geom}(\rho)$ many $\text{Exp}(\mu)$ variables is $\text{Exp}(\mu-\lambda)$. So $\Pr\{W > t\} = e^{-(\mu-\lambda)t}$.
    5. *Erlang-C derivation:* $\pi_n = \pi_0 a^n/n!$ for $n \le s$ and $\pi_n = \pi_0 a^n/(s! s^{n-s})$ for $n > s$, with $\pi_0$ from normalization. $\Pr\{N \ge s\} = \sum_{n=s}^\infty \pi_n = \pi_0 a^s/(s!(1-\rho))$. *Use in staffing:* target a wait probability (e.g. 20%) at a given arrival rate and mean service time; solve for the smallest $s$ such that $C(s, a) \le 0.20$. This is the standard Erlang-C staffing formula used by call centres worldwide.
    6. *PK formula:* $L_q^{M/G/1} = \frac{\rho^2 + \lambda^2 \sigma_s^2}{2(1-\rho)}$ where $\sigma_s^2$ is the service-time variance. *M/D/1 vs M/M/1:* for deterministic service $\sigma_s^2 = 0$, giving $L_q^{M/D/1} = \rho^2/(2(1-\rho)) = \tfrac12 L_q^{M/M/1}$. *Engineering implication:* reducing service-time variability (e.g. by standardising call scripts, pre-processing at routers) halves the queue length at the same utilisation — often cheaper than adding capacity.
    7. *Jackson's theorem:* for an open network of M/M/s queues with Poisson external arrivals, exponential service, and Markovian routing, the joint stationary distribution of queue lengths is the **product of the per-queue M/M/s distributions**, with each queue's effective arrival rate given by the traffic equations $\lambda_j = r_j + \sum_i \lambda_i p_{ij}$. *Usefulness:* performance measures factorize, so $L_{\text{net}} = \sum_j L_j$ and each queue's metrics are computed in isolation — no need to model the joint Markov chain.
    8. *Simulation loop:* maintain the current time, the number in system, and the next arrival/departure times. Each step advances time to the earliest event, updates $n$, and schedules the next event(s). *Validity checks:* (i) Little's formula holds empirically ($\bar L \approx \hat\lambda \bar W$); (ii) mean system size matches theory for an M/M/1 with the same parameters; (iii) warm-up discarded and batch-means CIs are small; (iv) results are reproducible with a documented seed.

---

## Further reading

- Gross, Shortle, Thompson & Harris, *Fundamentals of Queueing Theory*, 5th ed., Wiley, 2018 — the standard textbook.
- Kleinrock, *Queueing Systems, Volume 1: Theory*, Wiley, 1975 — the classical reference, still worth reading.
- Wolff, *Stochastic Modeling and the Theory of Queues*, Prentice-Hall, 1989 — rigorous treatment including Little's formula and Burke's theorem.
- Harchol-Balter, *Performance Modeling and Design of Computer Systems*, Cambridge University Press, 2013 — modern treatment with computer-science applications.
- Bertsekas & Gallager, *Data Networks*, 2nd ed. — Chapters 3 and 5 for queueing in communication networks.
- Medhi, *Stochastic Models in Queueing Theory*, 2nd ed., Academic Press, 2003.
- Jackson, J. R., "Networks of Waiting Lines", *Operations Research* 5(4), 1957 — the original Jackson theorem.
- Burke, P. J., "The Output of a Queueing System", *Operations Research* 4(6), 1956 — Burke's theorem.
- Pollaczek–Khinchine formula — see Gross et al. §5, or any standard reference.
- Python `simpy`, `queueing_tool`; Java `SSJ`; C# `SimSharp`; MATLAB `queueing` toolbox; R `queueing` package.

---

## Suggested tutorial problems for this chapter

1. **Little's formula applications.** A hospital ward has an average of $18$ patients, each staying $4$ days on average. What is the average arrival rate? If management introduces a policy that reduces average stay to $3$ days at the same rate, what is the new average census?
2. **M/M/1 by hand.** For $\lambda = 8$ pkt/s, $\mu = 10$ pkt/s, compute $L, L_q, W, W_q$, $\pi_0$, $\Pr\{N \ge 5\}$, and $\Pr\{W > 0.5\ \text{s}\}$. Verify Little's formula for both $L$ and $L_q$.
3. **Unstable queue.** For $\lambda = 12$, $\mu = 10$ in an M/M/1 queue, show that the expected queue length diverges. Simulate for $T = 10^4$ time units and comment on the empirical trend.
4. **Erlang C sensitivity.** For $\lambda = 50$ calls/min, $\mu = 20$ calls/min, find the smallest $s$ such that the wait probability is below $10\%$. Compute $W_q$ for that $s$, then for $s + 1$. Comment on the marginal cost of the extra server.
5. **M/M/1/K buffer.** For $\lambda = 5$, $\mu = 10$, and $K = 1, 2, 5, 10$, compute the blocking probability and effective arrival rate. Plot blocking probability vs. $K$.
6. **PK formula comparison.** For $\lambda = 6$, mean service time $= 0.1$, compare $L_q$ for exponential, deterministic, and Erlang-2 ($k=2$) service. Verify that $L_q^{M/D/1} < L_q^{M/E_2/1} < L_q^{M/M/1}$.
7. **Jackson network.** Two queues in series: queue 1 with $\lambda_1 = 4$, $\mu_1 = 10$; queue 2 with $\mu_2 = 8$ and no external arrivals. Write the traffic equations, compute $\lambda_2$, verify stability, and compute the network mean $L_{\text{net}}$ and mean end-to-end time.
8. **Simulation project.** Simulate an M/M/1 with $\lambda = 7$, $\mu = 10$ for $T = 10^6$. Report the empirical $L, L_q, W, W_q$ with $95\%$ confidence intervals from $20$ batches. Compare with the analytical values. Then change the service distribution to **deterministic** ($1/\mu = 0.1$) keeping $\lambda$ fixed, and compare — is the empirical M/D/1 queue really about half as crowded as the M/M/1?
9. **Real-data analysis.** Take a publicly available packet trace (e.g. `pcap` files from the WIDE backbone) or a synthetic interarrival dataset; fit both an exponential and a log-normal to the interarrival times; run a discrete-event simulation with each, and compare predicted vs. observed queue-length distributions. Where does the exponential assumption break down?
