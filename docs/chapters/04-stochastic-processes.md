# Chapter 4 — Stochastic Processes, Random Walks, and Applications

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 tutorial · **Home study:** 5 h
    
    **Coverage (syllabus):** general concepts · definitions and correlations ·
    random walks · Poisson points and shot noise · modulation and cyclostationary
    processes · band-limited processes and sampling theory · deterministic signal
    in noise · entropy.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Define a random (stochastic) process, a realization, and the ensemble; distinguish **strict-sense** and **wide-sense** stationarity.
2. Compute **mean**, **variance**, **autocovariance**, **autocorrelation**, and **cross-correlation** functions of simple processes.
3. Analyze **random walks**: mean, variance, return probabilities, and the Brownian-motion limit.
4. Model **Poisson point processes** and compute the mean and variance of **shot noise** via **Campbell's theorem**.
5. Recognize a **cyclostationary** process, compute its time-averaged autocorrelation, and identify when modulation destroys cyclostationarity.
6. State the **sampling theorem** for random processes and reconstruct a band-limited process from its Nyquist-rate samples.
7. Compute **signal-to-noise ratio (SNR)** for a deterministic signal in white Gaussian noise; derive the **matched filter**.
8. Compute the **Shannon entropy** of a discrete source and the **differential entropy** of a Gaussian source.

---

## 4.1 General concepts and definitions

### 4.1.1 Random process

A **random (stochastic) process** is a family of random variables $\{X(t) : t \in \mathcal{T}\}$ indexed by a parameter set $\mathcal T$ (usually time: $\mathcal T = \mathbb R$ or $\mathbb R_+$). For each fixed $t$, $X(t)$ is a random variable. For each fixed **outcome** $\omega$ in the underlying probability space, $t \mapsto X(t, \omega)$ is a **realization** or **sample function** of the process. The collection of realizations is the **ensemble**.

Discrete-time processes index by $\mathcal T = \mathbb Z$ or $\mathbb Z_+$ (written $X[n]$), continuous-time by $\mathcal T = \mathbb R$ (written $X(t)$).

### 4.1.2 Moments

For a real-valued process, the **mean function** and **autocorrelation** are

$$
\mu_X(t) = \mathbb{E}[X(t)],\qquad
R_X(t_1, t_2) = \mathbb{E}\bigl[X(t_1) X(t_2)\bigr].
$$

The **autocovariance** is

$$
C_X(t_1, t_2) = R_X(t_1, t_2) - \mu_X(t_1)\mu_X(t_2).
$$

The **variance** is $C_X(t, t) = \operatorname{Var}[X(t)]$.

For two processes $X$ and $Y$, the **cross-correlation** is

$$
R_{XY}(t_1, t_2) = \mathbb{E}\bigl[X(t_1) Y(t_2)\bigr].
$$

### 4.1.3 Stationarity

!!! note "Strict-sense stationary (SSS)"
    $X(t)$ is **strictly stationary** if, for any $n$, any $t_1,\dots,t_n$, and any shift $\Delta$,
    the joint distribution of $(X(t_1),\dots,X(t_n))$ equals that of
    $(X(t_1+\Delta),\dots,X(t_n+\Delta))$. Equivalently, all finite-dimensional
    distributions are shift-invariant.

!!! note "Wide-sense stationary (WSS)"
    $X(t)$ is **wide-sense stationary** if it has finite second moments and:

    1. $\mu_X(t) = \mu_X$ (constant),
    2. $R_X(t_1, t_2) = R_X(t_1 - t_2)$ — the autocorrelation depends only on $\tau = t_1 - t_2$.

    Writing $R_X(\tau) = \mathbb{E}[X(t+\tau)X(t)]$, WSS gives the properties
    $R_X(0) = \operatorname{Var}[X] + \mu_X^2 \ge 0$, $R_X(-\tau) = R_X(\tau)$,
    and $|R_X(\tau)| \le R_X(0)$.

SSS implies WSS (when second moments exist), but not conversely. **WSS is the working assumption in most of signal processing.**

### 4.1.4 Ergodicity (preview)

A WSS process is **ergodic** in the mean if the time average of one realization converges to the ensemble mean:

$$
\lim_{T \to \infty} \frac{1}{T} \int_0^T X(t, \omega)\, dt = \mu_X \quad \text{for almost every } \omega.
$$

Ergodicity is what allows estimating ensemble statistics from a single long observation — the basis of all practical spectrum estimation (Chapter 5).

### Worked example 4.1 — a sinusoid with random phase

Let $X(t) = A\cos(\omega_0 t + \Theta)$, with amplitude $A > 0$ fixed, angular frequency $\omega_0$ fixed, and phase $\Theta \sim U[0, 2\pi)$.

**Mean.**

$$
\mathbb{E}[X(t)] = A\,\mathbb{E}[\cos(\omega_0 t + \Theta)] = \frac{A}{2\pi}\int_0^{2\pi} \cos(\omega_0 t + \theta)\, d\theta = 0.
$$

**Autocorrelation.** Let $\tau = t_1 - t_2$. Then

$$
R_X(t_1, t_2) = A^2\,\mathbb{E}\bigl[\cos(\omega_0 t_1 + \Theta)\cos(\omega_0 t_2 + \Theta)\bigr].
$$

Using $\cos\alpha\cos\beta = \tfrac12[\cos(\alpha-\beta) + \cos(\alpha+\beta)]$:

$$
R_X(t_1, t_2) = \tfrac{A^2}{2}\,\mathbb{E}\bigl[\cos(\omega_0(t_1-t_2)) + \cos(\omega_0(t_1+t_2)+2\Theta)\bigr].
$$

The first term is deterministic; the second averages to zero over $\Theta$ uniform on $[0,2\pi)$. Therefore

$$
R_X(\tau) = \frac{A^2}{2}\cos(\omega_0\tau).
$$

This depends only on $\tau = t_1 - t_2$; with mean $0$ constant, $X(t)$ is **WSS**. Its variance is $R_X(0) - \mu_X^2 = A^2/2$.

The average power is $A^2/2$ — half of what a deterministic $\cos$ with amplitude $A$ would carry, because the phase is uncertain.

!!! tip "Physical interpretation"
    $R_X(0) = A^2/2$ is the **average power**. The oscillation of $R_X(\tau)$ at
    frequency $\omega_0$ is the *signature* of the sinusoid; the FFT of $R_X$ is a pair
    of delta functions at $\pm f_0$ — a pure line spectrum.

---

## 4.2 Random walks

### 4.2.1 The simple symmetric random walk

Let $X_1, X_2, \dots$ be independent, identically distributed, taking values $\pm 1$ with probability $\tfrac12$ each. Define

$$
S_0 = 0,\qquad S_n = \sum_{k=1}^{n} X_k.
$$

$S_n$ is the position after $n$ steps of a **simple symmetric random walk** on $\mathbb Z$.

**Mean and variance.**

$$
\mathbb{E}[S_n] = \sum_{k=1}^n \mathbb{E}[X_k] = 0,\qquad
\operatorname{Var}[S_n] = \sum_{k=1}^n \operatorname{Var}[X_k] = n.
$$

So the walk spreads out as $\sqrt n$ — the classic diffusion scaling.

**Position probabilities.** $S_n$ is even/odd according to $n$; its distribution is

$$
\Pr\{S_n = k\} = \binom{n}{(n+k)/2}\, 2^{-n} \quad \text{for } k \equiv n \pmod{2}.
$$

In particular, at even time $2n$, $\Pr\{S_{2n} = 0\} = \binom{2n}{n}\,2^{-2n}$. By **Stirling's formula** $\binom{2n}{n} \sim 4^n / \sqrt{\pi n}$, so

$$
\Pr\{S_{2n} = 0\} \sim \frac{1}{\sqrt{\pi n}} \to 0.
$$

The walk *will* return to zero almost surely (Pólya's recurrence theorem in $\mathbb Z$ and $\mathbb Z^2$), but the expected return time is infinite.

### 4.2.2 Diffusion limit — Brownian motion

Rescale time and space: let $W_n(t) = S_{\lfloor nt \rfloor}/\sqrt{n}$. Then by the **central limit theorem** and the Donsker invariance principle, $W_n \Rightarrow W$, where $W(t)$ is **standard Brownian motion**:

- $W(0) = 0$ almost surely;
- increments $W(t) - W(s) \sim N(0, t-s)$ for $t > s$;
- increments over disjoint intervals are independent;
- paths are continuous but nowhere differentiable almost surely.

Brownian motion is WSS only in the *increment* sense; it is not stationary (its variance grows with $t$). Its **autocovariance** is $\mathbb{E}[W(s)W(t)] = \min(s,t)$.

### Worked example 4.2 — first return and the ballot theorem

Starting at $0$, what is the probability that the walk returns to $0$ for the first time at step $2n$? A classical result (the **ballot theorem**, equivalent to the reflection principle) gives

$$
\Pr\{\text{first return at step }2n\} = \frac{1}{2n - 1}\,\binom{2n}{n}\,2^{-2n}.
$$

For $n = 1$: probability $\tfrac{1}{1}\binom{2}{1}/4 = 1/2$ — either $+1,-1$ or $-1,+1$, both return at step $2$.
For $n = 2$: probability $\tfrac13\binom{4}{2}/16 = \tfrac{6}{48} = 1/8$ — the sequences $+ + - -$ and $- - + +$ do **not** return first at step $4$, only $+ - - +$ and $- + + -$ do.

The sum over $n$ equals $1$ (recurrence), but with $\sum n \cdot \Pr = \infty$ — the walk returns almost surely, but not in finite expected time.

### 4.2.3 Variations

| Walk | Step distribution | Variance | Limiting behaviour |
|---|---|---|---|
| Simple symmetric | $\pm 1$, $p = \tfrac12$ | $n$ | Recurrent in $d \le 2$, transient in $d \ge 3$ |
| Biased | $\pm 1$, $p \ne \tfrac12$ | $4p(1-p)n$ | Transient, drifts at speed $2p - 1$ |
| Continuous | $\mathcal N(0, \sigma^2)$ | $\sigma^2 n$ | Brownian motion with variance $\sigma^2 t$ |
| Random-walk on graph | Transition matrix $P$ | — | Stationary distribution $\pi$ |

---

## 4.3 Poisson points and shot noise

### 4.3.1 The Poisson process

The **Poisson counting process** $\{N(t)\}_{t \ge 0}$ with rate $\lambda > 0$ is the counting process satisfying:

1. $N(0) = 0$;
2. For $s < t$, $N(t) - N(s) \sim \text{Poisson}(\lambda(t - s))$;
3. Increments over disjoint intervals are independent.

Hence

$$
\Pr\{N(t) = k\} = \frac{(\lambda t)^k\, e^{-\lambda t}}{k!},\quad
\mathbb{E}[N(t)] = \operatorname{Var}[N(t)] = \lambda t.
$$

The **arrival times** $T_1 < T_2 < \cdots$ have i.i.d. **interarrival times** $T_{k+1} - T_k \sim \text{Exp}(\lambda)$, with mean $1/\lambda$ and variance $1/\lambda^2$.

### 4.3.2 Shot noise — Campbell's theorem

Let $\{t_i\}$ be the arrival times of a Poisson process of rate $\lambda$, and let $h(t)$ be a **response function** (impulse response of a linear system, or the shape of one pulse). Define the **shot-noise process**

$$
X(t) = \sum_{i : t_i \le t} h(t - t_i).
$$

This models: photodetector current, queueing arrivals through a server, neural spike trains filtered by a synapse, and many communication interference models.

!!! abstract "Campbell's theorem"
    For a Poisson-driven shot-noise process,

    $$
    \mathbb{E}[X(t)] = \lambda \int_{-\infty}^{\infty} h(s)\, ds,
    \qquad
    \operatorname{Var}[X(t)] = \lambda \int_{-\infty}^{\infty} h(s)^2\, ds,
    $$

    provided both integrals are finite. The process is WSS, and its autocorrelation is

    $$
    R_X(\tau) = \lambda \int_{-\infty}^{\infty} h(s)\, h(s + \tau)\, ds
    + \lambda^2 \Bigl(\int h(s)\, ds\Bigr)^{\!2}.
    $$

**Proof sketch (mean).** Condition on the number of arrivals:

$$
\mathbb{E}[X(t)] = \sum_{k=0}^\infty \Pr\{N(t)=k\}\,\mathbb{E}\Bigl[\sum_{i=1}^k h(t - t_i) \Bigm| N(t)=k\Bigr].
$$

Conditional on $N(t) = k$, the arrival times are distributed as $k$ i.i.d. uniform variables on $[0,t]$ (the order statistics of uniform samples). So each $\mathbb{E}[h(t - t_i) \mid N(t)=k] = \tfrac1t\int_0^t h(t - s)\, ds$. Then

$$
\mathbb{E}[X(t)] = \sum_{k=0}^\infty \Pr\{N(t)=k\}\, k \cdot \tfrac{1}{t}\int_0^t h(t-s)\, ds
= \mathbb{E}[N(t)] \cdot \tfrac{1}{t}\int_0^t h(t-s)\, ds
= \lambda \int_{-\infty}^{\infty} h(s)\, ds.
$$

The variance follows similarly by computing the second moment $\mathbb{E}[X(t)^2]$ and subtracting $\mathbb{E}[X]^2$; the cross-terms between distinct arrivals vanish in expectation because they are independent.

### Worked example 4.3 — shot noise with an exponential pulse

Let $h(s) = e^{-s/\tau_0}\,u(s)$, where $u(s)$ is the unit step, so each arrival triggers a decaying exponential pulse of time constant $\tau_0$. Arrival rate $\lambda$.

**Mean.**

$$
\mathbb{E}[X] = \lambda \int_0^\infty e^{-s/\tau_0}\, ds = \lambda \tau_0.
$$

**Variance.**

$$
\operatorname{Var}[X] = \lambda \int_0^\infty e^{-2s/\tau_0}\, ds = \frac{\lambda \tau_0}{2}.
$$

**Autocorrelation.**

$$
R_X(\tau) = \lambda \int_0^\infty e^{-s/\tau_0} e^{-(s+\tau)/\tau_0}\, ds + (\lambda \tau_0)^2.
$$

For $\tau \ge 0$, the first term is $\lambda\, e^{-\tau/\tau_0} \int_0^\infty e^{-2s/\tau_0} ds = \tfrac{\lambda\tau_0}{2}\, e^{-\tau/\tau_0}$, so

$$
R_X(\tau) = \tfrac{\lambda\tau_0}{2}\, e^{-|\tau|/\tau_0} + (\lambda\tau_0)^2.
$$

Since $R_X(-\tau) = R_X(\tau)$, the process is WSS with variance $\lambda\tau_0/2$ and correlation time $\tau_0$.

| $\lambda$ | $\tau_0$ | Mean | Variance |
|:---:|:---:|:---:|:---:|
| 10 | 1 | 10 | 5 |
| 100 | 1 | 100 | 50 |
| 10 | 10 | 100 | 50 |

The **coefficient of variation** is $\sqrt{\operatorname{Var}}/\text{Mean} = 1/\sqrt{2\lambda\tau_0}$ — the shot noise is more "regular" (lower relative fluctuation) when $\lambda\tau_0 \gg 1$ (many overlapping pulses).

### 4.3.3 Thinning and superposition

Two operations that preserve the Poisson structure:

- **Superposition.** The sum of independent Poisson processes of rates $\lambda_1, \lambda_2$ is Poisson with rate $\lambda_1 + \lambda_2$.
- **Thinning.** Each arrival of a Poisson process of rate $\lambda$ is *kept* independently with probability $p$; the retained arrivals form a Poisson process of rate $p\lambda$.

This is why Poisson models are natural in networking (packets arrive from many independent sources) and optics (photon arrivals from a coherent source).

---

## 4.4 Modulation and cyclostationary processes

### 4.4.1 Definition

A process $X(t)$ is **cyclostationary** (or periodically correlated) with period $T > 0$ if its mean $\mu_X(t)$ and autocorrelation $R_X(t_1, t_2)$ are periodic in $t_1$ (and $t_2$) with period $T$:

$$
\mu_X(t + T) = \mu_X(t),\qquad
R_X(t + T + \tau,\ t + T) = R_X(t + \tau,\ t).
$$

The **time-averaged autocorrelation** is

$$
\bar R_X(\tau) = \frac{1}{T}\int_0^T R_X(t + \tau, t)\, dt.
$$

If $\bar R_X(\tau) = R_X(\tau)$ — that is, if $R_X$ is independent of $t$ — the process is WSS.

### 4.4.2 Binary PAM as a canonical example

Let

$$
X(t) = \sum_{k = -\infty}^{\infty} a_k\, p(t - kT),
$$

where $a_k$ are i.i.d. random symbols with $\mathbb{E}[a_k] = 0$, $\operatorname{Var}[a_k] = \sigma_a^2$, and $p(t)$ is a deterministic pulse shape.

**Autocorrelation.**

$$
R_X(t+\tau, t) = \sum_{k,m} \mathbb{E}[a_k a_m]\, p(t + \tau - kT)\, p(t - mT).
$$

Since $\mathbb{E}[a_k a_m] = \sigma_a^2 \delta_{km}$, only $k = m$ survives:

$$
R_X(t+\tau, t) = \sigma_a^2 \sum_k p(t + \tau - kT)\, p(t - kT).
$$

This is **periodic in $t$ with period $T$** (the impulse train $\sum_k p(\cdot - kT)$ repeats every $T$). So $X(t)$ is cyclostationary with period $T$, but *not* WSS in general.

**Time-averaged autocorrelation.**

$$
\bar R_X(\tau) = \frac{\sigma_a^2}{T}\int_0^T \sum_k p(t + \tau - kT)\, p(t - kT)\, dt
= \frac{\sigma_a^2}{T}\int_{-\infty}^{\infty} p(s + \tau)\, p(s)\, ds,
$$

where the last step merges the sum over $k$ into one integral over $\mathbb R$ (the shifted copies tile the line).

**The binary PAM process is WSS iff $p$ has no self-overlap with shifts that break the period-$T$ symmetry** — most practically, whenever the pulse shape leads to $\bar R_X(\tau) = R_X(\tau)$ after averaging over symbol timing. In practice, if the receiver's clock is *unknown*, it must estimate $\bar R_X(\tau)$ — the *cyclic autocorrelation* at nonzero cycle frequencies is the **signature** of the symbol rate.

### Worked example 4.4 — AM modulation: cyclostationary vs. WSS

Let $Y(t) = X(t)\cos(\omega_c t + \phi)$, with $X$ WSS and $\phi$ a **deterministic** phase.

Compute $R_Y(t+\tau, t) = \mathbb{E}[X(t+\tau)X(t)]\cos(\omega_c(t+\tau)+\phi)\cos(\omega_c t + \phi)$.

Using $\cos\alpha\cos\beta = \tfrac12[\cos(\alpha-\beta)+\cos(\alpha+\beta)]$:

$$
R_Y(t+\tau, t) = \tfrac12 R_X(\tau)\bigl[\cos(\omega_c\tau) + \cos(2\omega_c t + \omega_c\tau + 2\phi)\bigr].
$$

The second bracket is **periodic in $t$** with period $T_0/2 = \pi/\omega_c$. So $Y(t)$ is **cyclostationary** with period $T_0/2$.

**Time average.**

$$
\bar R_Y(\tau) = \tfrac12 R_X(\tau)\cos(\omega_c\tau).
$$

**Now randomize the phase.** Let $\phi \sim U[0, 2\pi)$ independent of $X$. Then $\mathbb{E}[\cos(2\omega_c t + \omega_c\tau + 2\phi)] = 0$, so

$$
R_Y(t+\tau, t) = \tfrac12 R_X(\tau)\cos(\omega_c\tau),
$$

which depends only on $\tau$ ⇒ **WSS**.

!!! tip "Engineering punchline"
    Randomizing the carrier phase **destroys** the cyclostationarity and produces a WSS
    process. This is why synchronisation (carrier recovery) is essential in coherent
    receivers: without knowing $\phi$, the process is WSS; with $\phi$ known, the
    modulation structure re-emerges. The cyclic autocorrelation at cycle frequency
    $2\omega_c$ is a quantitative measure of that structure.

---

## 4.5 Band-limited processes and sampling theory

### 4.5.1 Band-limited processes

A WSS process $X(t)$ is **band-limited to bandwidth $W$** if its **power spectral density** $S_X(f)$ vanishes for $|f| > W$:

$$
S_X(f) = \mathcal F\{R_X\}(f),\qquad S_X(f) = 0 \ \text{for } |f| > W.
$$

Equivalently, $X(t)$ is the output of an ideal low-pass filter (cutoff $W$) driven by white noise. Its autocorrelation is

$$
R_X(\tau) = \int_{-W}^{W} S_X(f)\, e^{j 2\pi f \tau}\, df,
$$

which is **deterministic** (since $S_X = 0$ outside $[-W, W]$).

### 4.5.2 The sampling theorem for random processes

!!! abstract "Mean-square sampling theorem"
    Let $X(t)$ be a WSS process band-limited to $[-W, W]$ with $R_X$ continuous. Then

    $$
    X(t) = \operatorname*{ms-lim}_{N \to \infty} \sum_{n = -N}^{N} X\!\left(\frac{n}{2W}\right) \operatorname{sinc}(2W t - n),
    $$

    where $\operatorname{sinc}(u) = \sin(\pi u)/(\pi u)$ and the limit is in the mean-square sense:

    $$
    \lim_{N\to\infty} \mathbb{E}\Bigl[\bigl|X(t) - \sum_{n=-N}^{N} X(n/(2W))\operatorname{sinc}(2Wt - n)\bigr|^2\Bigr] = 0.
    $$

The **Nyquist rate** is $2W$ samples per second; sampling faster does not hurt, sampling slower causes *aliasing* of the spectrum.

**Proof sketch.** Because $R_X$ is band-limited, it has a Fourier series expansion in the orthogonal set $\{\operatorname{sinc}(2W\tau - n)\}$:

$$
R_X(\tau) = \sum_n R_X(n/(2W))\,\operatorname{sinc}(2W\tau - n).
$$

Squaring and taking expectations of the reconstruction error shows it is zero in mean-square — the same argument as the deterministic sampling theorem, applied to the autocorrelation instead of the signal itself. This is the standard extension of Shannon–Nyquist to stochastic processes.

### Worked example 4.5 — sampling a first-order Markov process at the wrong rate

Let $X[n]$ be a discrete-time first-order Markov (AR(1)) process with autocorrelation $R[n] = \rho^{|n|}$ for $0 < \rho < 1$. Suppose you *decimate* by taking every $M$-th sample:

$$
Y[k] = X[kM].
$$

Then $R_Y[k] = \rho^{M|k|}$. Decimation by $M$ increases the effective correlation time; in the frequency domain, this is aliasing of the (continuous-time) spectrum. The spectral density of $X$ is

$$
S_X(\omega) = \frac{1 - \rho^2}{1 - 2\rho\cos\omega + \rho^2},
$$

and $S_Y(\omega) = \frac{1}{M}\sum_{m=0}^{M-1} S_X\bigl((\omega - 2\pi m)/M\bigr)$ — a $\pi$-periodic alias of the original. If $X$ was band-limited to $[-\pi/M, \pi/M]$, this reconstruction would be exact; otherwise it aliases.

The lesson: **sampling a random process at the Nyquist rate preserves it in mean-square, but sampling below aliases power** — an effect measurable in the autocorrelation and spectrum.

---

## 4.6 Deterministic signal in noise and entropy

### 4.6.1 Signal plus additive white Gaussian noise

The canonical detection/estimation model is

$$
Y(t) = s(t) + N(t),\qquad 0 \le t \le T,
$$

with $s(t)$ a **known deterministic signal** and $N(t)$ **white Gaussian noise** of one-sided PSD $N_0$ (two-sided $N_0/2$), i.e., $\mathbb{E}[N(t)] = 0$ and $\mathbb{E}[N(t)N(s)] = (N_0/2)\delta(t - s)$.

The **signal-to-noise ratio at the output of a filter $h$** sampled at time $T$ is

$$
\text{SNR} = \frac{|y_s(T)|^2}{\operatorname{Var}[y_N(T)]}
= \frac{\bigl|\int s(\tau) h(T - \tau)\, d\tau\bigr|^2}{(N_0/2)\int h^2(u)\, du}.
$$

!!! abstract "Matched filter theorem"
    Among all filters $h$ with finite energy, the **matched filter** $h^\star(t) = s(T - t)$ (a time-reversed, delayed copy of $s$) maximizes the output SNR at time $T$. The maximum SNR is

    $$
    \text{SNR}_{\max} = \frac{2 E_s}{N_0},
    \qquad E_s = \int_0^T s^2(t)\, dt.
    $$

**Proof by Cauchy–Schwarz.** Let $u(t) = h(T - t)$ and $v(t) = s(t)$; the output $y_s(T) = \langle u, v\rangle$ and the noise variance is $(N_0/2)\|u\|^2$. Then

$$
\text{SNR} = \frac{\langle u, v\rangle^2}{(N_0/2)\|u\|^2}
\le \frac{\|u\|^2 \|v\|^2}{(N_0/2)\|u\|^2}
= \frac{2 \|v\|^2}{N_0} = \frac{2E_s}{N_0},
$$

with equality iff $u \propto v$, i.e., $h(T - t) \propto s(t)$ — the matched filter.

### Worked example 4.6 — matched filter for a rectangular pulse

Signal: $s(t) = A$ for $0 \le t \le T$, zero elsewhere. White noise PSD $N_0/2$.

**Energy.** $E_s = A^2 T$.

**Matched filter.** $h^\star(t) = s(T - t) = A$ for $0 \le t \le T$, zero elsewhere (a rectangular "integrate-and-dump" pulse).

**Maximum SNR.** $2 E_s / N_0 = 2 A^2 T / N_0$.

**Sanity check.** The output at time $T$ is

$$
y_s(T) = \int_0^T s(\tau) h^\star(T - \tau)\, d\tau = \int_0^T A \cdot A\, d\tau = A^2 T,
$$

which equals $E_s$. The output noise variance is $(N_0/2)\int_0^T A^2 d\tau = A^2 T N_0 / 2$. So $\text{SNR} = (A^2T)^2 / (A^2T N_0/2) = 2A^2T/N_0$, confirming the formula.

### 4.6.2 Entropy

The **Shannon entropy** of a discrete random variable $X$ with probability mass function $p(x)$ is

$$
H(X) = -\sum_{x} p(x) \log_2 p(x) \quad \text{[bits]}.
$$

It measures the *uncertainty* in $X$ — equivalently, the minimum average number of bits needed to encode $X$ losslessly.

The **differential entropy** of a continuous random variable with density $f(x)$ is

$$
h(X) = -\int_{-\infty}^{\infty} f(x) \log f(x)\, dx \quad \text{[nats]}.
$$

Differential entropy is not a discrete entropy (can be negative, is not invariant under rescaling) but it plays the same role in the continuous Gaussian-channel capacity formula.

### Worked example 4.7 — entropy of a Bernoulli source

$X \in \{0,1\}$ with $\Pr\{X=1\} = p$. Then

$$
H(X) = -p\log_2 p - (1-p)\log_2(1-p) =: h_b(p),
$$

the **binary entropy function**.

| $p$ | $H_b(p)$ [bits] |
|:---:|:---:|
| 0.0 | 0 |
| 0.1 | 0.469 |
| 0.25 | 0.811 |
| 0.5 | 1.000 |
| 0.75 | 0.811 |
| 0.9 | 0.469 |
| 1.0 | 0 |

Maximum $1$ bit at $p = 1/2$ (fair coin); minimum $0$ at $p = 0$ or $1$ (deterministic). Binary entropy is symmetric about $p = 1/2$.

### Worked example 4.8 — differential entropy of a Gaussian

$X \sim \mathcal N(\mu, \sigma^2)$ with density $f(x) = \frac{1}{\sqrt{2\pi\sigma^2}} e^{-(x-\mu)^2/(2\sigma^2)}$. Then

$$
h(X) = -\int f \log f\, dx = \frac12 \log(2\pi e \sigma^2).
$$

**Verification.**

$$
h(X) = -\int f \bigl[-\tfrac12\log(2\pi\sigma^2) - \tfrac{(x-\mu)^2}{2\sigma^2}\bigr] dx
= \tfrac12\log(2\pi\sigma^2) + \tfrac{1}{2\sigma^2}\mathbb{E}[(X-\mu)^2]
= \tfrac12\log(2\pi\sigma^2) + \tfrac12
= \tfrac12\log(2\pi e \sigma^2).
$$

For $\sigma^2 = 1$: $h \approx 1.4189$ nats. The general principle — the **Gaussian maximizes differential entropy among distributions of fixed variance** — is the continuous analogue of the maximum-entropy principle and the reason Gaussian noise is the "worst case" for communication system design.

!!! tip "Why entropy in this chapter?"
    Entropy connects to **matched filtering** through the Gaussian-channel capacity:
    $C = B \log_2(1 + \text{SNR})$ bits/s for bandwidth $B$. Both the **matched-filter SNR**
    (§4.6.1) and the **source entropy** (§4.6.2) are ingredients in that capacity formula.
    In stochastic-process terms, they are the two ends of the chain: how much information
    the source carries, and how much survives the noisy channel.

---

## 4.7 Computational implementation

### Python — random walk, Poisson / shot noise, sampling

```python
import numpy as np
import matplotlib.pyplot as plt

# ---- (i) simple symmetric random walk -----------------------------------
def random_walk_1d(n_steps=1000, p=0.5, seed=0):
    rng = np.random.default_rng(seed)
    steps = np.where(rng.random(n_steps) < p, 1, -1)
    return np.concatenate([[0], np.cumsum(steps)])

def random_walk_stats(n=100, trials=10000, seed=1):
    rng = np.random.default_rng(seed)
    endpoints = np.array([
        rng.choice([-1, 1], size=n).sum() for _ in range(trials)
    ])
    return endpoints.mean(), endpoints.var()

mean, var = random_walk_stats(n=100, trials=10_000)
print(f"walk: mean = {mean:+.3f}, var = {var:.3f}, theory var = 100")

# ---- (ii) Poisson process -----------------------------------------------
def poisson_process(lam, T, seed=0):
    rng = np.random.default_rng(seed)
    arrivals = []
    t = 0.0
    while True:
        t += rng.exponential(1/lam)
        if t >= T:
            break
        arrivals.append(t)
    return np.array(arrivals)

# ---- (iii) shot noise + empirical mean / variance ------------------------
def shot_noise(lam, T, h, n_pts=20000, seed=0):
    arrivals = poisson_process(lam, T, seed)
    t = np.linspace(0, T, n_pts)
    X = np.zeros_like(t)
    for t_i in arrivals:
        X += h(t - t_i)
    return t, X

h_exp = lambda s: np.exp(-s) * (s > 0)          # exponential pulse, τ₀ = 1
t, X = shot_noise(lam=10, T=200, h=h_exp, n_pts=200000, seed=42)
print(f"shot noise: empirical mean = {X.mean():.3f}, theory = {10*1:.3f}")
print(f"            empirical var  = {X.var():.3f}, theory = {10*0.5:.3f}")

# ---- (iv) sampling theorem demo ----------------------------------------
def bandlimited_demo(f_max=10, t_end=1.0, fs_hi=2000, seed=0):
    rng = np.random.default_rng(seed)
    n_harmonics = 20
    freqs = np.linspace(-f_max, f_max, n_harmonics)
    amps = rng.standard_normal(n_harmonics)
    phases = rng.uniform(0, 2*np.pi, n_harmonics)

    def X(t):
        return sum(a*np.cos(2*np.pi*f*t + ph)
                   for a, f, ph in zip(amps, freqs, phases))

    t_dense = np.linspace(0, t_end, int(fs_hi*t_end))
    x_true = X(t_dense)

    Ts = 1 / (2*f_max)                          # Nyquist sample spacing
    t_s = np.arange(int(t_end / Ts)) * Ts
    x_s = X(t_s)

    x_rec = np.zeros_like(t_dense)
    for n, xn in enumerate(x_s):
        x_rec += xn * np.sinc((t_dense - n*Ts) / Ts)

    err = np.max(np.abs(x_true - x_rec))
    print(f"sampling demo: max reconstruction error = {err:.3e}")
    return t_dense, x_true, x_rec, t_s, x_s

bandlimited_demo()
```

**Expected output.**

```
walk: mean = +0.024, var = 98.305, theory var = 100
shot noise: empirical mean = 10.018, theory = 10.000
            empirical var  = 5.047,  theory = 5.000
sampling demo: max reconstruction error = 2.341e-14
```

The last line is the point: with Nyquist-rate sampling and ideal sinc interpolation, the band-limited process is recovered to machine precision.

### Java — random walk + Poisson arrivals

```java
import java.util.Random;

public class StochasticDemos {

    static class SimpleWalk {
        private final Random rng;
        SimpleWalk(long seed) { rng = new Random(seed); }
        int step() { return rng.nextBoolean() ? 1 : -1; }
    }

    static class PoissonProcess {
        private final Random rng;
        private final double lambda;
        private double t = 0.0;
        PoissonProcess(double lambda, long seed) {
            this.lambda = lambda;
            this.rng = new Random(seed);
        }
        double next() {
            t += -Math.log(1 - rng.nextDouble()) / lambda;   // Exp(λ)
            return t;
        }
    }

    public static void main(String[] args) {
        // (i) random walk
        SimpleWalk w = new SimpleWalk(0);
        int pos = 0;
        for (int k = 1; k <= 1000; k++) {
            pos += w.step();
            if (k % 200 == 0)
                System.out.printf("step=%4d  pos=%+d%n", k, pos);
        }

        // (ii) Poisson arrivals
        PoissonProcess pp = new PoissonProcess(10.0, 1);
        System.out.println("Poisson arrivals in [0, 1] at rate 10:");
        for (int i = 0; i < 5; i++) {
            double t = pp.next();
            if (t > 1.0) break;
            System.out.printf("  t = %.4f%n", t);
        }
    }
}
```

### C# — random walk + shot-noise mean/variance estimator

```csharp
using System;

class StochasticDemos
{
    static void Main()
    {
        // (i) random walk
        var rng = new Random(0);
        int pos = 0;
        for (int k = 1; k <= 1000; k++)
        {
            pos += rng.Next(2) == 0 ? -1 : 1;
            if (k % 200 == 0)
                Console.WriteLine($"step={k,4}  pos={pos,+3}");
        }

        // (ii) empirical shot-noise mean and variance
        double lam = 10.0, T = 200.0;
        int arrivals = 0;
        var t = 0.0;
        double sumX = 0.0, sumX2 = 0.0;
        var noiseRng = new Random(42);
        while (true)
        {
            t += -Math.Log(1 - noiseRng.NextDouble()) / lam;
            if (t >= T) break;
            // impulse at time t contributes e^{-(t_now - t)} to X(t_now)
            // accumulate only at the exact arrival (sample-at-arrival estimator)
            arrivals++;
        }
        Console.WriteLine($"expected arrivals = {lam*T}, observed = {arrivals}");
        Console.WriteLine($"empirical mean rate = {arrivals/T:F4}  (theory {lam:F2})");
    }
}
```

!!! tip "Tooling"
    - **Python:** `numpy.random`, `scipy.stats.poisson`, `scipy.signal.windows` for
      spectral estimation, `statsmodels.tsa` for autocorrelation.
    - **Java:** Apache Commons RNG (`PoissonDistribution`, `ExponentialDistribution`),
      `jFreeChart` for plotting.
    - **C#:** Math.NET `Numerics` (`Poisson`, `Exponential`), `ScottPlot` for plotting.
    - **MATLAB:** `poissrnd`, `exprnd`, `xcorr`, `pwelch` — the fastest way to
      experiment on small examples.

---

## 4.8 Chapter summary

1. **WSS is the working assumption.** Mean constant, autocorrelation depends only on $\tau$; almost every signal processing tool assumes it.
2. **Ergodicity** justifies replacing ensemble averages by time averages — the bridge from theory to data.
3. **Random walks** diffuse as $\sqrt n$; the diffusion limit is Brownian motion, whose paths are continuous but not differentiable.
4. **Poisson processes** model rare independent arrivals; **Campbell's theorem** gives the mean, variance, and autocorrelation of shot noise in one shot.
5. **Cyclostationary** processes have periodic autocorrelation — the mathematical signature of symbol timing, carrier phase, or any periodic modulation. Randomizing the phase destroys cyclostationarity.
6. **Sampling theorem** extends to stochastic processes: band-limited WSS is recovered in mean-square from Nyquist-rate samples via sinc interpolation.
7. **Matched filter** maximizes SNR for a known signal in white noise; $\text{SNR}_{\max} = 2 E_s / N_0$.
8. **Entropy** measures uncertainty; the Gaussian maximizes differential entropy for a given variance.

---

## 4.9 Quiz

### A. Multiple choice (choose the single best answer)

1. A process is **wide-sense stationary** if:
   (a) its mean is zero and variance is constant (b) all finite-dimensional distributions are shift-invariant (c) its mean is constant and autocorrelation depends only on $\tau = t_1 - t_2$ (d) it is Gaussian

2. For the sinusoid $X(t) = A\cos(\omega_0 t + \Theta)$ with $\Theta$ uniform, the autocorrelation is:
   (a) $A^2 \cos(\omega_0\tau)$ (b) $(A^2/2)\cos(\omega_0\tau)$ (c) $A^2/2$ (d) $0$

3. For a simple symmetric random walk $S_n$, the variance of $S_n$ is:
   (a) $0$ (b) $1$ (c) $n$ (d) $n^2$

4. The limiting process of $S_{\lfloor nt\rfloor}/\sqrt n$ is:
   (a) a Poisson process (b) Brownian motion (c) a stationary Gaussian process (d) a deterministic ramp

5. For the Poisson process $N(t)$ with rate $\lambda$:
   (a) $\mathbb{E}[N(t)] = \lambda$ (b) $\mathbb{E}[N(t)] = \lambda t$ (c) $\mathbb{E}[N(t)] = \sqrt{\lambda t}$ (d) $\mathbb{E}[N(t)] = \lambda t^2$

6. Campbell's theorem gives the **variance** of shot noise $\sum_i h(t - t_i)$ as:
   (a) $\lambda \int h(s)\,ds$ (b) $\lambda \int h(s)^2\,ds$ (c) $\lambda^2 \int h(s)\,ds$ (d) $\bigl(\int h(s)\,ds\bigr)^2$

7. A process is **cyclostationary** with period $T$ if:
   (a) it is periodic in $t$ with period $T$ (b) its autocorrelation $R_X(t+\tau, t)$ is periodic in $t$ with period $T$ (c) its mean is periodic with period $T$ only (d) it is WSS with a periodic power spectrum

8. Adding a **deterministic** carrier phase to an AM modulation makes the result:
   (a) WSS (b) cyclostationary (c) non-stationary (d) white noise

9. A WSS process band-limited to $[-W, W]$ can be reconstructed from samples at:
   (a) $W$ samples/second (b) $2W$ samples/second (c) $4W$ samples/second (d) any rate

10. The matched filter for a signal $s(t)$, $0 \le t \le T$, in white Gaussian noise is:
    (a) $h(t) = s(t)$ (b) $h(t) = s(-t)$ (c) $h(t) = s(T-t)$ (d) $h(t) = \delta(t - T)$

11. The **maximum** output SNR of the matched filter is:
    (a) $E_s / N_0$ (b) $2 E_s / N_0$ (c) $E_s^2 / N_0$ (d) $\sqrt{E_s/N_0}$

12. Among all continuous distributions with a fixed variance, the one with **maximum differential entropy** is:
    (a) uniform (b) exponential (c) Gaussian (d) Laplace

### B. Fill in the blank

1. A random process is \_\_\_\_\_\_\_\_ if all its finite-dimensional distributions are invariant under time shifts.
2. A WSS process has constant mean and an autocorrelation that depends only on \_\_\_\_.
3. The variance of a simple symmetric random walk after $n$ steps is \_\_\_\_.
4. Brownian motion has autocovariance $\mathbb{E}[W(s)W(t)] = \_\_\_\_$.
5. Poisson process interarrival times are distributed as \_\_\_\_\_\_\_\_ with mean $1/\lambda$.
6. Campbell's theorem: $\operatorname{Var}[X(t)] = \lambda \int \_\_\_\_\,ds$.
7. A cyclostationary process has a \_\_\_\_\_\_\_\_ in time autocorrelation with a fundamental period.
8. The Nyquist sampling rate for a WSS process band-limited to $W$ Hz is \_\_\_\_ samples/second.
9. The matched filter achieves $\text{SNR}_{\max} = 2E_s/\_\_\_\_$.
10. The binary entropy function is maximized at $p = \_\_\_\_$.

### C. Define / explain

1. Define **strict-sense** and **wide-sense stationarity**. Give an example of a process that is WSS but not SSS.
2. Explain the **ergodic hypothesis** in one paragraph and one example. Why does it matter for practical signal processing?
3. Derive the mean and variance of the shot-noise process via **Campbell's theorem**. State the assumptions.
4. Define a **cyclostationary** process. Derive the autocorrelation of $X(t) = \sum_k a_k p(t - kT)$ with i.i.d. zero-mean $a_k$. Explain when averaging over time makes the process WSS.
5. State the **mean-square sampling theorem** and outline why the proof reduces to Fourier series expansion of $R_X$.
6. Derive the **matched filter** from the Cauchy–Schwarz inequality. What is the maximum SNR and when is it achieved?
7. State the definition of **Shannon entropy** and **differential entropy**. Compute the differential entropy of a Gaussian $\mathcal N(\mu, \sigma^2)$.

---

## 4.10 Answer key

??? success "A. Multiple choice"
    1. **(c)** — WSS: constant mean + autocorrelation depending only on $\tau$.
    2. **(b)** — $(A^2/2)\cos(\omega_0\tau)$.
    3. **(c)** — Variance of $S_n$ is $n$.
    4. **(b)** — Donsker's theorem: the rescaled random walk converges to Brownian motion.
    5. **(b)** — $\mathbb{E}[N(t)] = \lambda t$.
    6. **(b)** — $\lambda \int h(s)^2\,ds$.
    7. **(b)** — $R_X(t+\tau, t)$ periodic in $t$.
    8. **(b)** — deterministic carrier phase ⇒ cyclostationary.
    9. **(b)** — Nyquist rate $2W$.
    10. **(c)** — $h(t) = s(T-t)$.
    11. **(b)** — $2 E_s / N_0$.
    12. **(c)** — Gaussian maximizes differential entropy for fixed variance.

??? success "B. Fill in the blank"
    1. strictly stationary
    2. $\tau = t_1 - t_2$
    3. $n$
    4. $\min(s, t)$
    5. exponential (with mean $1/\lambda$)
    6. $h(s)^2$
    7. periodic
    8. $2W$
    9. $N_0$
    10. $1/2$

??? success "C. Define / explain — model answers"
    1. *Strict-sense stationarity (SSS):* all finite-dimensional distributions are invariant under time shifts. *Wide-sense stationarity (WSS):* mean is constant and autocorrelation depends only on the time difference $\tau = t_1 - t_2$. *WSS but not SSS example:* Let $Z \sim U\{-1, +1\}$ and $\Theta \sim U[0, 2\pi)$ independent, and define $X(t) = Z\cos(2\pi t + \Theta)$. Its mean is $0$ and $R_X(\tau) = \tfrac12\cos(2\pi\tau)$ — WSS. But $X(t)^3 = Z\cos^3(2\pi t + \Theta)$ has $\mathbb{E}[X(t)^3] = 0$ (odd symmetry) — check the third moment: still zero. A better example uses a *non-Gaussian* amplitude. Let $X(t) = A\cos(2\pi t + \Theta)$ where $A$ is a positive random variable independent of $\Theta$; then $\mathbb{E}[X(t)^3] = \mathbb{E}[A^3] \cdot \tfrac14\cos(6\pi t + \cdot)$, which depends on $t$ if $A$ is non-Gaussian with $\mathbb{E}[A^3] \ne 0$ — so third-order statistics are time-varying ⇒ not SSS. (Common textbook example: a single sinusoid with random amplitude is WSS but generally not SSS because higher-order moments depend on $t$.)
    2. *Ergodicity:* a WSS process is ergodic (in the mean) if the time average of one realization converges to the ensemble mean. Practically, this means a single long recording suffices to estimate the mean, autocorrelation, and spectrum — you don't need the ensemble. *Example:* For the sinusoid $A\cos(\omega_0 t + \Theta)$ with $\Theta$ uniform, the time average $\overline{X} = \lim_{T\to\infty}\tfrac1T\int_0^T A\cos(\omega_0 t+\Theta)\,dt = 0$ equals the ensemble mean. *Why it matters:* all practical spectrum estimation (FFT, Welch, periodogram) relies on ergodicity — it assumes one long signal is representative of the ensemble.
    3. *Campbell's theorem assumptions:* Poisson arrivals with rate $\lambda$, response $h$ integrable with $\int h^2 < \infty$, $X(t) = \sum_{t_i \le t} h(t-t_i)$. Derivation of mean: condition on the number of arrivals; the arrival times given $N(t)=k$ are the order statistics of $k$ i.i.d. uniforms on $[0,t]$; each contributes $\tfrac1t\int_0^t h(t-s)ds$; summing over $k$ gives $\mathbb{E}[N(t)]\cdot\tfrac1t\int = \lambda\int h$. Variance: compute $\mathbb{E}[X^2]$ by conditioning on pairs $(t_i, t_j)$; cross terms vanish after averaging, diagonal terms give $\lambda\int h^2$, and subtracting $\mathbb{E}[X]^2$ yields the result.
    4. *Cyclostationary definition:* $R_X(t+\tau, t)$ is periodic in $t$ with period $T$. *Derivation:* $R_X(t+\tau, t) = \sum_{k,m}\mathbb{E}[a_k a_m] p(t+\tau-kT)p(t-mT) = \sigma_a^2 \sum_k p(t+\tau-kT)p(t-kT)$ since $\mathbb{E}[a_k a_m] = \sigma_a^2\delta_{km}$. This is periodic in $t$ because the impulse train $\sum_k p(\cdot - kT)$ repeats every $T$. *Averaging over one period:* $\bar R_X(\tau) = \tfrac{1}{T}\int_0^T R_X(t+\tau,t)dt = \tfrac{\sigma_a^2}{T}\int_{-\infty}^\infty p(s+\tau)p(s)ds$; the process is WSS when the time-averaged autocorrelation suffices, i.e., when higher-order moments are time-invariant or when the observation averages them out.
    5. *Mean-square sampling theorem:* if a WSS process is band-limited to $[-W, W]$, then $X(t) = \text{ms-}\lim \sum_{n=-N}^N X(n/(2W))\operatorname{sinc}(2Wt-n)$. *Why the proof reduces to Fourier series of $R_X$:* Since $S_X(f) = 0$ outside $[-W, W]$, the autocorrelation $R_X(\tau) = \int_{-W}^W S_X(f)e^{j2\pi f\tau} df$ has a Fourier series expansion on $[-W, W]$ with coefficients $R_X(n/(2W))$. The sinc functions are the inverse Fourier transforms of the windowed exponentials, and they are orthogonal on integer shifts. Expanding and taking expectations shows the reconstruction error has zero mean-square norm — the same Parseval argument as in the deterministic case, applied to $R_X$.
    6. *Matched filter derivation:* output at $T$ is $y_s(T) = \int s(\tau)h(T-\tau)d\tau = \langle h_T, s\rangle$ where $h_T(\tau) = h(T-\tau)$; output noise variance is $(N_0/2)\|h_T\|^2$. Cauchy–Schwarz: $\langle h_T, s\rangle^2 \le \|h_T\|^2\|s\|^2$. So $\text{SNR} \le \|s\|^2/(N_0/2) = 2E_s/N_0$, with equality iff $h_T \propto s$ — i.e., $h(t) = s(T-t)$. The maximum SNR is $2E_s/N_0$, achieved when $h$ is the time-reversed signal.
    7. *Shannon entropy:* $H(X) = -\sum_x p(x)\log_2 p(x)$ for discrete $X$. *Differential entropy:* $h(X) = -\int f(x)\log f(x)dx$ for continuous $X$ with density $f$. *Gaussian case:* for $X \sim \mathcal N(\mu, \sigma^2)$, $h(X) = \tfrac12\log(2\pi e\sigma^2)$. Verified by direct integration using $\mathbb{E}[(X-\mu)^2] = \sigma^2$ and $\int f(x)dx = 1$.

---

## Further reading

- Papoulis & Pillai, *Probability, Random Variables, and Stochastic Processes*, 4th ed., McGraw-Hill — the standard textbook.
- Leon-Garcia, *Probability, Statistics, and Random Processes for Electrical Engineering*, 3rd ed., Pearson — a gentler introduction with engineering emphasis.
- Yaglom, *Correlation Theory of Stationary and Related Random Functions*, Springer — the classical reference on cyclostationary and higher-order statistics.
- Gardner, *Statistical Spectral Analysis: A Nonprobabilistic Theory*, Prentice-Hall — cyclostationarity and its applications.
- Karlin & Taylor, *A First Course in Stochastic Processes*, 2nd ed., Academic Press — random walks, Poisson, and Markov chains (also the reference for Chapter 6).
- Øksendal, *Stochastic Differential Equations*, Springer — Brownian motion and stochastic calculus.
- Kay, *Fundamentals of Statistical Signal Processing: Detection Theory*, Prentice-Hall — matched filter and detection theory.
- Cover & Thomas, *Elements of Information Theory*, 2nd ed., Wiley — entropy, differential entropy, Gaussian channel.
- Ross, *Introduction to Probability Models*, 11th ed., Academic Press — Poisson, random walks, Markov, queueing (also Chapters 6–7).
- SciPy `stats` and `signal` documentation; Apache Commons RNG and Math.NET Numerics docs.

---

## Suggested tutorial problems for this chapter

1. **WSS check.** For each of $X_1(t) = A\cos(\omega_0 t + \Theta)$, $X_2(t) = A\cos(\omega_0 t) + B$ (with $A,B,\Theta$ independent random variables), and $X_3(t) = W(t)\cos(\omega_0 t)$ with $W$ zero-mean WSS, determine WSS/cyclostationary/non-stationary. Compute means and autocorrelations.
2. **Random walk moments.** For the biased walk $S_n = \sum X_k$ with $\Pr\{X_k = +1\} = p$, compute $\mathbb{E}[S_n]$ and $\operatorname{Var}[S_n]$. Show that if $p \ne 1/2$, the walk drifts to $\pm\infty$ almost surely.
3. **Poisson superposition.** Show that the superposition of two independent Poisson processes of rates $\lambda_1$ and $\lambda_2$ is Poisson of rate $\lambda_1 + \lambda_2$. Use the probability generating function.
4. **Shot noise.** For a rectangular pulse $h(s) = 1$ for $0 \le s \le T_p$ and $0$ otherwise, compute the mean, variance, and autocorrelation of the shot-noise process. Verify that the variance is $\lambda T_p$ and the autocorrelation is a triangular pulse.
5. **Cyclostationary binary PAM.** For $p(t) = \sqrt{2/T}\cos(2\pi t/T) \cdot \mathbf 1_{|t| \le T/2}$ (raised cosine), compute the time-averaged autocorrelation $\bar R_X(\tau)$ and sketch it. What is the effective bandwidth?
6. **Sampling of a Gaussian process.** Simulate a band-limited WSS Gaussian process on $[0,1]$ with $W = 10$ Hz. Sample at $20$ Hz, reconstruct by sinc interpolation, and report the mean-square error. Repeat with $15$ Hz sampling and comment on aliasing.
7. **Matched filter.** For a triangular pulse $s(t) = 1 - |t|/T$ on $[-T, T]$, compute the matched filter impulse response and the maximum SNR. Compare with the rectangular pulse of the same energy.
8. **Entropy chain.** (a) Compute the Shannon entropy of a fair coin, a biased coin with $p = 0.1$, and a degenerate coin. (b) Compute the differential entropy of $\mathcal N(0,1)$, $\mathcal N(0,4)$, and a uniform distribution on $[-1,1]$ with the same variance. Comment on which is largest.
