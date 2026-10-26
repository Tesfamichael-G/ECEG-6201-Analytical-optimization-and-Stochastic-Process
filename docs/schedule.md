# Schedule

Week-by-week plan for **ECEG-6201**, 18-week semester (including revision
and final-examination weeks).

!!! info "Reading before each lecture"
    The **Reading** column in the detailed breakdown below lists the sections to
    read *before* the lecture — approximately 45–60 minutes of preparation per
    lecture. The lectures assume you have at least skimmed the section; they
    focus on worked examples, derivations, and computer demonstrations, not on
    restating the text.

---

## At a glance

| Week | Chapter | Topic | Tutorial | Assessment |
|:---:|:---:|---|---|---|
| 1 | 1 | Preliminaries — characterization of minima, Pareto | — | — |
| 2 | 1 | Gradient / Hessian, optimality conditions, line search | T1 | A1 out |
| 3 | 2 | Unconstrained — direct search, Newton | — | — |
| 4 | 2 | Unconstrained — BFGS, Gauss–Newton, LM | T2 | A1 due / A2 out |
| 5 | 3 | Constrained — penalty, barrier | — | — |
| 6 | 3 | Constrained — augmented Lagrangian, QP | T3 | A2 due / A3 out |
| 7 | — | Review + mid-semester preparation | — | A3 due |
| 8 | — | **Mid-semester examination** | — | **30%** |
| 9 | 4 | Stochastic processes — WSS, correlations, random walks | — | A4 out |
| 10 | 4 | Poisson points, shot noise, cyclostationarity, sampling | T4 | A4 due / A5 out |
| 11 | 5 | Spectral representation — PSD, KL, factorization | — | — |
| 12 | 5 | Spectrum estimation, Wiener, Kalman | T5 | A5 due / A6 out |
| 13 | 6 | Markov chains — DTMC, classification, stationarity | — | — |
| 14 | 6 | CTMC, ergodicity, birth–death | T6 | A6 due / A7 out |
| 15 | 7 | Queuing theory — Little's formula, M/M/1 | — | — |
| 16 | 7 | M/M/s, Erlang-C, Jackson networks, simulation | T7 | A7 due |
| 17 | — | Revision week | — | — |
| 18 | — | **Final examination** | — | **50%** |

**Legend.** T = tutorial session (3 h) · A = assignment · — = no activity that week

---

## Detail by week

### Week 1 — Chapter 1, Part 1

- **Lecture.** Characterization of minima; local vs. global vs. Pareto;
  gradient and Hessian; slope and curvature.
- **Reading.** Chapter 1, §1.1–§1.5.
- **Preparation.** Review matrix algebra (eigenvalues, positive definiteness).

### Week 2 — Chapter 1, Part 2

- **Lecture.** Optimality conditions (unconstrained, equality, inequality, KKT);
  line-search subproblem; bracketing and golden-section; numerically vs.
  analytically derived gradients; adjoint variables.
- **Tutorial 1.** Worked examples 1.1–1.13, gradient check, golden-section
  implementation.
- **Reading.** Chapter 1, §1.6–§1.13.
- **Assignment 1** released: gradient checks, line search, adjoint derivation.

### Week 3 — Chapter 2, Part 1

- **Lecture.** Direct-search methods (pattern search, Nelder–Mead);
  second-order methods (Newton); damped Newton.
- **Reading.** Chapter 2, §2.1–§2.3.

### Week 4 — Chapter 2, Part 2

- **Lecture.** Quasi-Newton (BFGS); sums of squares; Gauss–Newton;
  Levenberg–Marquardt.
- **Tutorial 2.** BFGS update by hand; GN and LM on the exponential fit.
- **Reading.** Chapter 2, §2.4–§2.8.
- **Assignment 1 due; Assignment 2** released: implement Newton, BFGS, and LM.

### Week 5 — Chapter 3, Part 1

- **Lecture.** Penalty functions; logarithmic barrier; conditioning and
  duality-gap analysis; Slater's condition.
- **Reading.** Chapter 3, §3.1–§3.3.

### Week 6 — Chapter 3, Part 2

- **Lecture.** Augmented Lagrangian (method of multipliers); quadratic programming
  with equality constraints (KKT system); active-set preview.
- **Tutorial 3.** Penalty convergence, AL iteration, QP by KKT.
- **Reading.** Chapter 3, §3.4–§3.7.
- **Assignment 2 due; Assignment 3** released: penalty, AL, and QP.

### Week 7 — Review and mid-semester preparation

- **Lecture.** Review of Chapters 1–3; worked exam-style problems; common pitfalls.
- **Tutorial.** Q&A session; past mid-semester papers.
- **Reading.** Review all of Chapters 1–3.
- **Assignment 3 due.**

### Week 8 — Mid-semester examination

- **Exam.** 2 hours, closed book, **30% of final grade**. Covers Chapters 1–3.
- **No lecture this week.**

### Week 9 — Chapter 4, Part 1

- **Lecture.** Random (stochastic) processes; strict-sense and wide-sense
  stationarity; autocorrelation; cross-correlation; ergodicity preview.
- **Reading.** Chapter 4, §4.1–§4.2.
- **Assignment 4** released: WSS checks, random-walk moments.

### Week 10 — Chapter 4, Part 2

- **Lecture.** Poisson point processes; shot noise and Campbell's theorem;
  cyclostationarity; band-limited processes and the sampling theorem;
  matched filter and entropy.
- **Tutorial 4.** Shot-noise mean/variance, cyclostationarity of PAM,
  Nyquist-rate sampling.
- **Reading.** Chapter 4, §4.3–§4.8.
- **Assignment 4 due; Assignment 5** released: shot noise, matched filter,
  sampling.

### Week 11 — Chapter 5, Part 1

- **Lecture.** Wiener–Khinchin; PSD; Karhunen–Loève expansion;
  spectral factorization and the innovations representation.
- **Reading.** Chapter 5, §5.1–§5.3.

### Week 12 — Chapter 5, Part 2

- **Lecture.** Ergodicity; periodogram and Welch's method; LMMSE and
  orthogonality principle; Wiener filter; Kalman filter.
- **Tutorial 5.** Periodogram variance, LMMSE, causal Wiener, Kalman filter.
- **Reading.** Chapter 5, §5.4–§5.9.
- **Assignment 5 due; Assignment 6** released: LMMSE, Wiener, Kalman.

### Week 13 — Chapter 6, Part 1

- **Lecture.** Markov property; DTMC; transition matrix;
  Chapman–Kolmogorov; classification of states.
- **Reading.** Chapter 6, §6.1–§6.3.

### Week 14 — Chapter 6, Part 2

- **Lecture.** Stationary distribution; ergodicity; detailed balance;
  CTMC and generator matrix; birth–death chains.
- **Tutorial 6.** Classification; stationary distribution by power iteration
  and eigen-decomposition; CTMC.
- **Reading.** Chapter 6, §6.4–§6.7.
- **Assignment 6 due; Assignment 7** released: Markov-chain analysis.

### Week 15 — Chapter 7, Part 1

- **Lecture.** Kendall notation; Little's formula; birth–death framework;
  M/M/1 queue and performance measures; waiting-time distribution;
  Burke's theorem.
- **Reading.** Chapter 7, §7.1–§7.4.

### Week 16 — Chapter 7, Part 2

- **Lecture.** M/M/s queue; Erlang-C; M/M/1/K; PK formula;
  Jackson networks; simulation.
- **Tutorial 7.** Erlang-C staffing problem, M/M/1/K, discrete-event simulation.
- **Reading.** Chapter 7, §7.5–§7.9.
- **Assignment 7 due.**

### Week 17 — Revision

- **Lecture.** Review of Chapters 4–7; worked exam-style problems.
- **Tutorial.** Q&A session; past final-exam papers.

### Week 18 — Final examination

- **Exam.** 3 hours, closed book, **50% of final grade**. Covers the whole
  syllabus, weighted toward Chapters 4–7.

---

## Assessment timeline

| Week | Event | Weight |
|:---:|---|:---:|
| 2 | Assignment 1 released | — |
| 4 | Assignment 1 due · Assignment 2 released | 3% |
| 6 | Assignment 2 due · Assignment 3 released | 3% |
| 7 | Assignment 3 due | 3% |
| 8 | **Mid-semester examination** | **30%** |
| 9 | Assignment 4 released | — |
| 10 | Assignment 4 due · Assignment 5 released | 3% |
| 12 | Assignment 5 due · Assignment 6 released | 3% |
| 14 | Assignment 6 due · Assignment 7 released | 3% |
| 16 | Assignment 7 due | 2% |
| 18 | **Final examination** | **50%** |

*(Assignments total 20%; the per-assignment weights above are approximately equal.
The exact per-assignment weights will be announced when each assignment is released.)*

---

## Policies

!!! warning "Deadlines and attendance"
    - Late assignments lose **10% per calendar day**; after 3 days they score **zero**.
    - **Mid-semester** and **final** examinations are closed book;
      a non-programmable calculator is permitted.
    - Missing the mid-semester requires documented medical or family-emergency
      evidence; a make-up exam will be scheduled in **week 10**.
    - Attendance at tutorials is **strongly encouraged** — they are where most of
      the examinable worked examples are covered.

!!! tip "Holidays"
    National and religious holidays that fall on a lecture day will be announced
    on the [Home](index.md) page and, where possible, rescheduled. Check the
    announcements page each week.
