# Syllabus

Full course outline for **ECEG-6201 — Analytical Optimization and Stochastic
Process**. See the [Schedule](schedule.md) page for the week-by-week plan and
the [Assessment](assessment.md) page for grading details.

---

## Course identification

| Field | Value |
|---|---|
| **Course code** | ECEG-6201 |
| **Course title** | Analytical Optimization and Stochastic Process |
| **Degree program** | Master of Science in Communication Engineering |
| **Credit hours (ECTS)** | 3 (5) |
| **Semester** | Year 1, Semester I |
| **Status** | General Area Compulsory Course |
| **Pre-requisites** | None |
| **Instructor** | Dr. Tesfamichael Gebregziabher |
| **Contact** | *(office, email, office hours)* |

---

## Student workload

| Lecture | Tutorial | Practice / Laboratory | Home study |
|:---:|:---:|:---:|:---:|
| 2 h | 3 h | – | 5 h |

Total: ~10 hours per week × 18 weeks = **~180 hours** — consistent with 5 ECTS
(1 ECTS ≈ 30–36 hours of student effort).

---

## Course description

This course introduces the mathematical foundations and practical algorithms
of two complementary areas that underpin modern communication engineering:

- **Optimization** — how to formulate engineering design problems as mathematical
  programs and how to solve them numerically: unconstrained and constrained
  methods, gradient-based and derivative-free algorithms, penalty/barrier and
  augmented-Lagrangian approaches, quadratic programming.
- **Stochastic processes** — how to model random signals, noise, arrivals, and
  Markovian systems: stationarity and correlation, spectral representation,
  Karhunen–Loève expansion, innovations and factorization, Wiener and Kalman
  filtering, Markov chains, and queuing theory.

The two halves connect where they always do in practice: **estimation** (Chapter 5)
is optimization of a mean-square-error criterion under a stochastic model, and
**queuing** (Chapter 7) is a birth–death Markov chain whose stationary
distribution drives every performance measure.

---

## Course objectives and competences

By the end of the course, students will be able to:

1. **Familiarize** themselves with engineering optimization methods and know
   which to apply for a given problem.
2. **Formulate and solve** engineering design problems as optimization programs,
   choosing appropriate algorithms and justifying the choice.
3. **Understand the principles** of the modelling and analysis of random
   phenomena, distinguishing wide-sense from strict-sense stationarity and
   using correlation and spectral tools.
4. **Apply** the essential mathematical tools needed to analyze systems affected
   by random processes — including matched filtering, Wiener filtering, and
   Kalman filtering.
5. **Use** the algorithms required for the representation and estimation of
   random signals — Karhunen–Loève, spectral factorization, innovations
   representation, mean-square estimation.
6. **Analyze** Markov chains and elementary queuing systems, computing stationary
   distributions and performance measures.

---

## Learning outcomes

On successful completion, students will be able to:

| # | Outcome | Assessed by |
|---|---|---|
| LO1 | Characterize a minimum (local, global, Pareto), compute gradients and Hessians, apply first- and second-order optimality conditions | A1, Mid-sem |
| LO2 | Select and apply appropriate unconstrained algorithms (direct search, Newton, BFGS, Gauss–Newton, Levenberg–Marquardt) | A2, Mid-sem |
| LO3 | Formulate constrained problems and solve them via penalty, barrier, augmented-Lagrangian, and QP methods | A3, Mid-sem |
| LO4 | Analyse random processes: WSS/SSS, correlations, spectra, sampling, matched filter, entropy | A4, Final |
| LO5 | Apply spectral representation and estimation: KL, factorization, periodogram, LMMSE, Wiener, Kalman | A5, Final |
| LO6 | Classify and analyse Markov chains; compute stationary distributions; model as CTMC | A6, Final |
| LO7 | Analyse queuing systems: Little's formula, M/M/1, M/M/s, M/M/1/K, PK formula, Jackson networks | A7, Final |

---

## Teaching and learning methods

- **Lectures (2 h/week).** Motivate and derive each topic; work through at least
  one example per method; demonstrate the computational implementations on the
  laptop (Python is the primary language; Java and C# alternatives are provided
  in the notes).
- **Tutorials (3 h/week).** Small-group problem-solving sessions. The tutorial
  problems at the end of each chapter are worked on in class; solutions are
  posted to the [Solutions](solutions/index.md) page **after** the tutorial.
- **Assignments (7, one per chapter).** Longer problems combining theory and
  computation; graded with detailed feedback.
- **Independent study (5 h/week).** Reading, tutorial problems, assignment work,
  and running the computational examples.

---

## Course content

The syllabus comprises seven chapters, grouped into two halves.

### Part I — Optimization (Chapters 1–3)

**Chapter 1 — Optimization: Preliminaries**

- Characterization of a minimum; local vs. global minima; Pareto-optimality in multi-objective optimization
- Gradient and Hessian; slope and curvature
- Properties of linear and quadratic functions
- Optimality conditions: unconstrained; linearly constrained (equality, inequality); nonlinearly constrained (KKT)
- Line-search subproblem: bracketing, sectioning, golden-section
- Numerically derived gradients, analytically derived gradients, adjoint-variable approach

**Chapter 2 — Unconstrained Optimization**

- Non-derivative methods: polytope algorithm, Nelder–Mead simplex, pattern search (Hooke–Jeeves), AD-HOC methods
- Second-derivative methods: Newton, damped Newton
- First-derivative methods: quasi-Newton (BFGS), L-BFGS
- Sums of squares and nonlinear equations: Gauss–Newton, Levenberg–Marquardt

**Chapter 3 — Constrained Optimization**

- Penalty functions, barrier functions
- Augmented Lagrangian (method of multipliers)
- Quadratic programming with equality constraints (KKT system)
- QP with inequalities (active-set preview)

### Part II — Stochastic processes (Chapters 4–7)

**Chapter 4 — Stochastic Processes, Random Walks, and Applications**

- General concepts, definitions, correlations
- Random walks, Brownian-motion limit
- Poisson point processes; shot noise; Campbell's theorem
- Modulation and cyclostationary processes
- Band-limited processes and sampling theory
- Deterministic signal in noise; matched filter
- Shannon entropy and differential entropy

**Chapter 5 — Spectral Representation and Estimation**

- Power spectral density and Wiener–Khinchin
- Karhunen–Loève expansion
- Spectral factorization; innovations representation
- Ergodicity; periodogram; Welch's method
- Mean-square estimation and the orthogonality principle
- Wiener filter (non-causal and causal)
- Kalman filter; Riccati equation

**Chapter 6 — Markov Chains**

- Markov property; DTMC; transition matrix; Chapman–Kolmogorov
- Classification: recurrent/transient, periodic/aperiodic, communicating classes
- Stationary distribution; detailed balance; ergodic theorem
- Continuous-time Markov chains; generator matrix; Kolmogorov equations
- Birth–death chains

**Chapter 7 — Introduction to Queuing Theory**

- Kendall notation; elements of a queuing system
- Little's formula
- Birth–death process (CTMC)
- M/M/1: stationary distribution, performance measures, waiting-time distribution
- M/M/s: Erlang-C, performance measures
- M/M/1/K; Pollaczek–Khinchine (M/G/1); Jackson networks
- Simulation and data analysis of queuing systems

---

## Assessment overview

| Component | Weight | Timing |
|---|:---:|---|
| Assignments (7) | 20% | Weeks 4, 6, 7, 10, 12, 14, 16 |
| Mid-semester examination | 30% | Week 8 |
| Final examination | 50% | Week 18 |
| **Total** | **100%** | |

See the [Assessment](assessment.md) page for the grading rubric, submission
rules, late policy, and make-up examination procedure.

---

## Recommended textbook and references

**Primary references.**

- Boyd & Vandenberghe, *Convex Optimization*, Cambridge University Press, 2004.
- Papoulis & Pillai, *Probability, Random Variables, and Stochastic Processes*, 4th ed., McGraw-Hill, 2002.
- Nocedal & Wright, *Numerical Optimization*, 2nd ed., Springer, 2006.
- Gross et al., *Fundamentals of Queueing Theory*, 5th ed., Wiley, 2018.

**Secondary references.**

- Bertsekas, *Nonlinear Programming*, 3rd ed., Athena Scientific, 2016.
- Karlin & Taylor, *A First Course in Stochastic Processes*, 2nd ed., Academic Press, 1975.
- Ross, *Introduction to Probability Models*, 11th ed., Academic Press, 2014.
- Kay, *Fundamentals of Statistical Signal Processing: Estimation Theory*, Prentice-Hall, 1993.

Full references, free online versions, and tool recommendations are on the
[Resources](resources.md) page.

---

## Policies

- **Attendance.** Lectures and tutorials are not mandatory but strongly
  recommended — the tutorial sessions are where examinable worked examples
  are covered.
- **Academic integrity.** All submitted work must be your own. Collaboration on
  *understanding* is encouraged; collaboration on *submitted answers* is
  plagiarism. See the [Assessment](assessment.md) page for details.
- **Late work.** 10% per day, zero after 3 days.
- **Disability accommodations.** Contact the instructor at least one week
  before the first assessment to arrange accommodations.
- **Corrections.** To report a typo, an incorrect answer, or a broken link on
  this site, see the [Feedback](feedback.md) page.
