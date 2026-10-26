# ECEG-6201 — Analytical Optimization and Stochastic Processes

**MSc in Communication Engineering · Year 1, Semester I**

Welcome to the course website. This site hosts the complete lecture notes,
worked examples, computational implementations (Python / Java / C#), quizzes,
and tutorial problems for the two halves of this course: **optimization**
(Chapters 1–3) and **stochastic processes** (Chapters 4–7).

---

## Course at a glance

| | |
|---|---|
| **Course code** | ECEG-6201 |
| **Credit hours (ECTS)** | 3 (5) |
| **Semester** | Year 1, Semester I |
| **Status** | General Area Compulsory Course |
| **Pre-requisites** | None |
| **Instructor** | Dr. Tesfamichael Gebregziabher |
| **Contact** | *(office, email, office hours)* |

**Student workload (hours per week)**

| Lecture | Tutorial | Practice / Laboratory | Home study |
|:---:|:---:|:---:|:---:|
| 2 | 3 | – | 5 |

**Teaching & learning methods.** Lectures supported by tutorials and assignments.

---

## Learning objectives

By the end of this course you should be able to:

1. **Formulate and solve** unconstrained and constrained optimization problems using
   gradient, Newton, quasi-Newton, and penalty/barrier methods.
2. **Analyze** engineering problems as optimization problems and select an appropriate
   numerical method by considering smoothness, dimension, and cost per evaluation.
3. **Model and analyze** random phenomena using probability, autocorrelation,
   spectral representation, and mean-square estimation.
4. **Apply** stochastic-process tools (random walks, Poisson processes, shot noise,
   cyclostationarity, Wiener and Kalman filters) to communication-engineering systems.
5. **Analyze** Markov chains and elementary queuing systems and compute their
   stationary distributions and performance measures.

---

## Chapter overview

| Ch | Title | Focus |
|:---:|---|---|
| **1** | [Optimization — Preliminaries](chapters/01-optimization-preliminaries.md) | minima, gradients, Hessians, optimality conditions, line search, derivatives |
| **2** | [Unconstrained Optimization](chapters/02-unconstrained-optimization.md) | direct search, Newton, quasi-Newton (BFGS), Gauss–Newton, Levenberg–Marquardt |
| **3** | [Constrained Optimization](chapters/03-constrained-optimization.md) | penalty, barrier, augmented Lagrangian, quadratic programming |
| **4** | [Stochastic Processes & Random Walks](chapters/04-stochastic-processes.md) | random processes, correlations, random walks, Poisson points, shot noise, sampling, entropy |
| **5** | [Spectral Representation & Estimation](chapters/05-spectral-representation.md) | PSD, Karhunen–Loève, spectral factorization, ergodicity, Wiener filter, Kalman filter |
| **6** | [Markov Chains](chapters/06-markov-chains.md) | DTMC, CTMC, classification, stationarity, ergodicity |
| **7** | [Introduction to Queuing Theory](chapters/07-queuing-theory.md) | Kendall notation, Little's formula, M/M/1, M/M/s, Erlang-C, Jackson networks, simulation |

Each chapter contains:

- **Learning outcomes** and a chapter-at-a-glance summary
- **Numbered worked examples** with full derivations
- **Computational implementation** in **Python, Java, and C#**
- **Chapter summary** and **further reading**
- A **three-part quiz** — multiple choice, fill-in-the-blank, define/explain —
  with a hidden **answer key**
- **Suggested tutorial problems** for self-study

---

## Assessment

| Component | Weight |
|---|:---:|
| Assignments (7, one per chapter) | **20%** |
| Mid-semester examination (week 8) | **30%** |
| Final examination | **50%** |
| **Total** | **100%** |

See the [Assessment](assessment.md) page for submission rules, late policy,
and detailed grading rubric. For the week-by-week plan, see the
[Schedule](schedule.md).

---

## Recommended preparation

No prior background in optimization or stochastic processes is assumed,
but the following will make the course much easier:

- **Linear algebra** — matrix algebra, eigenvalues, positive-definite matrices
- **Calculus** — partial derivatives, gradients, Taylor expansion
- **Probability** — random variables, expectation, variance, common distributions
- **Python** — NumPy and SciPy are used throughout for the computational demos
  (Java and C# equivalents are also provided for students from those backgrounds)

If any of these feels shaky, review **Boyd & Vandenberghe, *Convex Optimization*,
Appendix A–C** (free PDF online) or **Papoulis & Pillai, Chapters 1–5** before
the first lecture.

---

## Tools and references

| Category | Recommended tools |
|---|---|
| **Optimization** | `scipy.optimize`, `cvxpy`, `ipopt`, `MOSEK` |
| **Stochastic processes** | `numpy.random`, `scipy.stats`, `statsmodels.tsa` |
| **Queuing theory** | `simpy`, `queueing_tool`, `SSJ` (Java), `SimSharp` (C#) |
| **Plotting** | `matplotlib`, `seaborn` |

Full references and textbook recommendations are on the
[Resources](resources.md) page.

---

## How to use this site

!!! info "Reading order"
    Chapters are designed to be read **in order** — each builds on the previous.
    Chapter 1 lays the mathematical foundations for optimization; Chapter 4
    introduces the probabilistic tools that Chapters 5–7 build on. Skipping
    ahead is possible but not recommended for a first pass.

!!! tip "Self-assessment"
    Every chapter has an **answer key** for its quizzes, so you can check your
    understanding immediately after working through the problems. The
    **tutorial problems** at the end of each chapter are designed as additional
    exercises — solutions are not published, but they are discussed in the
    weekly tutorial session.

!!! note "Corrections and feedback"
    For corrections, suggestions, or bug reports, please open an issue on the
    [GitHub repository](https://github.com/tesfamichael-g/ECEG-6201-Analytical-optimization-and-Stochastic-Process)
    or email the instructor.

---

## Quick links

| Link | Description |
|---|---|
| [Syllabus](syllabus.md) | Course outline, prerequisites, learning outcomes |
| [Schedule](schedule.md) | Week-by-week topics, readings, and assessment dates |
| [Chapters](chapters/01-optimization-preliminaries.md) | Lecture notes, worked examples, code |
| [Assessment](assessment.md) | Grading rubric, submission rules, late policy |
| [Resources](resources.md) | Textbooks, tools, and datasets |
