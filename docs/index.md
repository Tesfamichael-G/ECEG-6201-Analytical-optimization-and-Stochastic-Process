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

## What you will learn

<div class="grid cards" markdown>

- :material-target: **Formulate and solve optimization problems**

    ---

    Recognize an optimization problem, characterize its minima, and choose
    a gradient / Newton / quasi-Newton / penalty method appropriate to its
    smoothness, dimension, and cost per evaluation.

- :material-chart-bell-curve: **Model and analyze random phenomena**

    ---

    Work with random processes, correlations, spectra, and mean-square
    estimation. Apply stochastic-process tools to communication-engineering
    systems.

- :material-arrow-decision: **Analyze Markov chains and queues**

    ---

    Classify a Markov chain, compute its stationary distribution, and derive
    the performance measures of elementary queuing systems.

</div>

---

## Chapters

<div class="grid cards" markdown>

- :material-numeric-1-circle: **Optimization — Preliminaries**

    ---

    Minima, gradients, Hessians, optimality conditions (KKT), line search,
    numerical vs. analytic gradients, adjoint variables.

    [:octicons-arrow-right-24: Chapter 1](chapters/01-optimization-preliminaries.md)

- :material-numeric-2-circle: **Unconstrained Optimization**

    ---

    Direct search (Nelder–Mead, pattern search), Newton, quasi-Newton (BFGS),
    Gauss–Newton, Levenberg–Marquardt.

    [:octicons-arrow-right-24: Chapter 2](chapters/02-unconstrained-optimization.md)

- :material-numeric-3-circle: **Constrained Optimization**

    ---

    Penalty and barrier functions, augmented Lagrangian, quadratic programming
    (KKT system, active-set preview).

    [:octicons-arrow-right-24: Chapter 3](chapters/03-constrained-optimization.md)

- :material-numeric-4-circle: **Stochastic Processes & Random Walks**

    ---

    WSS vs. SSS, correlations, random walks, Poisson points, shot noise,
    cyclostationarity, sampling, matched filter, entropy.

    [:octicons-arrow-right-24: Chapter 4](chapters/04-stochastic-processes.md)

- :material-numeric-5-circle: **Spectral Representation & Estimation**

    ---

    Wiener–Khinchin, Karhunen–Loève, spectral factorization, ergodicity,
    periodogram/Welch, LMMSE, Wiener filter, Kalman filter.

    [:octicons-arrow-right-24: Chapter 5](chapters/05-spectral-representation.md)

- :material-numeric-6-circle: **Markov Chains**

    ---

    DTMC, classification (recurrent/transient/periodic), stationary distribution,
    ergodicity, CTMC and generator matrix.

    [:octicons-arrow-right-24: Chapter 6](chapters/06-markov-chains.md)

- :material-numeric-7-circle: **Introduction to Queuing Theory**

    ---

    Kendall notation, Little's formula, M/M/1, M/M/s (Erlang-C), M/M/1/K,
    Pollaczek–Khinchine, Jackson networks, discrete-event simulation.

    [:octicons-arrow-right-24: Chapter 7](chapters/07-queuing-theory.md)

</div>

---

## How to use this site

<div class="grid cards" markdown>

- :material-book-open-variant: **Read in order**

    ---

    Chapters are designed to be read in sequence — each builds on the previous.
    Chapter 1 lays the mathematical foundations for optimization; Chapter 4
    introduces the probabilistic tools that Chapters 5–7 build on.

- :material-check-circle-outline: **Check your understanding**

    ---

    Every chapter has a **three-part quiz** (MCQ, fill-in, define/explain)
    with an inline answer key. The **tutorial problems** at the end of each
    chapter are extra exercises; solutions are in the
    [Solutions](solutions.md) page.

- :material-laptop: **Run the code**

    ---

    Every chapter has Python, Java, and C# implementations of its main
    algorithms. Copy the snippets, run them, and compare against the analytical
    answers worked out by hand.

- :material-chart-bar: **See the schedule**

    ---

    The [Schedule](schedule.md) page has the week-by-week plan: lecture topics,
    tutorial sessions, reading assignments, and assessment dates.

</div>

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
Appendices A–C** (free PDF online) or **Papoulis & Pillai, Chapters 1–5** before
the first lecture. Both are linked from the [Resources](resources.md) page.

---

## Quick links

<div class="grid cards" markdown>

- :material-file-document: **[Syllabus](syllabus.md)**

    ---

    Course outline, prerequisites, learning outcomes, assessment weights.

- :material-calendar: **[Schedule](schedule.md)**

    ---

    Week-by-week topics, readings, assignment release/due dates, exam weeks.

- :material-book-multiple: **[Chapters](chapters/01-optimization-preliminaries.md)**

    ---

    Seven chapters of lecture notes, worked examples, code, and self-quizzes.

- :material-lightbulb-on: **[Solutions](solutions.md)**

    ---

    Key steps and final answers for the tutorial problems in every chapter.

- :material-clipboard-check: **[Assessment](assessment.md)**

    ---

    Grading rubric, submission rules, late policy, make-up exam procedure.

- :material-library: **[Resources](resources.md)**

    ---

    Textbooks, free online books and lecture notes, software, and benchmark
    datasets.

</div>
