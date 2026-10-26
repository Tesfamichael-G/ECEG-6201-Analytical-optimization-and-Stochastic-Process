# Welcome to ECEG-6201

## Analytical Optimization and Stochastic Process

**Degree Program:** Master of Science in Communication Engineering  
**Semester:** Year 1, Semester I | **Status:** General Area Compulsory Course  
**ECTS Credits:** 5 (3 Credit Hours)  

Welcome to the official course website for **ECEG-6201**. This course bridges the gap between deterministic engineering optimization and the analysis of random phenomena in communication systems.

---

## 🎯 Course Overview

=== "Optimization (Ch 1-3)"

    * Familiarize students with engineering optimization methods.
    * Equip students with different techniques for formulating and solving engineering problems.
    * Master unconstrained and constrained optimization (KKT conditions, Penalty/Barrier functions).

=== "Stochastic Processes (Ch 4-7)"

    * Understand the principles employed in the modelling and analysis of random phenomena.
    * Analyze systems affected by random processes and signals.
    * Master spectral representation, Wiener/Kalman filters, Markov Chains, and Queuing Theory.

---

## 📋 Course Information

| **Pre-requisites** | **Workload (Hours/Week)** | **Teaching Methods** |
| :--- | :--- | :--- |
| None | **Lecture:** 2 <br> **Tutorial:** 3 <br> **Home Study:** 5 | Lecture supported by tutorial and assignments. |

### 📊 Grading & Assessment

!!! warning "Evaluation Scheme"
    * **Assignments:** 20%
    * **Mid-semester Examination:** 30%
    * **Final Examination:** 50%

---

## 📅 16-Week Course Schedule

| Week | Chapter | Lecture Topics | Tutorial / Practice Activities | Deliverables |
| :---: | :--- | :--- | :--- | :--- |
| **1** | Ch 1 | **Optimization Preliminaries:** Characterization of minima, Local vs. Global, Pareto-optimality, Gradient & Hessian. | Review of multivariable calculus. Manual calculation of Gradients and Hessians. | |
| **2** | Ch 1 | **Optimality Conditions:** Unconstrained, Linear equality/inequality, Nonlinearly constrained (KKT). | Solving optimality condition problems. Identifying active constraints. | |
| **3** | Ch 1 & 2 | **Linear Search & Non-Derivative Methods:** Line search algorithms, Direct search: Polytope, Simplex. | Implementing 1D line search (Golden Section) and Nelder-Mead simplex in Python/MATLAB. | |
| **4** | Ch 2 | **Derivative-Based Unconstrained:** Newton’s method, Quasi-Newton methods (BFGS, DFP). | Step-by-step manual iteration of Newton's method. Comparing convergence rates. | **Assignment 1 Released** |
| **5** | Ch 2 & 3 | **Nonlinear Equations & Constrained Intro:** Gauss-Newton, Levenberg-Marquardt. Penalty/Barrier functions. | Solving nonlinear least squares. Formulating penalty functions. | |
| **6** | Ch 3 | **Advanced Constrained:** Augmented Lagrangian, Quadratic Programming (QP). | Formulating engineering problems as QPs. Solving QPs using `cvxpy`. | **Assignment 1 Due** |
| **7** | **Exam** | **Mid-Semester Examination** (Covers Chapters 1, 2, and 3). | Exam review session, Q&A, and post-exam discussion. | **Mid-Semester Exam (30%)** |
| **8** | Ch 4 | **Stochastic Processes Basics:** General concepts, Correlations, Random Walks, Poisson points. | Probability refresher. Simulating and visualizing 1D/2D Random Walks. | |
| **9** | Ch 4 | **Advanced Stochastic:** Modulation, Cyclostationary Processes, Sampling theory, Entropy. | Calculating autocorrelation and PSD. Entropy calculations. | **Assignment 2 Released** |
| **10** | Ch 5 | **Spectral Representation:** Factorization, Innovations, Karhunen-Loève (KL) expansions. | Working through KL expansion examples for simple processes. | |
| **11** | Ch 5 | **Estimation Theory:** Ergodicity, Mean square estimation, Wiener filter. | Designing a discrete-time Wiener filter for a given SNR scenario. | **Assignment 2 Due** |
| **12** | Ch 5 & 6 | **Kalman Filtering & Markov Chains Intro:** Kalman filter equations. Discrete/Continuous Markov chains. | Step-by-step prediction and update of a simple Kalman filter. | |
| **13** | Ch 6 & 7 | **Markov Properties & Queuing Basics:** Ergodic properties. Little’s formula, Birth-death processes. | Solving steady-state probabilities. Applying Little's Law ($L = \lambda W$). | **Assignment 3 Released** |
| **14** | Ch 7 | **Advanced Queuing & Simulation:** M/M/1 and M/M/s systems. Simulation and data analysis. | Coding an M/M/1 queue simulator. Analyzing simulated data. | **Assignment 3 Due** |
| **15** | Review | **Comprehensive Course Review:** Connecting Optimization with Stochastic Estimation. | Mock Final Exam, open problem-solving, and Q&A session. | |
| **16** | **Exam** | **Final Examination Period** (Covers all chapters, emphasis on Ch 4-7). | N/A | **Final Exam (50%)** |

---

## 📚 Quick Navigation

*   [📖 **Chapter 1:** Optimization Preliminaries](chapter_01.md)
*   *(More chapters will be added weekly as the course progresses)*
