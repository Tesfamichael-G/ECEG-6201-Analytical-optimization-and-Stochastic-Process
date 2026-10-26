# Resources

Textbooks, free online materials, software tools, and benchmark datasets for
**ECEG-6201**. All links open in a new tab; free resources are marked 🆓.

---

## Core textbooks

### Optimization (Chapters 1–3)

| Book | Chapters used | Notes |
|---|---|---|
| **Boyd & Vandenberghe**, *Convex Optimization*, Cambridge University Press, 2004 | 1–5, 9–11 | 🆓 PDF at [web.stanford.edu/~boyd/cvxbook](https://web.stanford.edu/~boyd/cvxbook/). Clearest treatment of optimality conditions and interior-point methods. |
| **Nocedal & Wright**, *Numerical Optimization*, 2nd ed., Springer, 2006 | 1–3, 6, 10, 16–17 | The standard reference for the algorithms in Chapters 1–3. |
| **Luenberger & Ye**, *Linear and Nonlinear Programming*, 4th ed., Springer, 2016 | 1–4 | Classical engineering treatment, lighter on convex analysis. |
| **Bertsekas**, *Nonlinear Programming*, 3rd ed., Athena Scientific, 2016 | 1, 3, 4 | The definitive reference on constrained optimization and duality. |
| **Bertsekas**, *Constrained Optimization and Lagrange Multiplier Methods*, Academic Press, 1982 | All | The original source for the augmented Lagrangian of Chapter 3. |
| **Fletcher**, *Practical Methods of Optimization*, 2nd ed., Wiley, 1987 | 3, 5, 10, 12 | Classic reference for QP active-set and Gauss–Newton methods. |

### Stochastic processes (Chapters 4–6)

| Book | Chapters used | Notes |
|---|---|---|
| **Papoulis & Pillai**, *Probability, Random Variables, and Stochastic Processes*, 4th ed., McGraw-Hill, 2002 | 8–11 | The standard stochastic-processes reference. Chapter 9 (spectra) is the source for Chapter 5. |
| **Leon-Garcia**, *Probability, Statistics, and Random Processes for Electrical Engineering*, 3rd ed., Pearson, 2008 | 5–9 | Gentle, engineering-oriented, with good examples. |
| **Karlin & Taylor**, *A First Course in Stochastic Processes*, 2nd ed., Academic Press, 1975 | 1–4, 12–14 | The classical reference for random walks, Poisson, and Markov chains (Chapters 4 and 6). |
| **Ross**, *Introduction to Probability Models*, 11th ed., Academic Press, 2014 | 2–6 | Excellent balance of theory and application; also the source for many queuing exercises. |
| **Norris**, *Markov Chains*, Cambridge University Press, 1997 | All | The definitive modern reference for Markov chains (Chapter 6). |
| **Levin, Peres & Wilmer**, *Markov Chains and Mixing Times*, AMS, 2009 | 1–4, 12 | Modern treatment with applications to MCMC. |
| **Kemeny & Snell**, *Finite Markov Chains*, Springer, 1976 | All | Matrix-theoretic treatment of finite chains. |
| **Kay**, *Fundamentals of Statistical Signal Processing: Estimation Theory*, Prentice-Hall, 1993 | 2–5, 12 | Rigorous treatment of LMMSE, Wiener, and Kalman. |
| **Hayes**, *Statistical Digital Signal Processing and Modeling*, Wiley, 1996 | 8, 10 | The definitive textbook for spectrum estimation (periodogram, Welch, AR, MUSIC). |
| **Stoica & Moses**, *Spectral Analysis of Signals*, Prentice-Hall, 2005 | 2–4 | Modern spectral analysis with rigorous treatment. |
| **Kailath, Sayed & Hassibi**, *Linear Estimation*, Prentice-Hall, 2000 | All | Unifies LMMSE, Wiener, and Kalman in one orthogonality framework. |
| **Anderson & Moore**, *Optimal Filtering*, Dover, 2005 | All | The classical Kalman and Wiener reference. |
| **Cover & Thomas**, *Elements of Information Theory*, 2nd ed., Wiley, 2006 | 2–3, 8–10 | Rigorous treatment of Shannon entropy, differential entropy, and the Gaussian channel. |

### Queuing theory (Chapter 7)

| Book | Chapters used | Notes |
|---|---|---|
| **Gross, Shortle, Thompson & Harris**, *Fundamentals of Queueing Theory*, 5th ed., Wiley, 2018 | 1–5 | The standard textbook. Chapter 5 covers the Pollaczek–Khinchine formula and M/G/1. |
| **Kleinrock**, *Queueing Systems, Volume 1: Theory*, Wiley, 1975 | 1–5 | The classical reference, still worth reading. |
| **Wolff**, *Stochastic Modeling and the Theory of Queues*, Prentice-Hall, 1989 | 1–6 | Rigorous treatment including Little's formula and Burke's theorem. |
| **Harchol-Balter**, *Performance Modeling and Design of Computer Systems*, Cambridge University Press, 2013 | 1–20 | Modern treatment with computer-science applications; excellent exercises. |
| **Bertsekas & Gallager**, *Data Networks*, 2nd ed., Prentice-Hall, 1992 | 3, 5 | Chapters 3 and 5 for queuing in communication networks. |
| **Medhi**, *Stochastic Models in Queueing Theory*, 2nd ed., Academic Press, 2003 | 1–4 | Comprehensive single-volume reference. |

---

## Free online textbooks and lecture notes 🆓

- **Boyd & Vandenberghe**, *Convex Optimization* — full PDF at
  [web.stanford.edu/~boyd/cvxbook](https://web.stanford.edu/~boyd/cvxbook/).
  Companion course: **Stanford EE364a** (lectures on YouTube).
- **Grinstead & Snell**, *Introduction to Probability* — full PDF at
  [math.dartmouth.edu/~prob/prob/prob.pdf](https://math.dartmouth.edu/~prob/prob/prob.pdf).
  Excellent for reviewing probability before Chapter 4.
- **MIT 6.041SC**, *Probabilistic Systems Analysis and Applied Probability* —
  [ocw.mit.edu/courses/6-041sc](https://ocw.mit.edu/courses/6-041sc-probabilistic-systems-analysis-and-applied-probability-fall-2013/).
  Lecture notes, problem sets, and video lectures.
- **MIT 6.262**, *Discrete Stochastic Processes* (Robert Gallager) —
  [ocw.mit.edu/courses/6-262](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/).
  Graduate-level stochastic processes; the closest match to Chapters 4–7.
- **Bertsekas**, free online books at
  [web.mit.edu/dimitrib/www/home.html](https://web.mit.edu/dimitrib/www/home.html).
  Includes *Convex Analysis and Optimization* and *Dynamic Programming and Optimal Control, Vol. II*.
- **Boyd**, *EE364a lecture slides and video* —
  [see.stanford.edu/see/courseinfo.aspx?coll=2d7120e8-8a46-4b7c-a9a0-25a05b1c0c9c](https://see.stanford.edu/see/courseinfo.aspx?coll=2d7120e8-8a46-4b7c-a9a0-25a05b1c0c9c)
  and the YouTube channel **Stanford Online**.
- **Nocedal & Wright**, *Numerical Optimization* — many chapters have
  freely available PDF preprints on the authors' homepages (Jorge Nocedal,
  Northwestern University; Stephen Wright, University of Wisconsin–Madison).

---

## Online courses and video lectures

| Course | Provider | Link |
|---|---|---|
| Convex Optimization (EE364a) | Stanford | [see.stanford.edu](https://see.stanford.edu/) |
| Probabilistic Systems Analysis (6.041SC) | MIT OCW | [ocw.mit.edu](https://ocw.mit.edu/courses/6-041sc-probabilistic-systems-analysis-and-applied-probability-fall-2013/) |
| Discrete Stochastic Processes (6.262) | MIT OCW | [ocw.mit.edu](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/) |
| Probabilistic Graphical Models | Stanford (Coursera) | [coursera.org](https://www.coursera.org/learn/probabilistic-graphical-models) |
| Matrix Methods in Data Analysis, Signal Processing, and Machine Learning (18.065) | MIT OCW | [ocw.mit.edu](https://ocw.mit.edu/courses/18-065-matrix-methods-in-data-analysis-signal-processing-and-machine-learning-spring-2018/) |
| Markov Chains and Mixing Times (short course) | Various | Search for the *Levin–Peres–Wilmer* companion materials |

---

## Software

### Python (primary language of the course)

| Package | Use |
|---|---|
| `numpy` | Arrays, linear algebra, `linalg.solve`, `linalg.eig`, `linalg.matrix_power` |
| `scipy.optimize` | `minimize` (Nelder–Mead, BFGS, Newton-CG, SLSQP, trust-constr), `least_squares` (LM, TRF), `minimize_scalar` (golden), `check_grad` |
| `scipy.linalg` | `expm` (matrix exponential for CTMC), `levinson` (Levinson–Durbin) |
| `scipy.signal` | `periodogram`, `welch`, `wiener` |
| `statsmodels.tsa` | Kalman filter, AR estimation |
| `simpy` | Discrete-event simulation framework |
| `queueing_tool` | Queueing-specific simulation |
| `cvxpy` | Convex optimization modeling (LP, QP, SOCP, SDP) |
| `filterpy`, `pykalman` | Kalman filter libraries |
| `markovchain`, `PyDTMC` | Markov-chain analysis |
| `matplotlib`, `seaborn` | Plotting |
| `jupyter` | Notebooks (recommended for the worked examples) |

Install everything with:

```bash
pip install numpy scipy matplotlib jupyter statsmodels simpy cvxpy filterpy pykalman
```

### Java

| Library | Use |
|---|---|
| **Apache Commons Math** | `Optimization` (BFGS, Nelder–Mead, LM), `KalmanFilter`, `RealMatrix` |
| **Hipparchus** | Modern successor to Commons Math (same API, bug fixes) |
| **JTransforms** | Fast FFT for spectral estimation |
| **SSJ** (*Stochastic Simulation in Java*, McGill) | Queueing, random-variable generation, DES |
| **SimJava** | General DES framework |
| **ojalgo** | QP, LP solvers |

### C#

| Library | Use |
|---|---|
| **Math.NET Numerics** | Linear algebra, FFT, distributions, `NelderMeadSimplex`, `BfgsMinimizer` |
| **SimSharp** | SimPy-inspired DES for .NET |
| **Google OR-Tools** | LP, MIP, QP solvers (wrappers for CPLEX / Gurobi / SCIP) |
| **ScottPlot** | Plotting |

### MATLAB / R / Julia

| Language | Toolbox |
|---|---|
| **MATLAB** | `fminunc` (unconstrained), `fmincon` (constrained), `quadprog` (QP), `pwelch`, `dsp.KalmanFilter`, `simEvents` |
| **R** | `optim`, `nloptr`, `queueing`, `markovchain` |
| **Julia** | `Optim.jl`, `JuMP`, `Distributions.jl`, `LinearAlgebra`, `Plots.jl` |

---

## Benchmark problems and datasets

### Optimization benchmarks

| Benchmark | Dimensions | Link |
|---|---|---|
| **Rosenbrock** | $n \ge 2$ | [Wikipedia](https://en.wikipedia.org/wiki/Rosenbrock_function) |
| **Rastrigin** | $n \ge 2$ | [Wikipedia](https://en.wikipedia.org/wiki/Rastrigin_function) |
| **Ackley** | $n \ge 2$ | [Wikipedia](https://en.wikipedia.org/wiki/Ackley_function) |
| **NIST StRD nonlinear regression** | Various | [itl.nist.gov/div898/strd/nls/nls_main.shtml](https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml) — 27 certified nonlinear least-squares problems |
| **NETLIB LP test set** | up to $10^5$ | [netlib.org/lp/data/](http://www.netlib.org/lp/data/) |
| **MIPLIB** | mixed-integer | [miplib.zib.de](https://miplib.zib.de/) |
| **CUTEst** | general nonlinear | [cuter.rl.ac.uk](https://www.cuter.rl.ac.uk/) |

### Stochastic-processes datasets

| Dataset | Use | Link |
|---|---|---|
| **CAIDA traces** | Network traffic, packet interarrival times | [caida.org](https://www.caida.org/) |
| **WIDE backbone traces** | Same | [wide.ad.jp](https://www.wide.ad.jp/) |
| **NASA HTTP logs** | Web-server arrivals | [ita.ee.lbl.gov/html/contrib/NASA-HTTP.html](http://ita.ee.lbl.gov/html/contrib/NASA-HTTP.html) |
| **Calgary HTTP logs** | Same | [ita.ee.lbl.gov/html/contrib/Calgary-HTTP.html](http://ita.ee.lbl.gov/html/contrib/Calgary-HTTP.html) |
| **UCI ML repository** | Time-series, sensor data | [archive.ics.uci.edu/ml](https://archive.ics.uci.edu/ml) |
| **PhysioNet** | Biomedical time series | [physionet.org](https://physionet.org/) |

---

## Mathematical references

| Resource | Use |
|---|---|
| **DLMF** — *Digital Library of Mathematical Functions* | [dlmf.nist.gov](https://dlmf.nist.gov/) — definitive reference for special functions |
| **MathWorld** | [mathworld.wolfram.com](https://mathworld.wolfram.com/) — concise definitions |
| **Wikipedia** | [en.wikipedia.org](https://en.wikipedia.org/) — good starting points for Rosenbrock, Rastrigin, KKT, Markov chains, M/M/1, etc. |
| **The Matrix Cookbook** | [math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf](https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf) — matrix identities, derivatives, useful for Chapter 1 |
| **Numerical Recipes** | [numerical.recipes](http://numerical.recipes/) — free 3rd edition PDF online; useful for the numerical algorithms in Chapters 1–2 |

---

## Getting started — installation quick reference

### Python (recommended)

```bash
# 1. Install Python 3.10 or later from python.org
# 2. Create a virtual environment (recommended)
python -m venv eceg-env
source eceg-env/bin/activate      # Windows: eceg-env\Scripts\activate

# 3. Install the core packages
pip install numpy scipy matplotlib jupyter statsmodels

# 4. Install optional packages as needed
pip install simpy cvxpy filterpy pykalman
```

### Java

```bash
# 1. Install JDK 17 or later (Adoptium, Amazon Corretto, or OpenJDK)
# 2. Add Apache Commons Math to your project
#    Maven: <dependency>
#             <groupId>org.apache.commons</groupId>
#             <artifactId>commons-math3</artifactId>
#             <version>3.6.1</version>
#           </dependency>
# 3. (Optional) Add Hipparchus for the newer API
```

### C#

```bash
# 1. Install .NET 8 SDK from dotnet.microsoft.com
# 2. Add Math.NET Numerics
dotnet add package MathNet.Numerics
# 3. (Optional) Add SimSharp for discrete-event simulation
dotnet add package SimSharp
```

### MATLAB / R / Julia

| Language | Installation |
|---|---|
| MATLAB | Requires a license (contact your university's IT services) |
| R | [cran.r-project.org](https://cran.r-project.org/); install packages with `install.packages(c("optim", "queueing"))` |
| Julia | [julialang.org](https://julialang.org/); install packages with `] add Optim, JuMP, Distributions` |

---

## Recommended first-week reading

Before the first lecture, read:

- **Boyd & Vandenberghe, §1.1–§1.4** — for the optimization view of Chapter 1.
- **Papoulis & Pillai, Chapter 8** — for the stochastic-process view of Chapter 4.
- **Ross, Chapter 2** — for the probability review.

Do the exercises in the first chapter of each book — they are the fastest way to
check that your background is sufficient.
