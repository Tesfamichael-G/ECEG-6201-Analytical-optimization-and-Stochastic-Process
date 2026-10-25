# Chapter 1: Optimization - Preliminary

**Weeks:** 1 - 3 | **Contact Hours:** 6 Lectures, 9 Tutorials

## 🎯 Learning Objectives

By the end of this chapter, students will be able to:

1. Formulate engineering problems as mathematical optimization models.
2. Distinguish between local and global minima, and understand Pareto-optimality.
3. Compute and interpret the Gradient and Hessian matrices for linear and quadratic functions.
4. Apply First-Order and Second-Order Necessary/Sufficient conditions for unconstrained optimization.
5. Formulate and solve constrained optimization problems using the Karush-Kuhn-Tucker (KKT) conditions.
6. Implement line search algorithms and understand the Adjoint Variable method for gradient computation.

---

## 1.1 Basics of Optimization

An optimization problem seeks to find the best element from a set of available alternatives. Mathematically, it is formulated as:

$$\begin{aligned} \min_{\mathbf{x}} \quad & f(\mathbf{x}) \\ \text{subject to} \quad & \mathbf{x} \in \mathcal{X} \end{aligned}$$

where $\mathbf{x} \in \mathbb{R}^n$ is the decision variable vector, $f: \mathbb{R}^n \rightarrow \mathbb{R}$ is the objective function, and $\mathcal{X}$ is the feasible set.

### 1.1.1 Characterization of a Minimum
* **Local Minimum:** A point $\mathbf{x}^*$ is a local minimum if there exists an $\epsilon > 0$ such that $f(\mathbf{x}^*) \le f(\mathbf{x})$ for all $\mathbf{x}$ in the $\epsilon$-neighborhood of $\mathbf{x}^*$ (i.e., $\Vert{}\mathbf{x} - \mathbf{x}^*\Vert{} \le \epsilon$).
* **Global Minimum:** A point $\mathbf{x}^*$ is a global minimum if $f(\mathbf{x}^*) \le f(\mathbf{x})$ for *all* feasible $\mathbf{x} \in \mathcal{X}$.

> **Engineering Insight:** In communication engineering (e.g., beamforming or resource allocation), we often deal with non-convex problems where algorithms only guarantee finding a *local* minimum. Convexity is the key to ensuring a local minimum is also a global minimum.

### 1.1.2 Pareto-Optimality in Multi-Objective Optimization
When multiple conflicting objectives exist (e.g., maximizing data rate while minimizing power consumption), we define a vector objective $\mathbf{F}(\mathbf{x}) = [f_1(\mathbf{x}), f_2(\mathbf{x}), \dots, f_k(\mathbf{x})]^T$.

A solution $\mathbf{x}^*$ is **Pareto-optimal** (or non-dominated) if there is no other feasible $\mathbf{x}$ that improves one objective without worsening at least one other objective. The set of all Pareto-optimal solutions forms the **Pareto Front**.

### 1.1.3 Calculus Tools: Gradient and Hessian
To find minima, we rely on differential calculus to measure the **slope** and **curvature** of the objective function.

* **Gradient (Slope):** The vector of first partial derivatives.
  $$\nabla f(\mathbf{x}) = \begin{bmatrix} \frac{\partial f}{\partial x_1} \\ \vdots \\ \frac{\partial f}{\partial x_n} \end{bmatrix}$$
  It points in the direction of the steepest ascent.

* **Hessian (Curvature):** The square matrix of second partial derivatives.
  $$\nabla^2 f(\mathbf{x}) = \mathbf{H}(\mathbf{x}) = \begin{bmatrix} \frac{\partial^2 f}{\partial x_1^2} & \cdots & \frac{\partial^2 f}{\partial x_1 \partial x_n} \\ \vdots & \ddots & \vdots \\ \frac{\partial^2 f}{\partial x_n \partial x_1} & \cdots & \frac{\partial^2 f}{\partial x_n^2} \end{bmatrix}$$
  It describes the local curvature of the function.

### 1.1.4 Properties of Linear and Quadratic Functions
Understanding these two classes of functions is critical, as many optimization algorithms approximate complex functions using them.

1. **Linear Functions:** $f(\mathbf{x}) = \mathbf{a}^T \mathbf{x} + b$
   * Gradient: $\nabla f(\mathbf{x}) = \mathbf{a}$
   * Hessian: $\nabla^2 f(\mathbf{x}) = \mathbf{0}$ (Zero curvature; no local minima).
2. **Quadratic Functions:** $f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T \mathbf{Q} \mathbf{x} + \mathbf{c}^T \mathbf{x} + d$ (where $\mathbf{Q}$ is symmetric).
   * Gradient: $\nabla f(\mathbf{x}) = \mathbf{Q}\mathbf{x} + \mathbf{c}$
   * Hessian: $\nabla^2 f(\mathbf{x}) = \mathbf{Q}$ (Constant curvature).
   * *Note:* If $\mathbf{Q}$ is Positive Definite ($\mathbf{Q} \succ 0$), the function is strictly convex and has a unique global minimum.

---

## 1.2 Optimality Conditions

Optimality conditions are mathematical equations/inequalities that a point must satisfy to be a candidate for a minimum.

### 1.2.1 Unconstrained Optimization
For $\min f(\mathbf{x})$ with no constraints:
* **First-Order Necessary Condition (FONC):** If $\mathbf{x}^*$ is a local minimum, then $\nabla f(\mathbf{x}^*) = \mathbf{0}$. (These are called *stationary points*).
* **Second-Order Necessary Condition (SONC):** If $\mathbf{x}^*$ is a local minimum, then $\nabla^2 f(\mathbf{x}^*) \succeq 0$ (Positive Semi-Definite).
* **Second-Order Sufficient Condition (SOSC):** If $\nabla f(\mathbf{x}^*) = \mathbf{0}$ and $\nabla^2 f(\mathbf{x}^*) \succ 0$ (Positive Definite), then $\mathbf{x}^*$ is a strict local minimum.

### 1.2.2 Constrained Optimization
When constraints are present, the gradient of the objective function does not have to be zero; it only needs to be zero *in the directions that remain feasible*.

#### Linear Equality Constraints
$\min f(\mathbf{x})$ subject to $\mathbf{A}\mathbf{x} = \mathbf{b}$.
We introduce the **Lagrangian**:
$$\mathcal{L}(\mathbf{x}, \boldsymbol{\lambda}) = f(\mathbf{x}) - \boldsymbol{\lambda}^T (\mathbf{A}\mathbf{x} - \mathbf{b})$$
The necessary conditions are $\nabla_\mathbf{x} \mathcal{L} = \mathbf{0}$ and $\nabla_{\boldsymbol{\lambda}} \mathcal{L} = \mathbf{0}$.

#### Nonlinearly Constrained Optimization (KKT Conditions)
For $\min f(\mathbf{x})$ subject to $g_i(\mathbf{x}) \le 0$ and $h_j(\mathbf{x}) = 0$.
The **Karush-Kuhn-Tucker (KKT)** conditions are the necessary conditions for a local minimum:

1. **Stationarity:**
   $$\nabla f(\mathbf{x}^*) + \sum \mu_i \nabla g_i(\mathbf{x}^*) + \sum \lambda_j \nabla h_j(\mathbf{x}^*) = \mathbf{0}$$
2. **Primal Feasibility:** 
   $$g_i(\mathbf{x}^*) \le 0, \quad h_j(\mathbf{x}^*) = 0$$
3. **Dual Feasibility:** 
   $$\mu_i \ge 0$$
4. **Complementary Slackness:** 
   $$\mu_i g_i(\mathbf{x}^*) = 0$$

---

## 1.3 Linear Search Subproblem & Derivatives

In iterative optimization (like Gradient Descent or Newton's method), we update our estimate as:
$$\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{d}_k$$
where $\mathbf{d}_k$ is the search direction and $\alpha_k$ is the step size. Finding $\alpha_k$ is the **Linear Search** problem: $\min_{\alpha > 0} \phi(\alpha) = f(\mathbf{x}_k + \alpha \mathbf{d}_k)$.

### 1.3.1 Bracketing and Sectioning Phase
1. **Bracketing:** Find three points $a < b < c$ such that $\phi(b) < \phi(a)$ and $\phi(b) < \phi(c)$. The minimum must lie between $a$ and $c$.
2. **Sectioning (e.g., Golden Section Search):** Iteratively narrow the bracket by evaluating the function at specific interior points based on the golden ratio, discarding the section that cannot contain the minimum.

### 1.3.2 Derivatives of the Objective Function
Algorithms require gradients. How do we get them?

**1. Numerically Derived Gradient (Finite Differences):**
If the analytical derivative is too complex, we approximate it.
* *Forward difference:* $\frac{\partial f}{\partial x_i} \approx \frac{f(\mathbf{x} + h\mathbf{e}_i) - f(\mathbf{x})}{h}$
* *Central difference (more accurate):* $\frac{\partial f}{\partial x_i} \approx \frac{f(\mathbf{x} + h\mathbf{e}_i) - f(\mathbf{x} - h\mathbf{e}_i)}{2h}$

**2. Analytically Derived Gradient & The Adjoint Variable Approach:**
In communication systems (e.g., inverse scattering, channel estimation), the objective function often depends on a state variable $\mathbf{u}$ (like an electromagnetic field), which is governed by a physical constraint (like Maxwell's equations) $\mathbf{F}(\mathbf{u}, \mathbf{m}) = \mathbf{0}$, where $\mathbf{m}$ are the design parameters.

* **Direct differentiation** requires solving a linear system for *every single parameter* in $\mathbf{m}$. If $\mathbf{m}$ has 10,000 variables, this is computationally impossible.
* **The Adjoint Variable Approach:** By introducing an adjoint state variable $\boldsymbol{\lambda}$ (the Lagrange multiplier for the physical constraint), we can compute the exact analytical gradient with respect to $\mathbf{m}$ by solving only **one** adjoint equation, regardless of the number of parameters. This reduces the computational complexity from $O(N_{params})$ to $O(1)$.

---

## 📝 Worked Example: Constrained Optimization

**Problem:** Minimize the power $f(x_1, x_2) = x_1^2 + x_2^2$ subject to the linear constraint $x_1 + x_2 = 1$.

**Solution using Lagrange Multipliers:**
1. Form the Lagrangian: $\mathcal{L}(x_1, x_2, \lambda) = x_1^2 + x_2^2 - \lambda(x_1 + x_2 - 1)$
2. Apply Stationarity (take derivatives and set to 0):
   * $\frac{\partial \mathcal{L}}{\partial x_1} = 2x_1 - \lambda = 0 \implies x_1 = \lambda / 2$
   * $\frac{\partial \mathcal{L}}{\partial x_2} = 2x_2 - \lambda = 0 \implies x_2 = \lambda / 2$
   * $\frac{\partial \mathcal{L}}{\partial \lambda} = -(x_1 + x_2 - 1) = 0 \implies x_1 + x_2 = 1$
3. Substitute $x_1$ and $x_2$ into the constraint:
   $(\lambda / 2) + (\lambda / 2) = 1 \implies \lambda = 1$
4. Calculate optimal $\mathbf{x}^*$:
   $x_1^* = 0.5, \quad x_2^* = 0.5$

*Check Hessian:* $\nabla^2 f = \begin{bmatrix} 2 & 0 \\ 0 & 2 \end{bmatrix} \succ 0$. The solution is a strict global minimum.

---

## 💻 Computational Implementation (Python)

Here is how we implement gradient checking and constrained optimization using Python's `scipy.optimize`.

```python
import numpy as np
from scipy.optimize import minimize

# 1. Define the Objective Function (Quadratic)
def objective(x):
    return x[0]**2 + x[1]**2

# 2. Define the Analytical Gradient
def gradient(x):
    return np.array([2*x[0], 2*x[1]])

# 3. Numerical Gradient Check (Central Difference)
def numerical_gradient(f, x, h=1e-5):
    grad = np.zeros_like(x)
    for i in range(len(x)):
        x_plus = x.copy(); x_plus[i] += h
        x_minus = x.copy(); x_minus[i] -= h
        grad[i] = (f(x_plus) - f(x_minus)) / (2*h)
    return grad

# Test Gradient Check
x_test = np.array([1.0, 2.0])
print("Analytical Grad:", gradient(x_test))
print("Numerical Grad: ", numerical_gradient(objective, x_test))

# 4. Solve the Constrained Problem
# Constraint: x1 + x2 = 1  =>  x1 + x2 - 1 = 0
cons = ({'type': 'eq', 'fun': lambda x: x[0] + x[1] - 1})

# Initial guess
x0 = np.array([0.0, 0.0])

# Run Optimization (SLSQP handles equality constraints well)
result = minimize(objective, x0, method='SLSQP', jac=gradient, constraints=cons)

print("\nOptimization Result:")
print("Optimal x:", result.x)
print("Minimum f(x):", result.fun)
