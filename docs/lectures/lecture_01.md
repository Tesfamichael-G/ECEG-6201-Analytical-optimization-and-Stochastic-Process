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

$$
\begin{aligned}
\min_{\mathbf{x}} \quad & f(\mathbf{x}) \\
\text{subject to} \quad & \mathbf{x} \in \mathcal{X}
\end{aligned}
$$

where $\mathbf{x} \in \mathbb{R}^n$ is the decision variable vector, $f: \mathbb{R}^n \rightarrow \mathbb{R}$ is the objective function, and $\mathcal{X}$ is the feasible set.

### 1.1.1 Characterization of a Minimum
* **Local Minimum:** A point $\mathbf{x}^*$ is a local minimum if there exists an $\epsilon > 0$ such that $f(\mathbf{x}^*) \le f(\mathbf{x})$ for all $\mathbf{x}$ in the $\epsilon$-neighborhood of $\mathbf{x}^*$ (i.e., $\|\mathbf{x} - \mathbf{x}^*\| \le \epsilon$).
* **Global Minimum:** A point $\mathbf{x}^*$ is a global minimum if $f(\mathbf{x}^*) \le f(\mathbf{x})$ for *all* feasible $\mathbf{x} \in \mathcal{X}$.

### 1.1.2 Pareto-Optimality in Multi-Objective Optimization
When multiple conflicting objectives exist, we define a vector objective $\mathbf{F}(\mathbf{x}) = [f_1(\mathbf{x}), f_2(\mathbf{x}), \dots, f_k(\mathbf{x})]^T$.

A solution $\mathbf{x}^*$ is **Pareto-optimal** if there is no other feasible $\mathbf{x}$ that improves one objective without worsening at least one other.

### 1.1.3 Calculus Tools: Gradient and Hessian
* **Gradient (Slope):**
  $$ \nabla f(\mathbf{x}) = \begin{bmatrix} \frac{\partial f}{\partial x_1} \\ \vdots \\ \frac{\partial f}{\partial x_n} \end{bmatrix} $$

* **Hessian (Curvature):**
  $$ \nabla^2 f(\mathbf{x}) = \mathbf{H}(\mathbf{x}) = \begin{bmatrix} \frac{\partial^2 f}{\partial x_1^2} & \cdots & \frac{\partial^2 f}{\partial x_1 \partial x_n} \\ \vdots & \ddots & \vdots \\ \frac{\partial^2 f}{\partial x_n \partial x_1} & \cdots & \frac{\partial^2 f}{\partial x_n^2} \end{bmatrix} $$
