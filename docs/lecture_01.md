---
sidebar_position: 1
---

# Chapter 1: Optimization - Preliminary

## 🎯 Learning Objectives

By the end of this chapter, students will be able to:
1. Formulate engineering problems as mathematical optimization models.
2. Distinguish between local and global minima, and understand Pareto-optimality.
3. Compute and interpret the Gradient and Hessian matrices.

## 1.1 Basics of Optimization

An optimization problem seeks to find the best element from a set of available alternatives:

$$
\begin{aligned}
\min_{\mathbf{x}} \quad & f(\mathbf{x}) \\
\text{subject to} \quad & \mathbf{x} \in \mathcal{X}
\end{aligned}
$$

### 1.1.1 Characterization of a Minimum
* **Local Minimum:** A point $\mathbf{x}^*$ is a local minimum if there exists an $\epsilon > 0$ such that $f(\mathbf{x}^*) \le f(\mathbf{x})$ for all $\mathbf{x}$ in the $\epsilon$-neighborhood.
* **Global Minimum:** A point $\mathbf{x}^*$ is a global minimum if $f(\mathbf{x}^*) \le f(\mathbf{x})$ for *all* feasible $\mathbf{x} \in \mathcal{X}$.

## 💻 Computational Implementation (Python)

```python
import numpy as np
from scipy.optimize import minimize

def objective(x):
    return x[0]**2 + x[1]**2

def gradient(x):
    return np.array([2*x[0], 2*x[1]])

cons = ({'type': 'eq', 'fun': lambda x: x[0] + x[1] - 1})
x0 = np.array([0.0, 0.0])

result = minimize(objective, x0, method='SLSQP', jac=gradient, constraints=cons)
print("Optimal x:", result.x)
