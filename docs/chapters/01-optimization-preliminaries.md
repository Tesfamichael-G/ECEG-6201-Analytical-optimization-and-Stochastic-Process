
# Chapter 1 — Optimization: Preliminaries

!!! abstract "Learning outcomes"
    By the end of this chapter you should be able to:
    
    - distinguish local, global, and Pareto-optimal solutions;
    - compute gradients and Hessians and interpret them geometrically;
    - state first- and second-order optimality conditions;
    - set up and solve a line-search subproblem.

## 1.1 Characterization of a minimum

A point \(x^\star\) is a **local minimizer** of \(f\) if there exists \(\epsilon > 0\) such that
\(f(x^\star) \le f(x)\) for all \(\|x - x^\star\| < \epsilon\).
If the inequality holds for all \(x\) in the feasible set, \(x^\star\) is a **global minimizer**.

## 1.2 Optimality conditions (unconstrained)

First-order necessary condition:

\[
\nabla f(x^\star) = 0
\]

Second-order necessary condition: \(\nabla^2 f(x^\star) \succeq 0\) (positive semidefinite).

## 1.3 Pareto-optimality

For multi-objective problems \(\min_x \big(f_1(x), \dots, f_m(x)\big)\), a point \(x^\star\) is **Pareto-optimal**
if no feasible \(x\) improves one objective without worsening another.

## 1.4 Line search

Given a descent direction \(d_k\), find \(\alpha_k\) solving

\[
\min_{\alpha > 0} \phi(\alpha) = f(x_k + \alpha d_k).
\]

Bracketing establishes an interval containing a local minimum; sectioning (golden-section,
bisection, or interpolation) then refines it.

!!! example "Exercise"
    Derive the adjoint-variable formula for the gradient of
    \(J(u) = \int_0^T g(x(t), u(t))\,dt\) subject to \(\dot x = F(x, u)\).
