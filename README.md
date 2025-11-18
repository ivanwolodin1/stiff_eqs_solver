# Solving Stiff Differential Equations with `py-pde`

This repository shows how to solve **stiff differential equations** using the  
[`py-pde`](https://py-pde.readthedocs.io/en/stable/) package (v.0.47.2).  
Two examples (so far) are included:

1. A **scalar stiff ODE** from *Hairer & Wanner (1996)*  
2. **Simple heat equation**
2. A **Custom 1D stiff PDE system** 

---

# Stiff ODE Example

We want to solve a stiff ODE:

$$\frac{dy}{dx} = -50(y - \cos x)$$

This equation appears in:

> Hairer, E. & Wanner, G.  
> *Solving Ordinary Differential Equations II. Stiff and Differential-Algebraic Problems.*  
> Springer, 1996.
  
Because of stiffness, explicit solvers require **very small time steps**.
This `py-pde` package is a handy programming interface that allows to solve these eqs. 

---

For instance:

```python
from pde import CartesianGrid, ScalarField, PDE

grid = CartesianGrid([[0, 1.5]], 100)
state = ScalarField(grid, data=0.0)

eq = PDE({'y': '-50 * y + 50 * cos(x)'})

result = eq.solve(
    state,
    t_range=[0, 20],
    dt=0.1,
    solver='scipy',
)
```

# System of stiff DE  
Th utimate goal is to solve the stiff system of eqs that we obtained and published in MS&T:

**Volodin, I., Alabuzhev, A.** *Linear Stability of Marangoni Convection in a Thin Film under Vertical Vibrations.* Microgravity Sci. Technol. 37, 48 (2025). [https://doi.org/10.1007/s12217-025-10201-8](https://doi.org/10.1007/s12217-025-10201-8)

It has the following form:

$$
\begin{align}
    \frac{\partial \bar{\xi}}{\partial \bar{t}}
    &= \frac{\partial}{\partial X} \left(
    \frac{\bar{\xi}^3}{3} \frac{\partial \bar{\Pi}}{\partial X}
    + \frac{\text{Ma}}{2} \bar{\xi}^2 \frac{\partial f}{\partial X} \right. \notag \\
    &\quad \left. - \frac{(b \Omega)^2}{2} \left(
    Q_1 \bar{\xi}^2 \left( \frac{\partial \bar{\xi}}{\partial X} \right)^3
    + Q_2 \bar{\xi}^2 \frac{\partial \bar{\xi}}{\partial X} \frac{\partial^2 \bar{\xi}}{\partial X^2}
    \right)
    \right),
\end{align}
$$

$$
\begin{align}
    \frac{\bar{\xi} \, \partial \bar{T}}{\partial \bar{t}}
    &= \frac{1}{\text{Pr}} \frac{\partial}{\partial X} \left( \bar{\xi} \frac{\partial \bar{T}}{\partial X} \right)
    - \frac{1}{2 \, \text{Pr}} \left( \frac{\partial \bar{\xi}}{\partial X} \right)^2
    - \frac{\beta}{\text{Pr}} f \notag \\
    &\quad + \left(
    \frac{\bar{\xi}^3}{3} \frac{\partial \bar{\Pi}}{\partial X}
    + \frac{\text{Ma}}{2} \bar{\xi}^2 \frac{\partial f}{\partial X}
    \right) \frac{\partial f}{\partial X} \notag \\
    &\quad + \frac{\partial}{\partial X} \left(
    \frac{\bar{\xi}^4}{8} \frac{\partial \bar{\Pi}}{\partial X}
    + \frac{\text{Ma}}{6} \bar{\xi}^3 \frac{\partial f}{\partial X}
    \right)
\end{align}
$$

where $\bar{\Pi} = G \bar{\xi} - C \frac{\partial^2 \bar{\xi}}{\partial X^2} + \frac{(b \Omega)^2}{2} \Re\{H\} - \varphi $, $f=\bar{T} - \bar{\xi}$ is the temperature perturbation at the film surface, and

$$
\begin{align*}
	H = - \frac{\partial}{\partial X} \left[ \bar{\xi} F(\theta \bar{\xi}) \frac{\partial \bar{\xi}}{\partial X} \right], F (\theta \bar{\xi}) = 1 - \frac{\tan(\theta \bar{\xi})}{\theta \bar{\xi}},
\end{align*}
$$
$$
\begin{align*}
	Q_1 = 3 \frac{4 \sinh \gamma \sin \gamma - \gamma (\psi_+ \varphi_+ - \psi_- \varphi_-)}{2 \gamma^2 \psi_+^2},
\end{align*}
$$
$$
\begin{align*}
	Q_2 = -\frac{1}{3} + \frac{11 \varphi_- - 3 \gamma \psi_-}{\gamma^3 \psi_+}, \gamma=\bar{\xi} \sqrt{2 \omega},
\end{align*}
$$
$$
\begin{align*}
	\psi_+ = \cosh \gamma + \cos \gamma,
	\psi_- = \cosh \gamma - \cos \gamma, \\
	\varphi_+ = \sinh \gamma + \sin \gamma,
	\varphi_- = \sinh \gamma - \sin \gamma.
\end{align*}
$$

The `py-pde` python package provides methods and classes also useful for solving these type of DE.  