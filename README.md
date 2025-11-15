# Solving Stiff Differential Equations with `py-pde`

This repository shows how to solve **stiff differential equations** using the  
[`py-pde`](https://py-pde.readthedocs.io/en/stable/) package (v.0.47.2).  
Two examples (so far) are included:

1. A **scalar stiff ODE** from *Hairer & Wanner (1996)*  
2. **Simple heat equation**
2. A **Custom 1D stiff PDE system** 

---

# 1. Stiff ODE Example

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

But the altimate goal is to solve the stiff system of eqs that we obtained and published in MS&T:

**Volodin, I., Alabuzhev, A.** *Linear Stability of Marangoni Convection in a Thin Film under Vertical Vibrations.* Microgravity Sci. Technol. 37, 48 (2025). [https://doi.org/10.1007/s12217-025-10201-8](https://doi.org/10.1007/s12217-025-10201-8)