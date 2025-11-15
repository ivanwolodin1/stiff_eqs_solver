from os import makedirs

import numpy as np

import matplotlib.pyplot as plt

from pde import (
    FieldCollection,
    PDEBase,
    CartesianGrid,
    ScalarField,
    MemoryStorage,
    PlotTracker,
)


class Stiff1DSystem(PDEBase):
    """Second equation from:
    Hairer, E. and Wanner, G. (1996) Solving Ordinary Differential Equations II,
    Stiff Differential-Algebraic Problems. 2nd Edition, Springer-Verlag, Berlin.
    """

    def __init__(
        self,
        alpha=0,
        A=0,
        B=0,
        bc={
            'x-': {'value': 1},
            'x+': {'value': 1},
        },
    ):
        super().__init__()
        self.bc = bc
        self.alpha = alpha
        self.A = A
        self.B = B

    # why do we have to pass t to this overloaded function, even if we don't use it?
    def evolution_rate(self, state, t=0):
        u, v = state

        u_t = (
            self.A
            + u**2 * v
            - (self.B + 1.0) * u
            + self.alpha * u.laplace(bc=self.bc)
        )
        v_t = self.B * u - u**2 * v + self.alpha * v.laplace(bc=self.bc)

        return FieldCollection([u_t, v_t])


Nx = 100
grid = CartesianGrid([[0, 1.0]], Nx)
x = grid.axes_coords[0]

# initial conditions
# https://py-pde.readthedocs.io/en/latest/examples_gallery/simple_pdes/pde_brusselator_expression.html#sphx-glr-examples-gallery-simple-pdes-pde-brusselator-expression-py
u0 = 1.0 + np.sin(2 * np.pi * x)
v0 = np.ones_like(x) * 3.0
u_initial = ScalarField(grid, u0)
v_initial = ScalarField(grid, v0)

state = FieldCollection([u_initial, v_initial])

storage = MemoryStorage()
trackers = [
    'progress',  # show progress bar during simulation
    'steady_state',  # abort when steady state is reached
    storage.tracker(interrupts=1),  # store data every simulation time unit
    PlotTracker(show=True),  # show images during simulation
]

eq = Stiff1DSystem(
    alpha=1.0 / 50.0,
    A=1.0,
    B=3.0,
)
result = eq.solve(
    state,
    t_range=10000,
    dt=0.001,
    tracker=trackers,
)


# --- Save steps ---
makedirs('stiff1D_system', exist_ok=True)

for i, step in enumerate(storage.data):
    y_data = step.data
    np.savetxt(f'stiff1D_system/y_{i}.txt', y_data)

result.plot()
