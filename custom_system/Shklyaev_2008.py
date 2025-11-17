# https://doi.org/10.1103/PhysRevE.77.036320
# Solution to 30(a)-31(b) eq.

from os import makedirs

import numpy as np
from numpy import sqrt, sin, sinh, cos, cosh, tan, real

import matplotlib.pyplot as plt

from pde import (
    PDEBase,
    CartesianGrid,
    ScalarField,
    MemoryStorage,
    PlotTracker,
)


class Shklyaev2008(PDEBase):
    """https://doi.org/10.1103/PhysRevE.77.036320
    Solution to 30(a)-31(b) eq.
    """

    def __init__(
        self,
        Ca=0,
        G0=0,
        b=0,
        omega=10,
        phi=0,
        bc={'x': 'periodic'},
    ):
        super().__init__()
        self.bc = bc
        self.Ca = Ca
        self.G0 = G0
        self.b = b
        self.omega = omega
        self.phi = phi
        self.alpha = (1j + 1) * (sqrt(2.0 * self.omega)) / 2.0

    def evolution_rate(self, h, t=0):
        # https://py-pde.readthedocs.io/en/latest/examples_gallery/advanced_pdes/pde_1d_class.html
        # on gradient() operator

        dx_h = h.gradient(bc=self.bc)[0]

        # second derivatives
        d2x_h = h.laplace(bc=self.bc)

        # --- Π ---
        f_alpha = 1 - (tan(self.alpha * h)) / (self.alpha * h)
        real_f_alpha = real(f_alpha)

        H = -(h * real_f_alpha * dx_h).gradient(bc=self.bc)[0]

        Pi = (
            self.phi / h**3
            - self.Ca * d2x_h
            + self.G0 * h
            + 0.5 * self.b**2 * self.omega**2 * H
        )

        # derivative Pi_x
        dx_Pi = Pi.gradient(bc=self.bc)[0]

        # --- Q ---
        gamma = sqrt(2 * self.omega) * h
        ksi_plus = cosh(gamma) + cos(gamma)
        ksi_minus = cosh(gamma) - cos(gamma)
        phi_plus = sinh(gamma) + sin(gamma)
        phi_minus = sinh(gamma) - sin(gamma)

        Q1 = 3.0 * (
            (
                4.0 * sinh(gamma) * sin(gamma)
                - gamma * (ksi_plus * phi_plus - ksi_minus * phi_minus)
            )
            / (2 * gamma * gamma * ksi_plus * ksi_plus)
        )

        Q2 = -1.0 / 3.0 + (11.0 * phi_minus - 3.0 * gamma * ksi_minus) / (
            gamma**3 * ksi_plus
        )

        Q = Q1 * h**2 * (dx_h) ** 3 + Q2 * h**3 * dx_h * d2x_h

        h_t = (
            (1.0 / 3.0) * h**3 * dx_Pi
            - 0.5 * self.b**2 * self.omega**2 * Q
        ).gradient(bc=self.bc)[0]

        return h_t


Nx = 100
grid = CartesianGrid([[0, 5]], [Nx], periodic=True)
x = grid.axes_coords[0]

# h0 = np.ones_like(x)
h0 = 1.0 + 0.0001 * np.sin(2 * np.pi / 5 * x)
h_initial = ScalarField(grid, h0)

storage = MemoryStorage()
trackers = [
    'progress',  # show progress bar during simulation
    # 'steady_state',  # abort when steady state is reached
    storage.tracker(interrupts=50),  # store data every simulation time unit
    PlotTracker(show=True),  # show images during simulation
]

eq = Shklyaev2008(
    Ca=0.0001,
    G0=0.00011,
    b=0.05,
    omega=10,
    phi=0.0001,
)
result = eq.solve(
    h_initial,
    t_range=100,
    dt=0.0001,
    tracker=trackers,
)

# --- Save steps ---
makedirs('Shklyaev2008', exist_ok=True)

for i, step in enumerate(storage.data):
    y_data = step.data
    np.savetxt(f'Shklyaev2008/h_{i}.txt', y_data)

result.plot()
