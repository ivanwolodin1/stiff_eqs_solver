# https://doi.org/10.1103/PhysRevE.77.036320
# Solution to 30(a)-31(b) eq.

from os import makedirs

import numpy as np
from numpy import sqrt, sin, sinh, cos, cosh, tan, real

# from numba import njit

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
        # self.alpha = (1j - 1) * (sqrt(2.0 * self.omega)) / 2.0
        self.alpha = -(1j - 1) * (sqrt(2.0 * self.omega)) / 2.0

        self.real_alpha = 1.0 - (sqrt(2) / (2 * sqrt(self.omega))) * (
            (
                sin(sqrt(2) * sqrt(self.omega) / 2)
                * cos(sqrt(2) * sqrt(self.omega) / 2)
                + sinh(sqrt(2) * sqrt(self.omega) / 2)
                * cosh(sqrt(2) * sqrt(self.omega) / 2)
            )
            / (
                cos(sqrt(2) * sqrt(self.omega) / 2) ** 2
                + cosh(sqrt(2) * sqrt(self.omega) / 2) ** 2
                - 1
            )
        )

    def evolution_rate(self, h, t=0):
        # https://py-pde.readthedocs.io/en/latest/examples_gallery/advanced_pdes/pde_1d_class.html
        # on gradient() operator

        dx_h = h.gradient(bc=self.bc)[0]

        # second derivatives
        d2x_h = h.laplace(bc=self.bc)

        # --- Π ---
        # f_alpha = 1.0 - (tan(self.alpha * h)) / (self.alpha * h)
        # real_f_alpha = real(f_alpha)
        real_f_alpha = self.real_alpha

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


def prepare_simulation_params():
    # simulation parameters

    # ======================================================================
    # Grid and time
    # ======================================================================

    Nx = 100
    L = 5
    t_range = 100
    dt = 0.000001

    grid = CartesianGrid([[0, L]], [Nx], periodic=True)
    x = grid.axes_coords[0]

    # ======================================================================
    # Physical parameters
    # ======================================================================

    film_thickness = 1e-7  # m
    density = 1000  # kg/m^3
    viscosity = 10**-6  # m^2/s
    gravitational_acceleration = 9.81  # m/s^2
    frequency = 1000  # Hz
    amplitude = film_thickness   # m
    surface_tension = 0.07  # N/m
    Hamaker_constant = 6.0 * 3.14 * 1e-21  # J

    # Ca = surface_tension * (film_thickness / (density * viscosity**2))
    # G0 = gravitational_acceleration * (film_thickness**3 / viscosity**2)
    # omega = frequency * (film_thickness**2 / viscosity)
    # b = amplitude * 10.0
    # phi = Hamaker_constant / (
    #     6.0 * 3.14 * density * viscosity**2 * film_thickness
    # )

    # typical values for:
    # film_thickness = 1e-7  # m
    # density = 1000  # kg/m^3
    # viscosity = 10**-6  # m^2/s
    # gravitational_acceleration = 9.81  # m/s^2
    # frequency = 1000  # Hz
    # amplitude = film_thickness * 10   # m
    # surface_tension = 0.07  # N/m
    # Hamaker_constant = 6.0 * 3.14 * 1e-21  # J

    # a=7.0, G0=1e-08, b=1e-06, omega=1e-05, phi=1e-05

    # stable:
    # Ca = 1e-4
    # G0 = 1.1e-4
    # b = 0.05
    # omega = 10
    # phi = 1e-4

    # ======================================================================

    G0 = 1e-8
    Ca = 7.0
    omega = 0.2
    phi = 1e-5

    V = 1.0
    b = sqrt((2 * V * Ca) / (omega**2))

    # Boundary conditions
    bc = {'x': 'periodic'}

    # Initial conditions
    h0 = 1.0 + 1e-4 * np.sin(2 * np.pi / L * x)
    h_initial = ScalarField(grid, h0)

    return {
        't_range': t_range,
        'dt': dt,
        'Ca': Ca,
        'G0': G0,
        'b': b,
        'omega': omega,
        'phi': phi,
        'bc': bc,
        'h_initial': h_initial,
    }


def run_simulation(
    t_range=0,
    dt=0.0,
    Ca=0,
    G0=0,
    b=0,
    omega=0,
    phi=0,
    bc={'x': 'periodic'},
    initial_conditions=None,
):

    # Storage and trackers
    storage = MemoryStorage()
    trackers = [
        'progress',
        # 'steady_state',
        storage.tracker(interrupts=1),
        PlotTracker(show=True),
    ]

    eq = Shklyaev2008(
        Ca=Ca,
        G0=G0,
        b=b,
        omega=omega,
        phi=phi,
        bc=bc,
    )

    result = eq.solve(
        initial_conditions,
        t_range=t_range,
        dt=dt,
        tracker=trackers,
    )

    return result, storage


def save_res_to_txt(storage, output_dir='Shklyaev2008'):
    makedirs(output_dir, exist_ok=True)

    for i, step in enumerate(storage.data):
        np.savetxt(f'{output_dir}/h_{i}.txt', step.data)


def main():
    params = prepare_simulation_params()
    t_range = params['t_range']
    dt = params['dt']
    Ca = params['Ca']
    G0 = params['G0']
    b = params['b']
    omega = params['omega']
    phi = params['phi']
    bc = params['bc']
    h_initial = params['h_initial']

    print(f'Ca={Ca}, G0={G0}, b={b}, omega={omega}, phi={phi}')

    result, storage = run_simulation(
        t_range=t_range,
        dt=dt,
        Ca=Ca,
        G0=G0,
        b=b,
        omega=omega,
        phi=phi,
        bc=bc,
        initial_conditions=h_initial,
    )

    save_res_to_txt(storage)

    result.plot()
    plt.show()


if __name__ == '__main__':
    main()
