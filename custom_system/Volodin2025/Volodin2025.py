# https://doi.org/10.1007/s12217-025-10201-8
# Solution to 57-57 eqs.

from os import makedirs

from matplotlib import pyplot as plt
from numpy import sqrt, sin, sinh, cos, cosh, pi, tan, real, savetxt

# from numba import njit

from pde import (
    PDEBase,
    CartesianGrid,
    ScalarField,
    MemoryStorage,
    PlotTracker,
    FieldCollection,
)


class Volodin2025(PDEBase):
    """https://doi.org/10.1007/s12217-025-10201-8
    Solution to 57-57 eqs.
    """

    def __init__(
        self,
        Ca=0,
        G0=0,
        b=0,
        omega=10,
        phi=0,
        Ma=0,
        Bi=0,
        Pr=0,
        bc={'x': 'periodic'},
    ):
        super().__init__()
        self.bc = bc
        self.Ca = Ca
        self.G0 = G0
        self.b = b
        self.omega = omega
        self.phi = phi
        self.Ma = Ma
        self.Bi = Bi
        self.Pr = Pr
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

    def evolution_rate(self, state, t=0):
        h, T = state
        f = T - h
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
            + 0.5 * self.Ma * h**2 * f.gradient(bc=self.bc)[0]
            - 0.5 * self.b**2 * self.omega**2 * Q
        ).gradient(bc=self.bc)[0]

        dx_T = T.gradient(bc=self.bc)[0]
        dx_f = f.gradient(bc=self.bc)[0]

        # (1/Pr) * d/dx ( h * dT/dx )
        term1 = (h * dx_T / self.Pr).gradient(bc=self.bc)[0]

        # - (1/(2Pr)) * (dh/dx)^2
        term2 = -(dx_h**2) / (2 * self.Pr)

        # - beta/Pr * f
        term3 = -(self.Bi / self.Pr) * f

        # ( h^3/3 * dPi/dx + Ma/2 * h^2 * df/dx ) * df/dx
        term4 = ((h**3 / 3) * dx_Pi + 0.5 * self.Ma * h**2 * dx_f) * dx_f

        # term5: d/dx( h^4/8 * dPi/dx + Ma/6 * h^3 * df/dx )
        inside = (h**4 / 8) * dx_Pi + (self.Ma / 6) * h**3 * dx_f
        term5 = inside.gradient(bc=self.bc)[0]

        # h * T_t = term1 + term2 + term3 + term4 + term5
        T_t = (term1 + term2 + term3 + term4 + term5) / h

        return FieldCollection([h_t, T_t])


def prepare_simulation_params():
    # simulation parameters

    # ======================================================================
    # Grid and time
    # ======================================================================

    Nx = 100
    L = 5
    t_range = 100
    dt = 0.0001

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

    Ca = surface_tension * (film_thickness / (density * viscosity**2))
    G0 = gravitational_acceleration * (film_thickness**3 / viscosity**2)
    omega = frequency * (film_thickness**2 / viscosity)
    b = amplitude * 10.0
    phi = Hamaker_constant / (
        6.0 * 3.14 * density * viscosity**2 * film_thickness
    )

    # omega *= 10000
    # phi *= 10000

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

    # # stable for dt = 1e-6
    # G0 = 1e-8
    # Ca = 7.0
    # omega = 0.2
    # phi = 1e-5
    # V = 1.0
    # b = sqrt((2 * V * Ca) / (omega**2))

    # Boundary conditions
    bc = {'x': 'periodic'}

    # Initial conditions
    h0 = 1.0 + 0.1 * sin(2 * pi / L * x)
    initial_h = ScalarField(grid, h0)

    T0 = 1.0 + 0.1 * sin(2 * pi / L * x)
    initial_T = ScalarField(grid, T0)

    return {
        't_range': t_range,
        'dt': dt,
        'Ca': Ca,
        'G0': G0,
        'b': b,
        'omega': omega,
        'phi': phi,
        'bc': bc,
        'initial_h': initial_h,
        'Ma': 1e-5,
        'Pr': 7.0,
        'Bi': 1e-6,
        'initial_T': initial_T,
    }


def run_simulation(
    t_range=0,
    dt=0.0,
    Ca=0,
    G0=0,
    b=0,
    omega=0,
    phi=0,
    Ma=0,
    Pr=0,
    Bi=0,
    bc={'x': 'periodic'},
    initial_h=None,
    initial_T=None,
):

    # Storage and trackers
    storage = MemoryStorage()
    trackers = [
        'progress',
        # 'steady_state',
        storage.tracker(interrupts=1),
        PlotTracker(show=True),
        # PlotTracker(
        #     show=True,
        #     title='Evolution of h(x) over time',
        #     transformation=lambda state, t: state[0],  # только h
        # ),
        # PlotTracker(
        #     show=True,
        #     title='Evolution of T(x) over time',
        #     transformation=lambda state, t: state[1],  # только T
        # ),
    ]

    eq = Volodin2025(
        Ca=Ca,
        G0=G0,
        b=b,
        omega=omega,
        phi=phi,
        Ma=Ma,
        Pr=Pr,
        Bi=Bi,
        bc=bc,
    )

    state = FieldCollection([initial_h, initial_T])

    result = eq.solve(
        state,
        t_range=t_range,
        dt=dt,
        tracker=trackers,
        solver='scipy',
        method='Radau',
        rtol=1e-6,  # relative tolerance
        atol=1e-9,  # absolute tolerance
    )

    return result, storage


def save_res_to_txt(storage, output_dir='res'):
    makedirs(output_dir, exist_ok=True)

    for i, step in enumerate(storage.data):
        h_field = step[0].data
        T_field = step[1].data
        savetxt(f'{output_dir}/h_{i}.txt', h_field)
        savetxt(f'{output_dir}/T_{i}.txt', T_field)

    print(f'Сохранено {len(storage.data)} шагов в {output_dir}/')


def plot_results(storage):
    x = storage[0].grid.cell_coords[:, 0]

    fig, axes = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

    # --- h ---
    axes[0].plot(
        x, storage[0][0].data, label=f'Initial step t={storage.times[0]:.2f}'
    )
    axes[0].plot(
        x, storage[-1][0].data, label=f'Final step t={storage.times[-1]:.2f}'
    )
    axes[0].set_ylabel('h(x)')
    axes[0].set_title('Evolution of h(x)')
    axes[0].legend()
    axes[0].grid(True)

    # --- T ---
    axes[1].plot(
        x, storage[0][1].data, label=f'Initial step t={storage.times[0]:.2f}'
    )
    axes[1].plot(
        x, storage[-1][1].data, label=f'Final step t={storage.times[-1]:.2f}'
    )
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('T(x)')
    axes[1].set_title('Evolution of T(x)')
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout()
    plt.show()


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

    Ma = params['Ma']
    Pr = params['Pr']
    Bi = params['Bi']

    initial_h = params['initial_h']
    initial_T = params['initial_T']

    print(
        f'Ca={Ca}, G0={G0}, b={b}, omega={omega}, phi={phi}, Ma={Ma}, Pr={Pr}, Bi={Bi}'
    )

    result, storage = run_simulation(
        t_range=t_range,
        dt=dt,
        Ca=Ca,
        G0=G0,
        b=b,
        omega=omega,
        phi=phi,
        Ma=Ma,
        Pr=Pr,
        Bi=Bi,
        bc=bc,
        initial_h=initial_h,
        initial_T=initial_T,
    )
    print('Simulation stopped')
    
    save_res_to_txt(storage, output_dir='res')
    plot_results(storage)


if __name__ == '__main__':
    main()
