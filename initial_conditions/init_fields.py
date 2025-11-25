from numpy import pi, cos
from pde import CartesianGrid, ScalarField


def prepare_initial_conditions(perturbation_amplitude, wave_k):
    Nx = 100
    wavelength = 2 * pi / wave_k
    grid = CartesianGrid([[0, wavelength]], [Nx], periodic=True)

    x = grid.axes_coords[0]
    h0 = 1 + perturbation_amplitude * cos(wave_k * x)
    T0 = 1 + perturbation_amplitude * cos(wave_k * x)

    return ScalarField(grid, h0), ScalarField(grid, T0)
