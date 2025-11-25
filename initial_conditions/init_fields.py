from numpy import pi, cos
from pde import CartesianGrid, ScalarField


def prepare_initial_conditions(perturbation_amplitude, wave_k):
    # Nx = 100
    # wavelength = 2 * pi / wave_k
    # grid = CartesianGrid([[0, wavelength]], [Nx], periodic=True)

    # x = grid.axes_coords[0]
    # h0 = 1 + perturbation_amplitude * cos(wave_k * x)
    # T0 = 1 + perturbation_amplitude * cos(wave_k * x)

    # return ScalarField(grid, h0), ScalarField(grid, T0)

    Nx = 100
    L = 1.0                      # длина области
    grid = CartesianGrid([[0, L]], [Nx], periodic=True)

    x_grid = grid.axes_coords[0]     # координаты в диапазоне [0,1]

    # Физическая длина волны определяется wave_k
    wavelength = 2 * pi / wave_k

    # Растягиваем координаты до физической длины волны
    x_phys = x_grid * wavelength

    h0 = 1 + perturbation_amplitude * cos(wave_k * x_phys)
    T0 = 1 + perturbation_amplitude * cos(wave_k * x_phys)

    return ScalarField(grid, h0), ScalarField(grid, T0)
