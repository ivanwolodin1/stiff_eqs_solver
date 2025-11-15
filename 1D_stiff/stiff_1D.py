import time
from os import makedirs

import matplotlib.pyplot as plt
from numpy import savetxt

from pde import CartesianGrid, ScalarField, PDE, MemoryStorage

storage = MemoryStorage()

supported_solvers = {
    'scipy': 'scipy',
    'adams–bashforth': 'adams–bashforth',
    'crank-nicolson': 'crank-nicolson',
    'explicit': 'explicit',
    'implicit': 'implicit',
}

grid = CartesianGrid([[0, 1.5]], 38)
state = ScalarField(grid, data=0.0)

# First example:
# Hairer, E. and Wanner, G. (1996) Solving Ordinary Differential Equations II,
# Stiff Differential-Algebraic Problems. 2nd Edition, Springer-Verlag, Berlin.
eq_to_solve = '-50 * y + 50 * cos(x)'
eq = PDE({'y': eq_to_solve})

t0 = time.time()
result_scipy = eq.solve(
    state,
    t_range=[0, 20],
    dt=0.03948,
    solver=supported_solvers['explicit'],
    tracker=['progress', storage.tracker(1)],
)
t_scipy = time.time() - t0

print(f'SciPy solver took {t_scipy:.3f} seconds.')

dir_to_save = 'stiff_1D'
makedirs(dir_to_save, exist_ok=True)

# --- Save each step ---
for i, step in enumerate(storage.data):
    y_data = step.data

    savetxt(f'{dir_to_save}/y_{i}.txt', y_data)
