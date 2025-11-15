from os import makedirs

import numpy as np

from pde import ScalarField, UnitGrid, PlotTracker, PDE, MemoryStorage

grid = UnitGrid([500, 100], periodic=[True, False])
state = ScalarField(grid, data=0.0)

bcs = {'y-': {'value': 0}, 'y+': {'value': 1}, 'x': 'periodic'}

eq = PDE(
    {'u': 'laplace(u)'},
    bc=bcs,
)

storage = MemoryStorage()
trackers = [
    'progress',  # show progress bar during simulation
    'steady_state',  # abort when steady state is reached
    storage.tracker(interrupts=100),  # store data every simulation time unit
    PlotTracker(show=True),  # show images during simulation
]

result = eq.solve(
    state,
    t_range=[0, 100000],
    # dt=1e-3,
    # solver='scipy',  # surprisingly it is better not to specify the solver! default works best
    tracker=trackers,
)

# --- Save steps ---
makedirs('laplace_periodic', exist_ok=True)

for i, step in enumerate(storage.data):
    y_data = step.data
    np.savetxt(f'laplace_periodic/y_{i}.txt', y_data)

result.plot()
