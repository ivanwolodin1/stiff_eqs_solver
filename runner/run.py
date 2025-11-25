from pde import MemoryStorage, FieldCollection
from pde.solvers.scipy import ScipySolverError

from pde_model.volodin_2025 import Volodin2025


def run_simulation(
    t_range, dt, Ca, G0, b, omega, phi, Ma, Pr, Bi, bc, initial_h, initial_T
):

    storage = MemoryStorage()
    trackers = [
        'progress',
        storage.tracker(interrupts=1),
    ]

    eq = Volodin2025(Ca, G0, b, omega, phi, Ma, Bi, Pr, bc)
    state = FieldCollection([initial_h, initial_T])

    try:
        result = eq.solve(
            state,
            t_range=t_range,
            solver='scipy',
            method='Radau',
            tracker=trackers,
            rtol=1e-6,
            atol=1e-9,
        )
    except ScipySolverError:
        print('Solver error: stopping simulation and returning last state')
        if len(storage.data) >= 2:
            return storage.data[-2], storage, False
        return state, storage, False

    return result, storage, True
