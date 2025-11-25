from pde import MemoryStorage, FieldCollection, PlotTracker
from pde.solvers.scipy import ScipySolverError

from pde_model.volodin_2025 import Volodin2025
from pde_model.shklyaev_2008 import Shklyaev2008, Shklyaev2008_VIII_B


def run_simulation(
    pde_model,
    t_range=None,
    dt=None,
    Ca=None,
    G0=None,
    b=None,
    omega=None,
    phi=None,
    Ma=None,
    Pr=None,
    Bi=None,
    V=None,
    k=None,
    bc=None,
    initial_h=None,
    initial_T=None,
):
    print(f'Running simulation for model: {pde_model}')
    print(f'Time range: {t_range}, dt: {dt}')
    print(
        f'Parameters: Ca={Ca}, G0={G0}, b={b}, omega={omega}, phi={phi}, Ma={Ma}, Bi={Bi}, Pr={Pr}, V={V}, k={k}'
    )
    print(f'Boundary conditions: {bc}')

    storage = MemoryStorage()
    trackers = [
        'progress',
        PlotTracker(show=True),
        storage.tracker(interrupts=1),
    ]

    print(f'PDE model={pde_model}')
    if pde_model == 'volodin_2025':
        eq = Volodin2025(
            Ca=Ca, G0=G0, b=b, omega=omega, phi=phi, Ma=Ma, Bi=Bi, Pr=Pr, bc=bc
        )
        state = FieldCollection([initial_h, initial_T])
    elif pde_model == 'shklyaev_2008':
        eq = Shklyaev2008(Ca=Ca, G0=G0, b=b, omega=omega, phi=phi, bc=bc)
        state = initial_h
    elif pde_model == 'shklyaev_2008_viii_b':
        eq = Shklyaev2008_VIII_B(V=V, G0=G0, k=k, omega=omega, bc=bc)
        state = initial_h
    else:
        raise ValueError(f'Unknown pde_model: {pde_model}')

    print('Starting solver...')
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
