from os import makedirs

from numpy import arange

from params import get_simulation_params
from initial_conditions.init_fields import prepare_initial_conditions
from runner.run import run_simulation
from io_utils.save_txt import save_res_to_txt
from io_utils.plotter import plot_results_to_file


def params_to_folder_name(params: dict) -> str:
    excluded = {'t_range', 'dt', 'bc'}
    parts = []

    for key, value in params.items():
        if key in excluded:
            continue
        parts.append(f'{key}={value}')

    return '_'.join(parts)


def main():
    params = get_simulation_params()
    amplitudes = [
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.1,
        1.2,
        1.3,
        1.4,
        1.5,
        1.6,
        1.7,
        1.8,
        1.9,
        2.0,
    ]
    omega_frs = [
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.1,
        1.2,
        1.3,
        1.4,
        1.5,
        1.6,
        1.7,
        1.8,
        1.9,
        2.0,
    ]

    initial_h, initial_T = prepare_initial_conditions(
        perturbation_amplitude=0.3,
        wave_k=0.9,
    )

    for om in omega_frs:
        for ampl in amplitudes:
            params['b'] = ampl
            params['omega'] = om
            print(f'parameters: {params}')

            result, storage, is_solved = run_simulation(
                initial_h=initial_h, initial_T=initial_T, **params
            )
            print('Simulation stopped')

            # подпапка с параметрами начальных условий
            root_output_dir = (
                f'results/ampl={ampl}_omfr={om}_is_solved={is_solved}'
            )
            makedirs(root_output_dir, exist_ok=True)

            save_res_to_txt(storage, root_output_dir, is_solved)
            plot_results_to_file(storage, root_output_dir, is_solved)

    # print(f'parameters: {params}')

    # # -------------------------------------------------------------
    # # Корневая папка: имя формируется ИЗ ПАРАМЕТРОВ
    # # -------------------------------------------------------------
    # params_folder = params_to_folder_name(params)
    # root_output_dir = f'results/{params_folder}'
    # makedirs(root_output_dir, exist_ok=True)
    # print(f'Root folder: {root_output_dir}')

    # initial_amplitude_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    # wave_k_values = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

    # for initial_amplitude in initial_amplitude_values:
    #     for wave_k in wave_k_values:

    #         initial_h, initial_T = prepare_initial_conditions(
    #             perturbation_amplitude=initial_amplitude,
    #             wave_k=wave_k,
    #         )

    #         result, storage, is_solved = run_simulation(
    #             initial_h=initial_h, initial_T=initial_T, **params
    #         )
    #         print('Simulation stopped')

    #         # подпапка с параметрами начальных условий
    #         output_dir = (
    #             f'{root_output_dir}/'
    #             f'k={wave_k:.2f}_a={initial_amplitude}_is_solved={is_solved}'
    #         )
    #         makedirs(output_dir, exist_ok=True)

    #         save_res_to_txt(storage, output_dir, is_solved)
    #         plot_results_to_file(storage, output_dir, is_solved)


if __name__ == '__main__':
    main()
