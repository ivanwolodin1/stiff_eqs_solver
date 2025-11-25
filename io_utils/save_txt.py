from os import makedirs
from numpy import savetxt


def save_res_to_txt(storage, output_dir=None, is_solved=None):
    makedirs(output_dir, exist_ok=True)

    # создаем поддиректории
    h_dir = f'{output_dir}/h'
    T_dir = f'{output_dir}/T'
    makedirs(h_dir, exist_ok=True)
    makedirs(T_dir, exist_ok=True)

    # сохраняем поля
    for i, step in enumerate(storage.data):
        h_field = step[0].data
        T_field = step[1].data

        savetxt(f'{h_dir}/h_{i}.txt', h_field)
        savetxt(f'{T_dir}/T_{i}.txt', T_field)

    print(f'Сохранено {len(storage.data)} шагов в {output_dir}/')

    # создаем текстовый файл статуса
    status_path = f'{output_dir}/status.txt'
    with open(status_path, 'w') as f:
        if is_solved is True:
            f.write('solved')
        elif is_solved is False:
            f.write('not_solved')
        else:
            f.write('unknown')

    print(f'Статус записан в {status_path}')
