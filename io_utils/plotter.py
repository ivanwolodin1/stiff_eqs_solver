import matplotlib.pyplot as plt
from os import makedirs


def plot_results_to_file(storage, output_dir=None, is_solved=True):
    makedirs(output_dir, exist_ok=True)

    # выбор индекса шага
    if is_solved is True:
        idx = -1                             # финальный шаг
    else:
        idx = (
            -2 if len(storage) > 1 else -1
        )   # предпоследний, если он существует

    x = storage[0].grid.cell_coords[:, 0]

    fig, axes = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

    # --- h ---
    axes[0].plot(
        x, storage[0][0].data, label=f'Initial t={storage.times[0]:.2f}'
    )
    axes[0].plot(
        x, storage[idx][0].data, label=f'Step {idx} t={storage.times[idx]:.2f}'
    )
    axes[0].set_ylabel('h(x)')
    axes[0].set_title('Evolution of h(x)')
    axes[0].legend()
    axes[0].grid(True)

    # --- T ---
    axes[1].plot(
        x, storage[0][1].data, label=f'Initial t={storage.times[0]:.2f}'
    )
    axes[1].plot(
        x, storage[idx][1].data, label=f'Step {idx} t={storage.times[idx]:.2f}'
    )
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('T(x)')
    axes[1].set_title('Evolution of T(x)')
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout()

    # путь для сохранения рисунка
    out_path = f'{output_dir}/plot.png'
    plt.savefig(out_path, dpi=200)

    # === 2) сохранение в папку plots/ ===
    plots_dir = 'plots'
    makedirs(plots_dir, exist_ok=True)

    # имя файла = имя папки output_dir + ".png"
    # например output_dir="a=0.1_k=0.50_is_solved=False" →
    # plots/a=0.1_k=0.50_is_solved=False.png
    output_name = output_dir.strip('/').split('/')[
        -1
    ]  # только последний сегмент
    plots_path = f'{plots_dir}/{output_name}.png'
    plt.savefig(plots_path, dpi=200)

    plt.close()

    print(f'График сохранён в:\n - {out_path}\n - {plots_path}')
