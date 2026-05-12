import matplotlib.pyplot as plt
from collections import defaultdict

files = [
    (
        'Ma = 0.0167',
        '../results/Ma_0167/amplitude_freq_extended/res/summary_table.txt',
    ),
    ('Ma = 0.00835', '../results/Ma_00835/summary_table.txt'),
    ('Ma = 0.0334', '../results/Ma_0334/summary_table.txt'),
]


def load_data(path):
    data = []
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            omfr = float(parts[0])
            ampl = float(parts[1])
            solved = parts[2] == 'True'
            data.append((omfr, ampl, solved))
    return data


def extract_neutral_curve(data):
    # группируем по частоте
    grouped = defaultdict(list)
    for ω, A, stable in data:
        grouped[ω].append((A, stable))

    neutral_ω = []
    neutral_A = []

    for ω, items in grouped.items():
        items.sort()  # сортировка по амплитуде

        # ищем первую смену стабильности
        for k in range(1, len(items)):
            prev_A, prev_stable = items[k - 1]
            A, stable = items[k]

            if prev_stable != stable:
                # точка перехода между двумя экспериментами
                mid_A = 0.5 * (prev_A + A)   # условная нейтральная точка
                neutral_ω.append(ω)
                neutral_A.append(mid_A)
                break

    # сортируем кривую по частоте
    sorted_pairs = sorted(zip(neutral_ω, neutral_A))
    return zip(*sorted_pairs) if sorted_pairs else ([], [])


plt.figure(figsize=(11, 8))

for label, path in files:
    data = load_data(path)
    ω, A = extract_neutral_curve(data)
    plt.plot(ω, A, label=f'Нейтральная кривая {label}', linewidth=2)

plt.xlabel('Частота')
plt.ylabel('Амплитуда')
plt.title('Нейтральные кривые для трёх карт')
plt.grid(True)
plt.legend()
plt.show()
