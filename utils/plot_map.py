import matplotlib.pyplot as plt

# --- читаем файл ---
data = []
path_to_file = (
    '../results/Ma_0167/amplitude_freq_extended/res/summary_table.txt'
)
# path_to_file = '../results/Ma_00835/summary_table.txt'
# path_to_file = '../results/Ma_0334/summary_table.txt'

with open(path_to_file) as f:
    next(f)  # пропускаем заголовок
    for line in f:
        parts = line.strip().split()
        if len(parts) != 3:
            continue
        omfr = float(parts[0])
        ampl = float(parts[1])
        solved = parts[2] == 'True'
        data.append((omfr, ampl, solved))

# --- создаём списки для графика ---
xs = [d[0] for d in data]
ys = [d[1] for d in data]
colors = ['green' if d[2] else 'red' for d in data]

# --- рисуем ---
plt.figure(figsize=(10, 8))
plt.scatter(xs, ys, c=colors, s=30)

plt.xlabel('Частота')
plt.ylabel('Амплитуда')
plt.title('Карта устойчивости')
plt.scatter([], [], c='green', label='Устойчиво')
plt.scatter([], [], c='red', label='Неустойчиво')
plt.legend(loc='best')
plt.grid(True)

plt.show()
