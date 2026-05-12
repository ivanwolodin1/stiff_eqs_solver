import os

directory = '../results/Ma_0334'   # ← путь к каталогу

data = []

for name in os.listdir(directory):
    full_path = os.path.join(directory, name)
    if not os.path.isdir(full_path):
        continue

    parts = name.split('_')
    parsed = {}

    for part in parts:
        # каждая часть теперь должна быть вида "key=value"
        if '=' not in part:
            continue
        key, value = part.split('=')
        parsed[key] = value
    # print(parsed)
    # exit(1)
    # проверим, что нужные ключи есть
    if not {'ampl', 'omfr', 'solved'}.issubset(parsed.keys()):
        print('Пропускаю, не хватает ключей:', name)
        continue

    ampl = float(parsed['ampl'])
    omfr = float(parsed['omfr'])
    solved = parsed['solved'] == 'True'

    data.append((omfr, ampl, solved))

print(data[:10])
print('Всего элементов:', len(data))

output_file = os.path.join(directory, 'summary_table.txt')

data_sorted = sorted(
    data, key=lambda x: (x[1], x[0])
)  # сортируем по ampl, затем по omfr

with open(output_file, 'w') as f:
    f.write('omfr\tampl\tbool\n')
    for omfr, ampl, solved in data_sorted:
        f.write(f'{omfr}\t{ampl}\t{solved}\n')

print('Файл создан:', output_file)
