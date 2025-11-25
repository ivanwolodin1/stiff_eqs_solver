import os

directory = "results/amplitude_frequency_map/apl_fre"   # ← сюда путь

data = []

for name in os.listdir(directory):
    if not os.path.isdir(os.path.join(directory, name)):
        continue

    parts = name.split("_")
    parsed = {}

    for part in parts:
        if part == 'is':
            continue
        key, value = part.split("=")
        parsed[key] = value

    print(parsed)
    ampl = float(parsed["ampl"])
    omfr = float(parsed["omfr"])
    solved = parsed["solved"] == "True"

    data.append((omfr, ampl, solved))

print(data[:10])
print("Всего элементов:", len(data))

output_file = "results/amplitude_frequency_map/apl_fre/summary_table.txt"

data_sorted = sorted(data, key=lambda x: (x[1], x[0]))

with open(output_file, "w") as f:
    f.write("omfr\tampl\tbool\n")
    for omfr, ampl, solved in data_sorted:
        f.write(f"{omfr}\t{ampl}\t{solved}\n")

print("Файл создан:", output_file)
