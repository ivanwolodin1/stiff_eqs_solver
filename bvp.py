import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# --- 1. Параметры задачи ---
L = 10.0          # Длина интервала
alpha = 2.1      # Параметр alpha
eps = 0.01       # Малый параметр epsilon

# --- 2. Функция ОДУ (fun) ---
# y[0] = h, y[1] = dh/dx, y[2] = Integral(h dx)
# p[0] = Const
def fun(x, y, p):
    y1 = y[0]  # Это массив h(x)
    y2 = y[1]
    Const = p[0]

    # ⚠️ ИСПРАВЛЕНИЕ: Используем np.where для векторизованной проверки
    # -----------------------------------------------------------

    # Вычисляем член f = (eps / y1)
    f = np.where(y1 > 0, eps / y1, 0.0)

    # Вычисляем правую часть (d^2h/dx^2).
    # Если y1 <= 0, d2h_dx2 должен быть очень большим (штраф),
    # чтобы алгоритм избегал эту область.

    # Расчет чистой правой части для h > 0:
    RHS_pure = (1 / alpha) * (Const + (1 / eps) * (f**3 - f**9))

    # Применяем штраф, если h <= 0
    # Используем np.where, чтобы выбрать RHS_pure, если h>0, или 1e10, если h<=0.
    d2h_dx2 = np.where(y1 > 0, RHS_pure, 1e10)

    # -----------------------------------------------------------

    # Вектор производных [dy1/dx, dy2/dx, dy3/dx]
    # y[1] и y[0] - это тоже массивы, соответствующие y2 и y1
    return np.vstack(
        [y2, d2h_dx2, y1]  # dy1/dx = y2  # dy2/dx = d2h/dx2
    )     # dy3/dx = y1


# --- 3. Функция краевых условий (bc) ---
# ya - решение в точке x=0, yb - решение в точке x=L, p - параметры
# Возвращает вектор невязок размером n+k = 3+1 = 4
def bc(ya, yb, p):
    # p[0] - наш Const

    # 4 условия:
    # 1. y2(0) = 0 (dh/dx(0) = 0)
    # 2. y2(L) = 0 (dh/dx(L) = 0)
    # 3. y3(0) = 0 (Начало интеграла)
    # 4. y3(L) = 1 (Интегральное условие: Integral(h dx) = 1)
    return np.array(
        [
            ya[1],  # ya[1] = y2(0)
            yb[1],  # yb[1] = y2(L)
            ya[2],  # ya[2] = y3(0)
            yb[2] - 1.0,
        ]
    )   # yb[2] = y3(L)


# --- 4. Начальное приближение ---
# Начальное предположение: h(x) = H_avg, dh/dx = 0, Integral(h dx) = H_avg * x
# Так как Integral(h dx) = 1 на интервале L=1, то h_avg должно быть ~ 1/L = 1.0
H_avg = 1.0
x_init = np.linspace(0, L, 100)
y_init = np.zeros((3, x_init.size))
y_init[0] = H_avg * np.ones_like(x_init)     # y1 = h(x) ~ 1
y_init[1] = 0.0 * np.ones_like(x_init)       # y2 = dh/dx ~ 0
y_init[2] = H_avg * x_init                   # y3 = Integral(h dx) ~ x

# Начальное предположение для Const (p)
p_init = [0.0]

# --- 5. Решение BVP ---
print(f'Поиск решения при alpha={alpha}, eps={eps}, L={L}...')
try:
    res = solve_bvp(fun, bc, x_init, y_init, p=p_init)
    print('\n--- Результат ---')
    print(f'Статус решения: {res.message}')
    print(f'Количество узлов: {res.x.size}')

    if res.success:
        # Извлечение найденного параметра Const
        Const_found = res.p[0]
        print(f'\nНайденное значение Const: {Const_found:.6f}')

        h_solution = res.sol(res.x)[0]

        # Проверка интегрального условия
        integral_check = np.trapz(h_solution, res.x)
        print(f'Проверка: Интеграл h(x) dx (численно): {integral_check:.6f}')

        # --- Построение графика ---
        plt.figure(figsize=(10, 5))
        plt.plot(res.x, h_solution, label='$h(x)$')
        plt.title(f'Решение краевой задачи: $\\int_0^L h(x) dx = 1$')
        plt.xlabel('$x$')
        plt.ylabel('$h$')
        plt.grid(True, linestyle='--')
        plt.legend()
        plt.show()

    else:
        print('Решение BVP не было найдено.')

except ValueError as e:
    print(f'Ошибка при решении: {e}')
    print(
        'Возможно, начальное приближение слишком далеко от истинного решения или задача некорректно поставлена.'
    )
