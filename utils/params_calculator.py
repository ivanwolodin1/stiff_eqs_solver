H = 10**-6
g = 9.81
rho = 1000
nu = 10**-6
khi = 0.14 * 10**-6
kappa = 10**-1
q = 10
A_prime = 6 * 3.14 * 10**-21
sigma = 72 * 10**-3
d_sigma_dT = -(10**-4)
dT_dz = 1 / H

frequency = 2 * 3.14 * 10**5
amplitude = H * 10  # amplitude is large


Ma = -(d_sigma_dT) * abs(dT_dz) / (6 * rho * nu**2) * H**2
Ga = g * H**3 / nu**2
Bi = q * H / kappa
Pr = nu / khi
Ca = sigma * H / (rho * nu**2)
A = A_prime / (6 * 3.14 * rho * nu**2 * H)
omega = frequency * H**2 / nu
b = amplitude / H

eps = 10**-2
print(
    f'{Ma=:.3e}, {Ga=:.3e}, {Bi=:.3e}, Pr={Pr}, Ca={Ca}, {A=:.3e}, b={b:.3e}, omega={omega:.3e}'
)

print('=================================================================')
print('After eps scaling:')

Ca = Ca * eps**2
Bi = Bi / eps**2

print(
    f'{Ma=:.3e}, {Ga=:.3e}, {Bi=:.3e}, Pr={Pr}, Ca={Ca}, {A=:.3e}, b={b:.3e}, omega={omega:.3e}'
)

print(f'V={(b**2*omega**2)/(2*Ca)}')
