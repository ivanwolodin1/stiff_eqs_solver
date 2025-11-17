import numpy as np

h = np.linspace(0.1, 2.0, 10)  # 10 elements from 0.1 to 2.0
omega = 2.0

alpha = (1j + 1) * (np.sqrt(2.0 * omega)) / 2.0

value = 1 - (np.tan(alpha * h)) / (alpha * h)

# real part
real_part = np.real(value)

print(real_part)
print(np.imag(1j + 1))
