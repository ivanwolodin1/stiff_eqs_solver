from pde import PDEBase, FieldCollection
from numpy import sqrt, tan, real, sin, sinh, cos, cosh


class Volodin2025(PDEBase):
    def __init__(self, Ca, G0, b, omega, phi, Ma, Bi, Pr, bc):
        super().__init__()
        self.Ca = Ca
        self.G0 = G0
        self.b = b
        self.omega = omega
        self.phi = phi
        self.Ma = Ma
        self.Bi = Bi
        self.Pr = Pr
        self.bc = bc
        self.alpha = -(1j - 1) * (sqrt(2.0 * self.omega)) / 2.0

    def evolution_rate(self, state, t=0):
        h, T = state
        f = T - h
        # https://py-pde.readthedocs.io/en/latest/examples_gallery/advanced_pdes/pde_1d_class.html
        # on gradient() operator

        dx_h = h.gradient(bc=self.bc)[0]

        # second derivatives
        d2x_h = h.laplace(bc=self.bc)

        # --- Π ---
        f_alpha = 1.0 - (tan(self.alpha * h)) / (self.alpha * h)
        real_f_alpha = real(f_alpha)
        # real_f_alpha = self.real_alpha

        H = -(h * real_f_alpha * dx_h).gradient(bc=self.bc)[0]

        Pi = (
            self.phi / h**3
            - self.Ca * d2x_h
            + self.G0 * h
            + 0.5 * self.b**2 * self.omega**2 * H
        )

        # derivative Pi_x
        dx_Pi = Pi.gradient(bc=self.bc)[0]

        # --- Q ---
        gamma = sqrt(2 * self.omega) * h
        ksi_plus = cosh(gamma) + cos(gamma)
        ksi_minus = cosh(gamma) - cos(gamma)
        phi_plus = sinh(gamma) + sin(gamma)
        phi_minus = sinh(gamma) - sin(gamma)

        Q1 = 3.0 * (
            (
                4.0 * sinh(gamma) * sin(gamma)
                - gamma * (ksi_plus * phi_plus - ksi_minus * phi_minus)
            )
            / (2 * gamma * gamma * ksi_plus * ksi_plus)
        )

        Q2 = -1.0 / 3.0 + (11.0 * phi_minus - 3.0 * gamma * ksi_minus) / (
            gamma**3 * ksi_plus
        )

        Q = Q1 * h**2 * (dx_h) ** 3 + Q2 * h**3 * dx_h * d2x_h

        h_t = (
            (1.0 / 3.0) * h**3 * dx_Pi
            + 0.5 * self.Ma * h**2 * f.gradient(bc=self.bc)[0]
            - 0.5 * self.b**2 * self.omega**2 * Q
        ).gradient(bc=self.bc)[0]

        dx_T = T.gradient(bc=self.bc)[0]
        dx_f = f.gradient(bc=self.bc)[0]

        # (1/Pr) * d/dx ( h * dT/dx )
        term1 = (h * dx_T / self.Pr).gradient(bc=self.bc)[0]

        # - (1/(2Pr)) * (dh/dx)^2
        term2 = -(dx_h**2) / (2 * self.Pr)

        # - beta/Pr * f
        term3 = -(self.Bi / self.Pr) * f

        # ( h^3/3 * dPi/dx + Ma/2 * h^2 * df/dx ) * df/dx
        term4 = ((h**3 / 3) * dx_Pi + 0.5 * self.Ma * h**2 * dx_f) * dx_f

        # term5: d/dx( h^4/8 * dPi/dx + Ma/6 * h^3 * df/dx )
        inside = (h**4 / 8) * dx_Pi + (self.Ma / 6) * h**3 * dx_f
        term5 = inside.gradient(bc=self.bc)[0]

        # h * T_t = term1 + term2 + term3 + term4 + term5
        T_t = (term1 + term2 + term3 + term4 + term5) / h

        return FieldCollection([h_t, T_t])
