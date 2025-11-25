from pde import PDEBase
from numpy import sqrt, tan, real, sin, sinh, cos, cosh


class Shklyaev2008(PDEBase):
    """https://doi.org/10.1103/PhysRevE.77.036320
    Solution to 30(a)-31(b) eq.
    """

    def __init__(
        self,
        Ca=0,
        G0=0,
        b=0,
        omega=10,
        phi=0,
        bc={'x': 'periodic'},
    ):
        super().__init__()
        self.bc = bc
        self.Ca = Ca
        self.G0 = G0
        self.b = b
        self.omega = omega
        self.phi = phi
        # self.alpha = (1j - 1) * (sqrt(2.0 * self.omega)) / 2.0
        self.alpha = -(1j - 1) * (sqrt(2.0 * self.omega)) / 2.0

        # self.real_alpha = 1.0 - (sqrt(2) / (2 * sqrt(self.omega))) * (
        #     (
        #         sin(sqrt(2) * sqrt(self.omega) / 2)
        #         * cos(sqrt(2) * sqrt(self.omega) / 2)
        #         + sinh(sqrt(2) * sqrt(self.omega) / 2)
        #         * cosh(sqrt(2) * sqrt(self.omega) / 2)
        #     )
        #     / (
        #         cos(sqrt(2) * sqrt(self.omega) / 2) ** 2
        #         + cosh(sqrt(2) * sqrt(self.omega) / 2) ** 2
        #         - 1
        #     )
        # )

    def evolution_rate(self, h, t=0):
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
        eps = 10000000000
        Pi = (
            self.phi / h**3
            - self.Ca * eps * eps * d2x_h
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
            - 0.5 * self.b**2 * self.omega**2 * Q
        ).gradient(bc=self.bc)[0]

        return h_t


class Shklyaev2008_VIII_B(PDEBase):
    """https://doi.org/10.1103/PhysRevE.77.036320
    Solution to 30(a)-31(b) eq.
    """

    def __init__(
        self,
        V,
        k,
        G0,
        omega,
        bc,
    ):
        super().__init__()
        self.bc = bc
        self.V = V
        self.k = k
        self.G0 = G0
        self.omega = omega
        self.alpha = -(1j - 1) * (sqrt(2.0 * self.omega)) / 2.0

    def evolution_rate(self, h, t=0):
        # https://py-pde.readthedocs.io/en/latest/examples_gallery/advanced_pdes/pde_1d_class.html
        # on gradient() operator

        dx_h = h.gradient(bc=self.bc)[0]

        # second derivatives
        d2x_h = h.laplace(bc=self.bc)

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

        # --- Π ---
        f_alpha = 1.0 - (tan(self.alpha * h)) / (self.alpha * h)
        real_f_alpha = real(f_alpha)
        # real_f_alpha = self.real_alpha

        H = -(h * real_f_alpha * dx_h).gradient(bc=self.bc)[0]
        P = (
            1.0 / (3.0 * h**3)
            + self.G0 * h
            + self.k * self.k * (self.V * H - d2x_h)
        )

        dx_P = P.gradient(bc=self.bc)[0]

        h_t = (
            self.k
            * self.k
            * (h**3 * dx_P - 3.0 * self.V * self.k * self.k * Q).gradient(
                bc=self.bc
            )[0]
        )

        return h_t
