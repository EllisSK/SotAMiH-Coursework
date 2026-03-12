import concurrent.futures
from .elements import Interface, Cell, Vector2D

class Solver:
    def __init__(self) -> None:
        pass

    def calculate_flux(self, l_cell: Cell, r_cell: Cell) -> Vector2D:
        raise NotImplementedError

class HLLSolver(Solver):
    def __init__(self) -> None:
        super().__init__()

    def _calculate_flux_vector(self, h, hu):
        g = 9.81
        dry_tolerance = 1e-6

        if h > dry_tolerance:
            u = hu / h
        else:
            u = 0.0

        F_vector = Vector2D(0, 0)

        F_vector[0] = h * u
        F_vector[1] = (h * (u**2)) + (0.5 * g * (h**2))

        return F_vector


    def _calculate_signal_velocities(self, l_cell: Cell, r_cell: Cell) -> tuple[float, float]:
        g = 9.81
        dry_tolerance = 1e-6

        h_L = l_cell.Q_vector[0]
        h_R = r_cell.Q_vector[0]

        u_L = l_cell.Q_vector[1] / h_L if h_L > dry_tolerance else 0.0
        u_R = r_cell.Q_vector[1] / h_R if h_R > dry_tolerance else 0.0

        s_L = min((u_L - ((g * h_L) ** 0.5)), (u_R - ((g * h_R) ** 0.5)))
        s_R = max((u_L + ((g * h_L) ** 0.5)), (u_R + ((g * h_R) ** 0.5)))

        return s_L, s_R

    def calculate_flux(self, l_cell: Cell, r_cell: Cell):
        S_l, S_r = self._calculate_signal_velocities(l_cell, r_cell)

        if S_l >= 0:
            hll_flux = self._calculate_flux_vector(
                l_cell.Q_vector[0], l_cell.Q_vector[1]
            )
        elif S_r <= 0:
            hll_flux = self._calculate_flux_vector(
                r_cell.Q_vector[0], r_cell.Q_vector[1]
            )
        else:
            F_l = self._calculate_flux_vector(l_cell.Q_vector[0], l_cell.Q_vector[1])
            F_r = self._calculate_flux_vector(r_cell.Q_vector[0], r_cell.Q_vector[1])

            Q_l = l_cell.Q_vector
            Q_r = r_cell.Q_vector

            hll_flux = ((S_r * F_l) - (S_l * F_r) + (S_l * S_r * (Q_r - Q_l))) / (
                S_r - S_l
            )

        return hll_flux