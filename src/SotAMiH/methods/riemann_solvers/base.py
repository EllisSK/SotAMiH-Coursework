from abc import ABC, abstractmethod

class RiemannSolver(ABC):
    @abstractmethod
    def solve(self, q_L, q_R, physics):
        """Returns the numerical flux at the interface."""
        pass