from abc import ABC, abstractmethod

class RiemannSolver(ABC):
    @abstractmethod
    def solve(self):
        pass