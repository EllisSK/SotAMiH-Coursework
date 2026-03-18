from abc import ABC, abstractmethod

class TemporalIntegrator(ABC):
    @abstractmethod
    def reconstruct(self, state, mesh):
        pass