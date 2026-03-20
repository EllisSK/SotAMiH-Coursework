from abc import ABC, abstractmethod

class TemporalIntegrator(ABC):
    @abstractmethod
    def integrate(self, mesh, spatial, dt):
        pass