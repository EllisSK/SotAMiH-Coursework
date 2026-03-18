from abc import ABC, abstractmethod

class SpatialReconstruction(ABC):
    @abstractmethod
    def reconstruct(self, state, mesh):
        """Returns the Left (q_L) and Right (q_R) states at cell interfaces."""
        pass