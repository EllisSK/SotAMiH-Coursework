from .base import SpatialReconstruction

class MUSCL(SpatialReconstruction):
    def __init__(self, limiter_function):
        self.limiter = limiter_function

    def reconstruct(self, state, mesh):
        pass