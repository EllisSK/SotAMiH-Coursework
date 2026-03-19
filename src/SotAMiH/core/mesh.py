from abc import ABC

class Mesh(ABC):
    pass

class Mesh1D(Mesh):
    def __init__(self, length: float, resolution: float) -> None:
        super().__init__()
        self.length = length
        self.dx = resolution
        self.N = length / resolution