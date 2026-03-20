import numpy as np

from abc import ABC, abstractmethod
from .boundaries import BoundaryCondition

class Mesh(ABC):
    @abstractmethod
    def apply_boundary_conditions(self, boundary_conditions: dict[str, BoundaryCondition]):
        pass

class Mesh1D(Mesh):
    def __init__(self, length: float, resolution: float) -> None:
        super().__init__()
        self.length = length
        self.dx = resolution
        self.N = int(length / resolution)

        #Create 2D array that is as long as the domain + 2 ghost cells
        self.Q_array = np.zeros((self.N+2, 2))
        self.F_array = np.zeros((self.N+2, 2))

    def apply_boundary_conditions(self, boundary_conditions: dict[str, BoundaryCondition]):
        lb = boundary_conditions["left_boundary"]
        rb = boundary_conditions["right_boundary"]

        lb.apply(
            interior_slice=self.Q_array[1], 
            ghost_slice=self.Q_array[0], 
        )

        rb.apply(
            interior_slice=self.Q_array[-2], 
            ghost_slice=self.Q_array[-1], 
        )