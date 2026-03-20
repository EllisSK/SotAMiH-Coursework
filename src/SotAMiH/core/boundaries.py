import numpy as np

from abc import ABC, abstractmethod

class BoundaryCondition:
    @abstractmethod
    def apply(self, interior_slice: np.ndarray, ghost_slice: np.ndarray):
        pass

class ReflectiveBoundary(BoundaryCondition):
    def apply(self, interior_slice: np.ndarray, ghost_slice: np.ndarray, normal_idx: int = 1):
        #Copy the depth and mirror q
        ghost_slice[:] = interior_slice[:]
        ghost_slice[..., normal_idx] = -interior_slice[..., normal_idx]

class TransmissiveBoundary(BoundaryCondition):
    def apply(self, interior_slice: np.ndarray, ghost_slice: np.ndarray):
        #Copy the interior cell values to the ghost cell
        ghost_slice[:] = interior_slice[:]

class FixedFlowBoundary(TransmissiveBoundary):
    def __init__(self, target_q: float):
        self.target_q = target_q

    def apply(self, interior_slice: np.ndarray, ghost_slice: np.ndarray, normal_idx: int = 1):
        #Copy h and set a fixed q in the normal idx direction (future proofing for 2D)
        ghost_slice[:] = interior_slice[:]
        ghost_slice[..., normal_idx] = self.target_q

class FixedDepthBoundary(TransmissiveBoundary):
    def __init__(self, target_h: float):
        self.target_h = target_h

    def apply(self, interior_slice: np.ndarray, ghost_slice: np.ndarray):
        #Copy q and set a fixed depth
        ghost_slice[:] = interior_slice[:]
        ghost_slice[..., 0] = self.target_h