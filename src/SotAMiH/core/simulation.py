import numpy as np

from .mesh import Mesh
from.boundaries import BoundaryCondition

from SotAMiH.methods.spatial import SpatialReconstruction
from SotAMiH.methods.temporal import TemporalIntegrator
from SotAMiH.methods.riemann_solvers import RiemannSolver
from SotAMiH.physics import Physics

class Simulation:
    def __init__(self, mesh: Mesh, physics: Physics, spatial: SpatialReconstruction, temporal: TemporalIntegrator, riemann: RiemannSolver, bcs: dict[str, BoundaryCondition]):
        self.mesh = mesh
        self.physics = physics
        self.spatial = spatial
        self.temporal = temporal
        self.riemann = riemann
        self.bcs = bcs

    def run(self, end_time: int | float | None = None, convergance_threshold: int | float | None = None, record_times: list[int | float] | None = None):
        if end_time and convergance_threshold:
            raise ValueError("Please provide only 1 end condition.")
        elif not end_time and not convergance_threshold:
            raise ValueError("Please provide 1 end condition.")
        
        self.t = 0
        self.max_change = 0
        
        record_time = False

        if not end_time:
            end_time = np.inf

        if not convergance_threshold:
            convergance_threshold = np.inf

        if record_times:
            self.ttr = record_times.copy()

        while (self.t < end_time) and (self.max_change < convergance_threshold):
            dt = self.physics.dynamic_timestep(self.mesh.Q_array)

            if (self.t + dt) > end_time:
                dt = end_time - self.t
                record_time = True
            elif (self.t + dt) > self.ttr[0]:
                dt = self.ttr[0] - self.t
                record_time = True
                self.ttr.pop(0)
            else:
                record_time = False

            self.mesh.apply_boundary_conditions(self.bcs)

            self.temporal.integrate(self.mesh, self.spatial, dt)

            self.t += dt

            if record_time:
                pass