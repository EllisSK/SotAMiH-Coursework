from .mesh import Mesh
from .state import StateContainer
from.boundaries import BoundaryCondition


from SotAMiH.methods.spatial import SpatialReconstruction
from SotAMiH.methods.temporal import TemporalIntegrator
from SotAMiH.methods.riemann_solvers import RiemannSolver
from SotAMiH.physics import Physics

class Simulation:
    def __init__(self, mesh: Mesh, state: StateContainer, physics: Physics, spatial: SpatialReconstruction, temporal: TemporalIntegrator, riemann: RiemannSolver, bcs: BoundaryCondition):
        self.mesh = mesh
        self.state = state
        self.physics = physics
        self.spatial = spatial
        self.temporal = temporal
        self.riemann = riemann
        self.bcs = bcs