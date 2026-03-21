import numpy as np

from pathlib import Path

from src.SotAMiH.core.simulation import Simulation
from src.SotAMiH.core.mesh import Mesh1D
from src.SotAMiH.physics.shallow_water import ShallowWater1D
from src.SotAMiH.methods.spatial.MUSCL import MUSCL1D
from src.SotAMiH.methods.temporal.range_kutta import RK2
from src.SotAMiH.methods.riemann_solvers.hll import HLLSolver
from src.SotAMiH.core.boundaries import ReflectiveBoundary, VariableConservedBoundary, TransmissiveBoundary, VariableDepthBoundary

def test_case_1():
    def bed_fn(x):
        return (10) + (40 * x / 14000) + (10 * np.sin(np.pi * ((4 * x / 14000) - (1/2))))
    
    def initial_cond(x):
        return 65, 0
    
    mesh = Mesh1D(14000, 280, initial_cond, bed_fn)
    physics = ShallowWater1D(mesh.dx)
    spatial = MUSCL1D()
    temporal = RK2()
    riemann = HLLSolver()

    bcs = {
        "left_boundary" : ReflectiveBoundary(),
        "right_boundary" : ReflectiveBoundary()
    }
    
    sim = Simulation(mesh, physics, spatial, temporal, riemann, bcs)

    sim.run(5000)


def test_case_2():
    def bed_fn(x):
        return (10) + (40 * x / 14000) + (10 * np.sin(np.pi * ((4 * x / 14000) - (1/2))))
    
    def initial_cond(x):      
        return 60.5, 0
    
    def b_ht(t):
        return 64.5 - (4 * np.sin(np.pi * (((4 * t) / 86400) + (0.5))))

    def lb_qt(t):
        def bed_fn(x):
            return (10) + (40 * x / 14000) + (10 * np.sin(np.pi * ((4 * x / 14000) - (1/2))))

        def b_ht(t):
            return 64.5 - (4 * np.sin(np.pi * (((4 * t) / 86400) + (0.5))))

        return (b_ht(0) - bed_fn(0)) * (((np.pi * (0-14000)) / (5400 * b_ht(t))) * np.cos(np.pi * (((4 * t) / 86400) + (0.5))))

    def rb_qt(t):
        def bed_fn(x):
            return (10) + (40 * x / 14000) + (10 * np.sin(np.pi * ((4 * x / 14000) - (1/2))))

        def b_ht(t):
            return 64.5 - (4 * np.sin(np.pi * (((4 * t) / 86400) + (0.5))))

        return (b_ht(14000) - bed_fn(14000)) * (((np.pi * (14000-14000)) / (5400 * b_ht(t))) * np.cos(np.pi * (((4 * t) / 86400) + (0.5))))
    
    mesh = Mesh1D(14000, 280, initial_cond, bed_fn)
    physics = ShallowWater1D(mesh.dx)
    spatial = MUSCL1D()
    temporal = RK2()
    riemann = HLLSolver()

    bcs = {
        "left_boundary" : VariableConservedBoundary(b_ht, lb_qt),
        "right_boundary" : VariableConservedBoundary(b_ht, rb_qt)
    }
    
    sim = Simulation(mesh, physics, spatial, temporal, riemann, bcs)

    sim.run(7552.13)

def test_case_3():
    def bed_fn(x):
        return np.zeros_like(x)
    
    def initial_cond(x):
        Q = np.zeros((len(x), 2), dtype=float)
        Q[:, 0] = np.where(x < 25, 1.0, 0.1)
        return Q
    
    mesh = Mesh1D(50, 0.25, initial_cond, bed_fn)
    physics = ShallowWater1D(mesh.dx)
    spatial = MUSCL1D()
    temporal = RK2()
    riemann = HLLSolver()

    bcs = {
        "left_boundary" : TransmissiveBoundary(),
        "right_boundary" : TransmissiveBoundary()
    }
    
    sim = Simulation(mesh, physics, spatial, temporal, riemann, bcs)

    sim.run(5)

def test_case_4():
    def bed_fn(x):
        zb = np.zeros_like(x)
        zb = np.where(abs(x-750) < 187.5, 8.0, 0.0)
        return zb
    
    def b_ht(t):
        return 20 - (4 * np.sin(np.pi * (((4 * t) / 86400) + (0.5))))
    
    def initial_cond(x):      
        return 16, 0
    
    mesh = Mesh1D(1500, 7.5, initial_cond, bed_fn)
    physics = ShallowWater1D(mesh.dx)
    spatial = MUSCL1D()
    temporal = RK2()
    riemann = HLLSolver()

    bcs = {
        "left_boundary" : VariableDepthBoundary(b_ht),
        "right_boundary" : ReflectiveBoundary()
    }
    
    sim = Simulation(mesh, physics, spatial, temporal, riemann, bcs)

    sim.run(32400, record_times=[10800, 32400])

def main():
    test_case_4()


if __name__ == "__main__":
    main()
