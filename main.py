import numpy as np

from src.SotAMiH.solver.flumes import VariableBedlume
from src.SotAMiH.solver.solvers import HLLSolver

def test_case_1():
    def bed_fn(x, L):
        return (10) + (40 * x / L) + (10 * np.sin(np.pi * ((4 * x / L) - (1/2))))
    
    flume = VariableBedlume(14000, 1, 280, 0, HLLSolver(), bed_fn)

    flume.solve_flow(0, True)

def test_case_2():
    def bed_fn(x, L):
        return (10) + (40 * x / L) + (10 * np.sin(np.pi * ((4 * x / L) - (1/2))))
    
    flume = VariableBedlume(14000, 1, 280, 0, HLLSolver(), bed_fn)

    flume.solve_flow(0, True)

def test_case_3():
    pass

def test_case_4():
    pass

def main():
    test_case_1()


if __name__ == "__main__":
    main()
