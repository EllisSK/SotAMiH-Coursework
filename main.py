import numpy as np

from pathlib import Path

from src.SotAMiH.solver.flumes import VariableBedChannel
from src.SotAMiH.solver.solvers import HLLSolver

def test_case_1():
    def bed_fn(x, L):
        return (10) + (40 * x / L) + (10 * np.sin(np.pi * ((4 * x / L) - (1/2))))
    
    def initial_depth(x):
        return 65 -((10) + (40 * x / 14000) + (10 * np.sin(np.pi * ((4 * x / 14000) - (1/2)))))
    
    channel = VariableBedChannel(14000, 1, 280, 0, HLLSolver(), bed_fn, initial_depth_func= initial_depth, time_recorders=[1000, 2000, 3000, 4000, 5000], max_time=5000)

    channel.solve_flow(live_view= True)

    channel.write_results(Path("exports/testcase1/results.csv"))

def test_case_2():
    def bed_fn(x, L):
        return -0.01*x
    
    def initial_depth(x):
        if x < 10:
            return 2
        else:
            return 10e-6
    
    channel = VariableBedChannel(140, 1, 0.5, 0.01, HLLSolver(), bed_function= bed_fn, initial_depth_func=initial_depth)

    channel.solve_flow(live_view= True)

def test_case_3():
    pass

def test_case_4():
    pass

def main():
    test_case_1()


if __name__ == "__main__":
    main()
