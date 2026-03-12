from .elements import Cell, Interface, Vector2D
from .solvers import Solver, HLLSolver
from .timestepers import calculate_timestep, first_order_time_marching

import matplotlib.pyplot as plt
import time as t

from typing import Callable

class Flume:
    def __init__(self, length: float, width: float, resolution: float, mannings_n: float, solver: Solver) -> None:
        self.n_cells = int(length / resolution) + 2

        self.length = length
        self.width = width

        self.resolution = resolution
        self.mannings_n = mannings_n
        self.g = 9.81

        self.cells: list[Cell] = []
        self.interfaces: list[Interface] = []

        self.solver = solver

    def _calculate_source_vector(self, cell: Cell, l_cell: Cell, r_cell: Cell):
        g = self.g
        h = cell.Q_vector[0]

        dry_tolerance = 1e-6

        if h > dry_tolerance:
            u = cell.Q_vector[1] / h
            S_fric = - (g * (self.mannings_n ** 2) * abs(u) * u) / (h ** (1/3))
        else:
            u = 0.0
            S_fric = 0.0


        S_vector = Vector2D(0, 0)

        S_slope = -g * h * ((r_cell.elevation - l_cell.elevation) / (2 * self.resolution))

        S_vector[1] = S_slope + S_fric

        return S_vector
    
    def _get_depth_profile(self):
        n_cells = len(self.cells)
        profile = []
        for i in range(n_cells):
            if (i != 0) and (i != n_cells-1):
                profile.append(self.cells[i].Q_vector[0])
        return profile
    
    def _get_velocity_profile(self):
        n_cells = len(self.cells)
        profile = []
        for i in range(n_cells):
            if (i != 0) and (i != n_cells-1):
                h = self.cells[i].Q_vector[0]
                if h != 0:
                    u = self.cells[i].Q_vector[1] / h
                else:
                    u = 0
                profile.append(u)
        return profile

    def solve_flow(self, flow: float, live_view: bool = False, max_iterations: int = 1000000, min_iterations: int = 100, convergance_definition: float = 10E-9):
        raise NotImplementedError

class EmptyFlume(Flume):
    def __init__(self, length: float, width: float, resolution: float, mannings_n: float, solver: Solver) -> None:
        super().__init__(length, width, resolution, mannings_n, solver)

        for i in range(self.n_cells):
            self.cells.append(Cell(0))
            if i < self.n_cells-1:
                self.interfaces.append(Interface())

    def _set_boundary_conditions(self, q: float):
        self.cells[0].Q_vector[0] = self.cells[0].Q_vector[1]
        self.cells[0].Q_vector[1] = q

        self.cells[-1].Q_vector[0] = self.cells[-2].Q_vector[0]
        self.cells[-1].Q_vector[1] = self.cells[-2].Q_vector[1]

    def solve_flow(self, flow: float, live_view: bool = False, max_iterations: int = 1000000, min_iterations: int = 1000, convergance_definition: float = 10E-9):
        unit_flow = flow / self.width
        
        time = 0
        max_h_change = 10E100
        old_profile = self._get_depth_profile()

        if live_view:
            plt.ion()
            fig, ax = plt.subplots()

            line, = ax.plot([], [], 'b-') 
            ax.set_xlabel('Cell Index')
            ax.set_ylabel('Depth Profile (m)')
            ax.set_title(f'Empty Flume Flow. Set Flow: {flow*1000} l/s')
            time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, fontsize=12, verticalalignment='top')
            ax.set_ylim(0, 0.8)
            line.set_xdata(list(range(len(old_profile))))

        iterations = 0
        while (iterations < min_iterations) or ((iterations < max_iterations) and (max_h_change > convergance_definition)):
            timestep = calculate_timestep(self.cells, self.resolution)
            self._set_boundary_conditions(unit_flow)

            for i in range(len(self.interfaces)):
                interface = self.interfaces[i]
                l_cell = self.cells[i]
                r_cell = self.cells[i+1]
                interface.F_vector = self.solver.calculate_flux(l_cell, r_cell)

            for i in range(1, len(self.cells)-1):
                cell = self.cells[i]
                l_cell = self.cells[i-1]
                r_cell = self.cells[i+1]
                
                S_vector = self._calculate_source_vector(cell, l_cell, r_cell)
                F_left = self.interfaces[i-1].F_vector
                F_right = self.interfaces[i].F_vector

                cell.Q_vector = first_order_time_marching(F_left, cell.Q_vector, F_right, S_vector, timestep, self.resolution)

            new_profile = self._get_depth_profile()
            max_h_change = max(abs(new - old) for new, old in zip(new_profile, old_profile))
            old_profile = new_profile

            if live_view :
                line.set_ydata(new_profile)
                
                ax.relim()
                ax.autoscale_view(scaley=False)
                
                fig.canvas.draw()
                fig.canvas.flush_events()

                time_text.set_text(f"t = {time:.2f}")

            time += timestep
            iterations += 1

        return self._get_depth_profile()

class WeirFlume(Flume):
    def __init__(self, length: float, width: float, resolution: float, mannings_n: float, solver: Solver, weir_psoition: float) -> None:
        if weir_psoition >= length - resolution:
            raise ValueError("The weir is out of bounds!")
        
        super().__init__(length, width, resolution, mannings_n, solver)

    def _set_boundary_conditions(self, q: float):
        self.cells[0].Q_vector[0] = self.cells[0].Q_vector[1]
        self.cells[0].Q_vector[1] = q

        self.cells[-1].Q_vector[0] = self.cells[-2].Q_vector[0]
        self.cells[-1].Q_vector[1] = self.cells[-2].Q_vector[1]

    def solve_flow(self, flow: float, live_view: bool = False, max_iterations: int = 1000000, min_iterations: int = 1000, convergance_definition: float = 10E-9):
        unit_flow = flow / self.width
        
        time = 0
        max_h_change = 10E100
        old_profile = self._get_depth_profile()

        iterations = 0
        while (iterations < min_iterations) or ((iterations < max_iterations) and (max_h_change > convergance_definition)):
            print(iterations)
            
            timestep = calculate_timestep(self.cells, self.resolution)
            self._set_boundary_conditions(unit_flow)

            for i in range(len(self.interfaces)):
                interface = self.interfaces[i]
                l_cell = self.cells[i]
                r_cell = self.cells[i+1]
                interface.F_vector = self.solver.calculate_flux(l_cell, r_cell)

            for i in range(1, len(self.cells)-1):
                cell = self.cells[i]
                l_cell = self.cells[i-1]
                r_cell = self.cells[i+1]
                
                S_vector = self._calculate_source_vector(cell, l_cell, r_cell)
                F_left = self.interfaces[i-1].F_vector
                F_right = self.interfaces[i].F_vector

                cell.Q_vector = first_order_time_marching(F_left, cell.Q_vector, F_right, S_vector, timestep, self.resolution)

            new_profile = self._get_depth_profile()
            max_h_change = max(abs(new - old) for new, old in zip(new_profile, old_profile))
            old_profile = new_profile

            time += timestep
            iterations += 1

        return self._get_depth_profile()

class SluiceFlume(Flume):
    def __init__(self, length: float, width: float, resolution: float, mannings_n: float, solver: Solver) -> None:
        super().__init__(length, width, resolution, mannings_n, solver)

class LeakyBarrierFlume(Flume):
    def __init__(self, length: float, width: float, resolution: float, mannings_n: float, solver: Solver) -> None:
        super().__init__(length, width, resolution, mannings_n, solver)

class VariableBedlume(Flume):
    def __init__(self, length: float, width: float, resolution: float, mannings_n: float, solver: Solver, bed_function: Callable) -> None:
        super().__init__(length, width, resolution, mannings_n, solver)

        self.z_profile = []

        for i in range(1, self.n_cells-1):
            z = bed_function(((i-1) * resolution) + (resolution / 2), length)
            self.z_profile.append(z)

        for i in range(self.n_cells):
            z = bed_function(((i-1) * resolution) + (resolution / 2), length)
            self.cells.append(Cell(z))
            if i < self.n_cells-1:
                self.interfaces.append(Interface())

    def _set_boundary_conditions(self, q: float):
        self.cells[0].Q_vector[0] = self.cells[1].Q_vector[0]
        self.cells[0].Q_vector[1] = -self.cells[1].Q_vector[1]

        self.cells[-1].Q_vector[0] = self.cells[-2].Q_vector[0]
        self.cells[-1].Q_vector[1] = -self.cells[-2].Q_vector[1]

    def _cell_update(self, cell):
        pass

    def solve_flow(self, flow: float, live_view: bool = False, max_iterations: int = 1000000, min_iterations: int = 1000, convergance_definition: float = 10E-9):
        unit_flow = flow / self.width
        
        time = 0
        max_h_change = 10E100
        old_profile = self._get_depth_profile()

        if live_view:
            plt.ion()
            fig, ax = plt.subplots()

            line_bed, = ax.plot([], [], 'k-', linewidth=2, label='Bed Profile')
            line_water, = ax.plot([], [], 'b-', linewidth=2, label='Water Surface')

            ax.set_xlabel('Cell Index')
            ax.set_ylabel('Depth Profile (m)')
            ax.set_title(f'Empty Flume Flow. Set Flow: {flow*1000} l/s')
            time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, fontsize=12, verticalalignment='top')
            ax.set_ylim(0, 0.8)
            
            min_z = min(self.z_profile)
            max_z = max(self.z_profile)
            ax.set_ylim(min_z - 0.1, max_z + 65) 
            ax.legend(loc='upper right')

            x_data = list(range(len(old_profile)))
            line_bed.set_data(x_data, self.z_profile)
            line_water.set_xdata(x_data)

        iterations = 0
        while (iterations < min_iterations) or ((iterations < max_iterations) and (max_h_change > convergance_definition)):
            timestep = calculate_timestep(self.cells, self.resolution)
            self._set_boundary_conditions(unit_flow)

            for i in range(len(self.interfaces)):
                interface = self.interfaces[i]
                l_cell = self.cells[i]
                r_cell = self.cells[i+1]
                interface.F_vector = self.solver.calculate_flux(l_cell, r_cell)

            for i in range(1, len(self.cells)-1):
                cell = self.cells[i]
                l_cell = self.cells[i-1]
                r_cell = self.cells[i+1]
                
                S_vector = self._calculate_source_vector(cell, l_cell, r_cell)
                F_left = self.interfaces[i-1].F_vector
                F_right = self.interfaces[i].F_vector

                cell.Q_vector = first_order_time_marching(F_left, cell.Q_vector, F_right, S_vector, timestep, self.resolution)

            new_profile = self._get_depth_profile()
            max_h_change = max(abs(new - old) for new, old in zip(new_profile, old_profile))
            old_profile = new_profile

            if live_view :
                water_surface = [h + z for h, z in zip(new_profile, self.z_profile)]
                line_water.set_ydata(water_surface)
                
                ax.relim()
                ax.autoscale_view(scaley=False)
                
                fig.canvas.draw()
                fig.canvas.flush_events()

                time_text.set_text(f"t = {time:.2f}")

            time += timestep
            iterations += 1

            if time == 5000:
                print("Depth Profile:", new_profile)
                print("Velocity Profile:", self._get_velocity_profile())

        return self._get_depth_profile()