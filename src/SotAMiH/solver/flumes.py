from .elements import Cell, Interface, Vector2D
from .solvers import Solver, HLLSolver
from .timestepers import calculate_timestep, first_order_time_marching

import matplotlib.pyplot as plt
import numpy as np
import csv

from typing import Callable
from pathlib import Path

class Channel:
    def __init__(self, length: float, width: float, resolution: float, mannings_n: float, solver: Solver, initial_depth_func: Callable | None = None, time_recorders: list[float | int] | None = None, max_time: int | float = np.inf) -> None:
        self.n_cells = int(length / resolution) + 2

        self.length = length
        self.width = width

        self.resolution = resolution
        self.mannings_n = mannings_n
        self.g = 9.81

        self.cells: list[Cell] = []
        self.interfaces: list[Interface] = []

        self.solver = solver
        self.time_recorders = time_recorders
        self.recorded_times = []

        self.max_time = max_time

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

    def solve_flow(self, live_view: bool = False, max_iterations: int = 1000000, min_iterations: int = 100, convergance_definition: float = 10E-9):
        raise NotImplementedError

class VariableBedChannel(Channel):
    def __init__(self, length: float, width: float, resolution: float, mannings_n: float, solver: Solver,  bed_function: Callable, initial_depth_func: Callable | None = None, time_recorders: list[float | int] | None = None, max_time: int | float = np.inf) -> None:
        super().__init__(length, width, resolution, mannings_n, solver, initial_depth_func, time_recorders, max_time)

        self.z_profile = []

        for i in range(1, self.n_cells-1):
            z = bed_function(((i-1) * resolution) + (resolution / 2), length)
            self.z_profile.append(z)

        for i in range(self.n_cells):
            x = ((i-1) * resolution) + (resolution / 2)
            z = bed_function(x, length)
            cell = Cell(z)
            if initial_depth_func:
                cell.Q_vector[0] = initial_depth_func(x)
            self.cells.append(cell)
            if i < self.n_cells-1:
                self.interfaces.append(Interface())

    def _set_boundary_conditions(self):
        self.cells[0].Q_vector[0] = self.cells[1].Q_vector[0]
        self.cells[0].Q_vector[1] = -self.cells[1].Q_vector[1]

        self.cells[-1].Q_vector[0] = self.cells[-2].Q_vector[0]
        self.cells[-1].Q_vector[1] = -self.cells[-2].Q_vector[1]

    def solve_flow(self, live_view: bool = False, max_iterations: int = 1000000, min_iterations: int = 1000, convergance_definition: float = 10E-9):        
        time = 0
        max_h_change = 10E100
        old_profile = self._get_depth_profile()

        if live_view:
            plt.ion()
            fig, ax = plt.subplots()

            line_bed, = ax.plot([], [], 'k-', linewidth=2, label='Bed Profile')
            line_water, = ax.plot([], [], 'b-', linewidth=2, label='Water Surface')

            ax.set_xlabel('Cell Index')
            ax.set_ylabel('Height Above Datum (m)')
            time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, fontsize=12, verticalalignment='top')
            ax.set_ylim(0, 0.8)
            
            min_z = min(self.z_profile)
            max_z = max(self.z_profile)

            ax.set_ylim(min_z - 0.1, max_z+20) 
            ax.legend(loc='upper right')

            x_data = list(range(len(old_profile)))
            line_bed.set_data(x_data, self.z_profile)
            line_water.set_xdata(x_data)

        if self.time_recorders:
            times_to_capture = self.time_recorders.copy()
        else:
            times_to_capture = None
        record_step = False

        iterations = 0
        while (time <= self.max_time) and ((iterations < min_iterations) or ((iterations < max_iterations) and (max_h_change > convergance_definition))):
            timestep = calculate_timestep(self.cells, self.resolution)

            if time + timestep > self.max_time:
                timestep = self.max_time - time

            if (times_to_capture) and (time + timestep >= times_to_capture[0]):
                timestep = times_to_capture[0] - time
                record_step = True
                times_to_capture.pop(0)
            else:
                record_step = False
            
            self._set_boundary_conditions()

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

            if live_view :
                water_surface = [h + z for h, z in zip(new_profile, self.z_profile)]
                line_water.set_ydata(water_surface)
                
                ax.relim()
                ax.autoscale_view(scaley=False)
                
                fig.canvas.draw()
                fig.canvas.flush_events()

                time_text.set_text(f"t = {time:.2f}")

            velocity_profile = self._get_velocity_profile()

            if record_step:
                self.recorded_times.append((new_profile, velocity_profile))

            iterations += 1

        if live_view:
            plt.ioff()
            plt.show()
    
    def write_results(self, path: Path):
            path.parent.mkdir(parents=True, exist_ok=True)
            
            x_coords = [i * self.resolution + (self.resolution / 2) for i in range(len(self.z_profile))]
            
            with open(path, 'w', newline='') as f:
                writer = csv.writer(f)
                
                writer.writerow(['# Channel Configuration'])
                writer.writerow([
                    f'# Length: {self.length}', 
                    f'Width: {self.width}', 
                    f'Resolution: {self.resolution}', 
                    f'Mannings_n: {self.mannings_n}'
                ])
                
                writer.writerow(['Time', 'Distance_x', 'Bed_z', 'Depth_h', 'Velocity_u'])
                
                if self.time_recorders:
                    for i, (depths, velocities) in enumerate(self.recorded_times):
                        t = self.time_recorders[i]
                        for x, z, h, u in zip(x_coords, self.z_profile, depths, velocities):
                            writer.writerow([t, x, z, h, u])