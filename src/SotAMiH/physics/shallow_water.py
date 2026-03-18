import numpy as np

from.base import Physics

class ShallowWater1D(Physics):
    def __init__(self, dx):
        self.g = 9.81
        self.dx = dx

    def flux(self, h, q):
        u = np.divide(q, h, where=h > 0)

        return h * u, (np.power(q, 2) / h) + (0.5 * self.g * np.power(h, 2))
        
    def max_wave_speed(self, h, q):
        u = np.divide(q, h, where=h > 0)
        a = np.sqrt(self.g * h)
        
        return np.max(np.abs(u) + a)
    
    def dynamic_timestep(self, h, q):
        max_speed = self.max_wave_speed(h, q)
        return self.dx / self.max_wave_speed