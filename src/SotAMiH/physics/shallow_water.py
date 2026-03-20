import numpy as np

from.base import Physics

class ShallowWater1D(Physics):
    def __init__(self, dx):
        self.g = 9.81
        self.dx = dx

    def flux(self, Q_array):
        h = Q_array[:, 0]
        q = Q_array[:, 1]
        u = np.divide(q, h, where=h > 0)
        
        F_array = np.empty_like(Q_array)
        F_array[:, 0] = h * u
        F_array[:, 1] = (np.power(q, 2) / h) + (0.5 * self.g * np.power(h, 2))
        
        return F_array
        
    def max_wave_speed(self, Q_array):
        h = Q_array[:, 0]
        q = Q_array[:, 1]
        u = np.divide(q, h, where=h > 0)
        a = np.sqrt(self.g * h)
        return np.max(np.abs(u) + a)
    
    def dynamic_timestep(self, Q_array) -> float:
        max_speed = self.max_wave_speed(Q_array)
        #Shits result down 1 bit ensuring fp errors never take Cr over 1
        return np.nextafter(self.dx / max_speed, -np.inf)