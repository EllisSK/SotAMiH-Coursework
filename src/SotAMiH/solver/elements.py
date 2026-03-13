from __future__ import annotations

class Vector2D:
    def __init__(self, x: int | float, y: int | float) -> None:
        self.x = float(x)
        self.y = float(y)

    def __getitem__(self, index: int):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        else:
            raise IndexError("Vector2D index out of range")
    
    def __setitem__(self, index: int, value: int | float):
        if index == 0:
            self.x = float(value)
        elif index == 1:
            self.y = float(value)
        else:
            raise IndexError("Vector2D index out of range")
    
    def __add__(self, other: Vector2D) -> Vector2D:
        return Vector2D(self.x + other.x, self.y + other.y)
        
    def __sub__(self, other: Vector2D) -> Vector2D:
        return Vector2D(self.x - other.x, self.y - other.y)

    def __mul__(self, other: Vector2D | int | float) -> Vector2D:
        if isinstance(other, Vector2D):
            return Vector2D(self.x * other.x, self.y * other.y)
        else:
            return Vector2D(self.x * other, self.y * other)
        
    def __rmul__(self, other: Vector2D | int | float) -> Vector2D:
        return self.__mul__(other)
        
    def __truediv__(self, other: int | float):
        return Vector2D(self.x / other, self.y / other)
    

class Cell:
    def __init__(self, elevation: float) -> None:
        self.elevation: float = elevation
        self.Q_vector: Vector2D = Vector2D(0, 0)

class Interface:
    def __init__(self) -> None:
        self.F_vector: Vector2D = Vector2D(0, 0)