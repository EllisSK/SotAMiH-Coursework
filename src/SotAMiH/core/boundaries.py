from abc import ABC

class BoundaryCondition:
    pass

class ReflectiveBoundary(BoundaryCondition):
    pass

class TransmissiveBoundary(BoundaryCondition):
    pass

class FixedFlowBoundary(TransmissiveBoundary):
    pass

class FixedDepthBoundary(TransmissiveBoundary):
    pass