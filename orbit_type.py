from enum import Enum

class OrbitType(Enum):
    circular = 0
    elliptic = 1
    parabolic = 2
    hyperbolic = 3

    def __init__(self,val):
        self.val = val

    @property
    def isclosed(self):
        return self.val in [0,1]

    @property
    def isopen(self):
        return self.val in [2,3]

    @staticmethod
    def from_eccentricity(e,δe=1e-5):
        if e < δe:
            return OrbitType.circular
        elif e < (1 - δe):
            return OrbitType.elliptic
        elif e < (1 + δe):
            return OrbitType.parabolic
        else:
            return OrbitType.hyperbolic

circular = OrbitType.circular
elliptic = OrbitType.elliptic
hyperbolic = OrbitType.hyperbolic
parabolic = OrbitType.parabolic
