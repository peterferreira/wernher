from .orbit import Orbit

class Vessel(object):
    def __init__(self,vessel):
        self.orbit = Orbit(vessel.orbit)


