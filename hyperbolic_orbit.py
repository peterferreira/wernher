from .base_orbit import *
from .elliptic_orbit import EllipticOrbit

class HyperbolicOrbit(EllipticOrbit,BaseOrbit):


    def __init__(self,orbit):
        #super().__init__(orbit)
        BaseOrbit.__init__(self,orbit)


    ### Specific Angular Momentum (h) vs Semi-Major Axis (a)
    @property
    def specific_angular_momentum(self):
        if not hasattr(self,'_specific_angular_momentum'):
            e = self.eccentricity
            a = self.semi_major_axis
            μ = self.body.gravitational_parameter
            h = sqrt(a * μ * (e**2 - 1))
            self._specific_angular_momentum = h
        return self._specific_angular_momentum


    @property
    def semi_major_axis(self):
        if not hasattr(self,'_semi_major_axis'):
            h = self.specific_angular_momentum
            e = self.eccentricity
            μ = self.body.gravitational_parameter
            a = sqrt(a * μ * (e**2 - 1))
            self._semi_major_axis = a
        return self._semi_major_axis


    @property
    def true_anomaly_at_infinity(self):
        '''η'''
        e = self.eccentricity
        return arccos(-1/e)

    @property
    def turning_angle(self):
        '''δ'''
        e = self.eccentricity
        return 2 * arcsin(1/e)

    @property
    def impact_parameter(self):
        '''b'''
        a = self.semi_major_axis
        δ = self.turning_angle
        return -a / tan(δ/2)
