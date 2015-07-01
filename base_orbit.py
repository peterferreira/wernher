import numpy as np
from scipy import optimize as opt

π = np.pi
inf = np.inf
nan = np.nan

sqrt = np.sqrt
sin = np.sin
cos = np.cos
tan = np.tan
arcsin = np.arcsin
arccos = np.arccos
arctan = np.arctan
arctan2 = np.arctan2
sinh = np.sinh
cosh = np.cosh

from .celestial_body import CelestialBody
from .orbit_type import OrbitType

class BaseOrbit(object):
    _orbital_elements = '''
        _eccentricity
        inclination
        longitude_of_ascending_node
        argument_of_periapsis
        _semi_major_axis
        _mean_anomaly_at_epoch
        _specific_angular_momentum
        _true_anomaly_at_epoch
    '''.split()
    _state_vectors = '''
        position_at_epoch
        velocity_at_epoch
    '''.split()
    _intermediate_values  = '''
        _semi_minor_axis
        _apoapsis
        _periapsis
        _longitude_of_periapsis
        _mean_motion
        _period
        _apoapsis_altitude
        _periapsis_altitude
        _transform
        _orbit_type
    '''.split()

    def __init__(self,**kwargs):
        for k,v in kwargs.items():
            setattr(self,k,v)

    '''
    def clear_attributes(self,attrs):
        for a in attrs:
            if hasattr(self,a):
                delattr(self,a)

    def clear_orbital_elements(self):
        self.clear_attributes(BaseOrbit._orbital_elements)

    def clear_state_vectors(self):
        self.clear_attributes(BaseOrbit._state_vectors)

    def clear_intermediate_results(self):
        self.clear_attributes(BaseOrbit._intermediate_values)

    def __setattr__(self,name,value):
        self.clear_intermediate_results()
        if name in BaseOrbit._orbital_elements:
            self.clear_state_vectors()
        elif name in BaseOrbit._state_vectors:
            self.clear_orbital_elements()
        super().__setattr__(name,value)
    '''
    @property
    def apoapsis(self):
        if not hasattr(self,'_apoapsis'):
            try:
                e = self.eccentricity
                a = self.semi_major_axis
                ap = a * (1 + e)
            except AttributeError:
                ap_alt = self.apoapsis_altitude
                R = self.body.equatorial_radius
                ap = ap_alt - R
            self._apoapsis = ap
        return self._apoapsis

    @apoapsis.setter
    def apoapsis(self,ap):
        self._apoapsis = ap
        if hasattr(self,'eccentricity'):
            if hasattr(self,'semi_major_axis'):
                del self.semi_major_axis
        if hasattr(self.body,'equatorial_radius'):
            if hasattr(self,'apoapsis_altitude'):
                del self.apoapsis_altitude

    @property
    def periapsis(self):
        if not hasattr(self,'_periapsis'):
            try:
                e = self.eccentricity
                a = self.semi_major_axis
                pe = a * (1 - e)
            except AttributeError:
                pe_alt = self.periapsis_altitude
                R = self.body.equatorial_radius
                pe = pe_alt - R
            self._periapsis = pe
        return self._periapsis

    @periapsis.setter
    def periapsis(self,pe):
        self._periapsis = pe
        if hasattr(self,'eccentricity'):
            if hasattr(self,'semi_major_axis'):
                del self.semi_major_axis
        if hasattr(self.body,'equatorial_radius'):
            if hasattr(self,'periapsis_altitude'):
                del self.periapsis_altitude

    @property
    def eccentricity(self):
        if not hasattr(self,'_eccentricity'):
            pe = self.periapsis
            ap = self.apoapsis
            self._eccentricity = (ap - pe) / (ap + pe)
        return self._eccentricity

    @eccentricity.setter
    def eccentricity(self,e):
        self._eccentricity = e

    @property
    def transform(self):
        '''P⃗,Q⃗,W⃗'''
        if not hasattr(self,'_transform'):
            i = self.inclination
            Ω = self.longitude_of_ascending_node
            ω = self.argument_of_periapsis

            si = sin(i)
            ci = cos(i)
            sΩ = sin(Ω)
            cΩ = cos(Ω)
            sω = sin(ω)
            cω = cos(ω)

            P = (cΩ * cω - sΩ * ci * sω,
                 sΩ * cω + cΩ * ci * sω,
                 si * sω)
            Q = (-cΩ * sω - sΩ * ci * cω,
                 -sΩ * sω + cΩ * ci * cω,
                 si * cω)
            W = (si*cω, si*sΩ, -si*sΩ)

            P = np.array(P)
            Q = np.array(Q)
            W = np.array(W)
            self._transform = (P,Q,W)
        return self._transform

'''
    @property
    def eccentricity(self):
        return _eccentricity
    @eccentricity.setter
    def eccentricity(self,e):
        assert e>0, 'Eccentricity must be positive. Value given: '+str(e)
        self._eccentricity = e

    @property
    def longitude_of_ascending_node(self):
        return _longitude_of_ascending_node
    @property
    def argument_of_periapsis(self):
        return _argument_of_periapsis
    @property
    def semi_major_axis(self):
        return _semi_major_axis
    @property
    def specific_angular_momentum(self):
        return _specific_angular_momentum
    @property
    def mean_anomaly_at_epoch(self):
        return _mean_anomaly_at_epoch
    @property
    def true_anomaly_at_epoch(self):
        return _true_anomaly_at_epoch
    @property
    def position_at_epoch(self):
        return _position_at_epoch
    @property
    def velocity_at_epoch(self):
        return _velocity_at_epoch
    @property
    def semi_minor_axis(self):
        return _semi_minor_axis
    @property
    def apoapsis(self):
        return _apoapsis
    @property
    def periapsis(self):
        return _periapsis
    @property
    def longitude_of_periapsis(self):
        return _longitude_of_periapsis
    @property
    def mean_motion(self):
        return _mean_motion
    @property
    def period(self):
        return _period
    @property
    def apoapsis_altitude(self):
        return _apoapsis_altitude
    @property
    def periapsis_altitude(self):
        return _periapsis_altitude
    @property
    def transform(self):
        return _transform

    @property
    def orbit_type(self):
        if not hasattr(self,'_orbit_type'):
            e = self.eccentricity
            δe = self.eccentricity_tol
            self._orbit_type = OrbitType.from_eccentricity(e,δe)
        return self._orbit_type

    @property
    def speed_at_periapsis(self):
        a = self.semi_major_axis
        pe = self.periapsis
        μ = self.body.gravitational_parameter
        return sqrt(μ * (2/pe - 1/a))

    @property
    def escape_speed_at_periapsis(self):
        pe = self.periapsis
        μ = self.body.gravitational_parameter
        return sqrt(2*μ / pe)

    @property
    def excess_speed_at_periapsis(self):
        vpe = self.speed_at_periapsis
        vesc = self.escape_speed_at_periapsis
        return vpe - vesc
'''
