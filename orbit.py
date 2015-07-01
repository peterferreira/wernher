from enum import Enum

import numpy as np
from scipy import optimize as opt

π = np.pi

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

class OrbitType(Enum):
    circular = 0
    elliptic = 1
    parabolic = 2
    hyperbolic = 3

    @property
    def isclosed(self):
        return (self.circular or self.elliptic)

    @property
    def isopen(self):
        return (self.parabolic or self.hyperbolic)

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

class Orbit(object):
    _orbital_elements = '''
        eccentricity
        inclination
        longitude_of_ascending_node
        argument_of_periapsis
        semi_major_axis
        mean_anomaly_at_epoch
    '''.split()
    _state_vectors = '''
        position_at_epoch
        velocity_at_epoch
    '''.split()
    _intermediate_values  = '''\
        b ap pe ϖ n T
        ap_alt pe_alt
        _transform _orbit_type
    '''.split()

    def __init__(self, orbit):
        ### Orbital Elements (fundamental contants of the orbit)
        self.eccentricity                = orbit.eccentricity
        self.inclination                 = orbit.inclination
        self.longitude_of_ascending_node = orbit.longitude_of_ascending_node
        self.argument_of_periapsis       = orbit.argument_of_periapsis
        self.semi_major_axis             = orbit.semi_major_axis
        self.mean_anomaly_at_epoch       = orbit.mean_anomaly_at_epoch

        ### Reference time for mean anomaly
        self.epoch = orbit.epoch

        ### Orbiting this CelestialBody
        self.body = CelestialBody(orbit.body)

        ### tolerance for determining orbit type
        self.eccentricity_tol = 1e-5

    def clear_attributes(self,attrs):
        for a in attrs:
            if hasattr(self,a):
                delattr(self,a)

    def clear_orbital_elements(self):
        self.clear_attributes(Orbit._orbital_elements)

    def clear_state_vectors(self):
        self.clear_attributes(Orbit._state_vectors)

    def clear_intermediate_results(self):
        self.clear_attributes(Orbit._intermediate_values)

    def __setattr__(self,name,value):
        '''clear intermediate values everytime any attribute is set'''
        self.clear_intermediate_results()
        if name in Orbit._orbital_elements:
            self.clear_state_vectors()
        elif name in Orbit._state_vectors:
            self.clear_orbital_elements()
        super(self.__class__,self).__setattr__(name,value)

    @property
    def orbit_type(self):
        if not hasattr(self,'_orbit_type'):
            e = self.eccentricity
            δe = self.eccentricity_tol
            self._orbit_type = OrbitType.from_eccentricity(e,δe)
        return self._orbit_type

    ### Specific Angular Momentum (h) vs Semi-Major Axis (a)
    @property
    def specific_angular_momentum(self):
        if not hasattr(self,'h'):
            e = self.eccentricity
            a = self.semi_major_axis
            μ = self.body.gravitational_parameter
            self.h = sqrt(a * μ * abs(1 - e**2))
        return self.h

    @specific_angular_momentum.setter
    def specific_angular_momentum(self,h):
        self.clear_intermediate_results()
        if hasattr(self,'a'):
            del self.a
        self.h = h

    @property
    def semi_major_axis(self):
        if not hasattr(self,'a'):
            h = self.specific_angular_momentum
            e = self.eccentricity
            μ = self.body.gravitational_parameter
            self.a = sqrt(h * μ * abs(1 - e**2))
        return self.a

    @semi_major_axis.setter
    def semi_major_axis(self,a):
        self.clear_intermediate_results()
        if hasattr(self,'h'):
            del self.h
        self.a = a




    ### Mean Anomaly at Epoch (M0) vs True Anomaly at Epoch (θ0)
    @property
    def mean_anomaly_at_epoch(self):
        if not hasattr(self,'M0'):
            e = self.eccentricity
            θ0 = self.true_anomaly_at_epoch
            self.M0 = 2 * arctan2(sqrt(1+e)*cos(θ0/2),sqrt(1-e)*sin(θ0/2)) \
                    - ((e*sqrt(1-e**2)*sin(θ0)) / (1+e*cos(θ0)))
        return self.M0

    @mean_anomaly_at_epoch.setter
    def mean_anomaly_at_epoch(self,M0):
        self.clear_intermediate_results()
        if hasattr(self,'θ0'):
            del self.θ0
        self.M0 = M0

    @property
    def true_anomaly_at_epoch(self):
        if not hasattr(self,'θ0'):
            self.θ0 = self.true_anomaly()
        return self.θ0

    @true_anomaly_at_epoch.setter
    def true_anomaly_at_epoch(self,θ0):
        self.clear_intermediate_results()
        if hasattr(self,'M0'):
            del self.M0
        self.θ0 = θ0



    ### derived constants unique to hyperbolic trajectories
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



    ### derived constants of the orbit
    @property
    def semi_minor_axis(self):
        if not hasattr(self,'b'):
            a = self.semi_major_axis
            e = self.eccentricity
            self.b = a * sqrt(abs(1 - e**2))
        return self.b

    @property
    def apoapsis(self):
        '''ap: distance from center of body'''
        if not hasattr(self,'ap'):
            e = self.eccentricity
            a = self.semi_major_axis
            self.ap = a * (1 + e)
        return self.ap

    @property
    def periapsis(self):
        '''pe: distance from center of body'''
        if not hasattr(self,'pe'):
            e = self.eccentricity
            a = self.semi_major_axis
            self.pe = a * (1 - e)
        return self.pe

    @property
    def longitude_of_periapsis(self):
        if not hasattr(self,'ϖ'):
            Ω = self.longitude_of_ascending_node
            ω = self.argument_of_periapsis
            self.ϖ = (Ω + ω) % (2*π)
        return self.ϖ

    @property
    def mean_motion(self):
        if not hasattr(self,'n'):
            e = self.eccentricity
            a = self.semi_major_axis
            μ = self.body.gravitational_parameter
            otype = self.orbit_type
            if otype in [hyperbolic,parabolic]:
                a *= -1
            self.n = sqrt(μ / a**3)
        return self.n

    @property
    def period(self):
        if not hasattr(self,'T'):
            n = self.mean_motion
            self.T = 2*π / n
        return self.T

    @property
    def apoapsis_altitude(self):
        '''ap: distance from eq radius (surface) of body'''
        if not hasattr(self,'ap_alt'):
            ap = self.apoapsis
            R = self.body.equatorial_radius
            self.ap_alt = ap - R
        return self.ap_alt

    @property
    def periapsis_altitude(self):
        '''pe: distance from eq radius (surface) of body'''
        if not hasattr(self,'pe_alt'):
            pe = self.periapsis
            R = self.body.equatorial_radius
            self.pe_alt = pe - R
        return self.pe_alt

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

    ### tranformation vectors to Cartesian coordinates
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

    ### predictions of the orbit
    def mean_anomaly(self,t=None):
        M0 = self.mean_anomaly_at_epoch
        if t is None:
            return M0
        else:
            otype = self.orbit_type
            if otype is parabolic:
                h = self.specific_angular_momentum
                μ = self.body.gravitational_parameter
                return M0 + (μ**2 / h**3) * (t - t0)
            else:
                t0 = self.epoch
                n = self.mean_motion
                return M0 + n * (t - t0)

    def eccentric_anomaly(self,t=None):
        otype = self.orbit_type
        if otype is circular:
            pass
        elif otype is elliptic:
            e = self.eccentricity
            M = self.mean_anomaly(t)
            def fn(E,e,M):
                return E - e * sin(E) - M
            if hasattr(t,'__iter__'):
                return np.array([opt.newton(fn, π, args=(e,m)) \
                                 for m in M])
            else:
                return opt.newton(fn, π, args=(e,M))
        elif otype is hyperbolic:
            e = self.eccentricity
            M = self.mean_anomaly(t)
            def fn(F,e,M):
                return e * sinh(F) - F - M
            if hasattr(t,'__iter__'):
                return np.array([opt.brentq(fn, -2*π, 2*π, args=(e,m)) \
                                 for m in M])
            else:
                return opt.brentq(fn, -2*π, 2*π, args=(e,M))
        else:
            pass


    def true_anomaly(self,t=None):
        otype = self.orbit_type
        if otype is circular:
            pass
        elif otype is elliptic:
            e = self.eccentricity
            E = self.eccentric_anomaly(t)
            return π - 2 * arctan2(sqrt(1-e) * cos(E/2),
                                   sqrt(1+e) * sin(E/2))
        elif otype is hyperbolic:
            e = self.eccentricity
            F = self.eccentric_anomaly(t)
            return π - 2 * arctan2(sqrt(e+1) * cosh(F/2),
                                   sqrt(e-1) * sinh(F/2))
        else:
            pass

    def flight_path_angle(self,t=None):
        e = self.eccentricity
        ν = self.true_anomaly(t)
        return arctan2(e * sin(ν), 1 + e * cos(ν))

    def zenith_angle(self,t=None):
        ϕ = self.flight_path_angle(t)
        return π/2 - ϕ

    def radius(self,t=None):
        e = self.eccentricity
        a = self.semi_major_axis
        if self.orbit_type in [elliptic,circular]:
            try:
                E = self.eccentric_anomaly(t)
                return a * (1 - e * cos(E))
            except:
                ν = self.true_anomaly(t)
                return a * (1 - e**2) / (1 + e * cos(ν))
        else:
            F = self.eccentric_anomaly(t)
            return a * (1 - e * cosh(F))

    def speed(self,t=None):
        μ = self.body.gravitational_parameter
        a = self.semi_major_axis
        r = self.radius(t)
        return sqrt(μ * (2/r - 1/a))

    def mean_longitude(self,t=None):
        M = self.mean_anomaly(t)
        ϖ = self.longitude_at_periapsis
        return M + ϖ

    def true_longitude(self,t=None):
        ν = self.true_anomaly(t)
        ϖ = self.longitude_of_periapsis
        return ν + ϖ

    def universal_anomaly(self,t=None):
        e = self.eccentricity
        if e < 1:
            a = self.semi_major_axis
            E = self.eccentric_anomaly(t)
            return sqrt(a) * E
        elif e > 1:
            a = self.semi_major_axis
            F = self.eccentric_anomaly(t)
            return sqrt(-a) * F
        else:
            h = self.specific_angular_momentum
            θ = self.true_anomaly(t)
            μ = self.body.graviational_parameter
            return (h / μ) * tan(θ/2)

    def position(self,t=None):
        '''celestial position'''
        a = self.semi_major_axis
        b = self.semi_minor_axis
        e = self.eccentricity
        P,Q,W = self.transform
        if hasattr(t,'__iter__'):
            P = P.reshape((3,1))
            Q = Q.reshape((3,1))
            W = W.reshape((3,1))
        if e < 1:
            E = self.eccentric_anomaly(t)
            if hasattr(t,'__iter__'):
                E = E.reshape((1,len(E)))
            return a * (cos(E) - e) * P + b * sin(E) * Q
        else:
            F = self.eccentric_anomaly(t)
            if hasattr(t,'__iter__'):
                F = F.reshape((1,len(F)))
            return a * (cosh(F) - e) * P - b * sinh(F) * Q

    def velocity(self,t=None):
        '''celestial velocity'''
        n = self.mean_motion
        a = self.semi_major_axis
        b = self.semi_minor_axis
        e = self.eccentricity
        P,Q,W = self.transform
        if hasattr(t,'__iter__'):
            P = P.reshape((3,1))
            Q = Q.reshape((3,1))
            W = W.reshape((3,1))
        if e < 1:
            E = self.eccentric_anomaly(t)
            dEdt = n / (1 - e * cos(E))
            if hasattr(t,'__iter__'):
                E = E.reshape((1,len(E)))
            return - a * sin(E) * dEdt * P + b * cos(E) * dEdt * Q
        else:
            F = self.eccentric_anomaly(t)
            dFdt = n / (1 - e * cosh(F))
            if hasattr(t,'__iter__'):
                F = F.reshape((1,len(F)))
            return a * sinh(F) * dFdt * P - b * cosh(F) * dFdt * Q

    ### predictions with respect to the surface of the orbited body
    def latitude(self,t=None):
        '''(-π/2,π/2]'''
        x,y,z = self.position(t)
        r = self.radius(t)
        return arcsin(z/r)

    def longitude(self,t=None):
        '''celestial longitude [0,2π)'''
        t = t if t is not None else self.epoch
        x,y,z = self.position(t)
        λ0 = self.body.right_ascension(t)
        return ((arctan2(y,x) - λ0 + π/2) % (2*π)) - π

    def surface_speed(self,t=None):
        v = self.speed(t)
        vo = self.body.equatorial_rotational_speed
        #return v - vo
        return None

