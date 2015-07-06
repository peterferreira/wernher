import numpy as np
from scipy import optimize as opt

π = np.pi
inf = np.inf
nan = np.nan

sqrt = lambda x: np.sqrt(np.float(x))
sin = np.sin
cos = np.cos
tan = np.tan
arcsin = np.arcsin
arccos = np.arccos
arctan = np.arctan
arctan2 = np.arctan2
sinh = np.sinh
cosh = np.cosh
arctanh = np.arctanh

from .celestial_body import CelestialBody
from .locked_property import cached_property, locked_property, \
                             property_alias, LockError
from .orbit_type import OrbitType, circular, elliptic, hyperbolic, \
                        parabolic
from .stumpff import stumpff_c, stumpff_s

class Orbit(object):
    '''
    Fundamental Elements
                orbit_type
        e       eccentricity
        i       inclination
        Ω       longitude_of_ascending_node
        ω       argument_of_periapsis
        a       semi_major_axis
        h       specific_angular_momentum
        M0      mean_anomaly_at_epoch
        θ0      true_anomaly_at_epoch
        χ0      universal_anomaly_at_epoch

    State Vectors
        x0      position_at_epoch
        v0      velocity_at_epoch

    Independent Parameters
        t0  epoch

    Intermediate Parameters
        r0      radius_at_epoch
        b       semi_minor_axis
        ap      apoapsis
        pe      periapsis
        ϖ       longitude_of_periapsis
        n       mean_motion
        T       period
        ap_alt  apoapsis_altitude
        pe_alt  periapsis_altitude
        vap     speed_at_apoapsis
        vpe     speed_at_periapsis

        P,Q,W   transform

    Celestial Body
        R       equatorial_radius
        μ       gravitational_parameter

    '''
    # Fundamental Elements
    e      = property_alias('eccentricity')
    i      = property_alias('inclination')
    Ω      = property_alias('longitude_of_ascending_node')
    ω      = property_alias('argument_of_periapsis')
    a      = property_alias('semi_major_axis')
    h      = property_alias('specific_angular_momentum')
    M0     = property_alias('mean_anomaly_at_epoch')
    θ0     = property_alias('true_anomaly_at_epoch')
    χ0     = property_alias('univeral_anomaly_at_epoch')

    # State Vectors
    x0     = property_alias('position_at_epoch')
    v0     = property_alias('velocity_at_epoch')

    # Reference Frame Parameters
    t0     = property_alias('epoch')

    # Intermediate Parameters
    r0     = property_alias('radius_at_epoch')
    b      = property_alias('semi_minor_axis')
    ap     = property_alias('apoapsis')
    pe     = property_alias('periapsis')
    ϖ      = property_alias('longitude_of_periapsis')
    n      = property_alias('mean_motion')
    T      = property_alias('period')
    ap_alt = property_alias('apoapsis_altitude')
    pe_alt = property_alias('periapsis_altitude')
    vap    = property_alias('speed_at_apoapsis')
    vpe    = property_alias('speed_at_periapsis')

    def __init__(self,**kwargs):
        if 'body' in kwargs:
            self.body = CelestialBody(**kwargs.pop('body'))
        for k,v in kwargs.items():
            setattr(self,k,v)

    def _try_calc(self,fnlist):
        try:
            return fnlist.pop(0)(self)
        except LockError:
            return self._try_calc(fnlist)
        except IndexError:
            raise LockError from None

    @locked_property
    def orbit_type(self):
        e = self.eccentricity
        otype = OrbitType.from_eccentricity(e)
        return otype

    @locked_property
    def eccentricity(self):

        def _1(self):
            h = self.specific_angular_momentum
            r = self.radius_at_epoch
            θ = self.true_anomaly_at_epoch
            μ = self.body.gravitational_parameter
            e = (h**2 / (μ * r) - 1) / cos(θ)
            return e

        def _2(self):
            pe = self.periapsis
            ap = self.apoapsis
            e = (ap - pe) / (ap + pe)
            return e

        def _3(self):
            r = self.radius_at_epoch
            v = self.speed_at_epoch
            vt = self.tangent_speed_at_epoch
            μ = self.body.gravitational_parameter
            e = sqrt(1 + (r * vt**2 / μ) * ((r * v**2 / μ) - 2))
            return e

        def _4(self):
            h = self.specific_angular_momentum
            pe = self.periapsis
            μ = self.body.gravitational_parameter
            e = h**2 / (μ * pe) - 1
            return e

        def _5(self):
            r = self.radius_at_epoch
            v = self.speed_at_epoch
            θ = self.true_anomaly_at_epoch
            μ = self.body.gravitational_parameter

            # solving for e:
            # A * e^2 + B * e + C = 0
            A = μ / r
            B = (2 * μ / r - v**2) * cos(θ)
            C =  μ / r - v**2

            if (B**2 - 4*A*C) < 0:
                raise LockError

            x = -B / (2*A)
            y = sqrt(B**2 - 4*A*C) / (2*A)

            if y > 0:
                e = x + y
            else:
                e = x - y

            return e

        def _6(self):
            otype = self.orbit_type
            if otype is circular:
                e = 0
            elif otype is parabolic:
                e = 1
            else:
                raise LockError from None
            return e

        e = self._try_calc([_1,_2,_3,_4,_5,_6])
        return e

    @locked_property
    def semi_major_axis(self):
        if self.orbit_type is parabolic:
            return inf
        try:
            e = self.eccentricity
            h = self.specific_angular_momentum
            μ = self.body.gravitational_parameter
            a = h**2 / (μ * (1 - e**2))
        except LockError:
            try:
                pe = self.periapsis
                ap = self.apoapsis
                a = (pe + ap) / 2
            except LockError:
                r = self.radius_at_epoch
                v = self.speed_at_epoch
                μ = self.body.gravitational_parameter
                a = 1 / (2/r - (v**2)/μ)
        return a

    @locked_property
    def specific_angular_momentum(self):
        try:
            if self.orbit_type is parabolic:
                raise LockError from None
            e = self.eccentricity
            a = self.semi_major_axis
            μ = self.body.gravitational_parameter
            h = sqrt(abs(a) * μ * abs(1 - e**2))
        except LockError:
            try:
                pe = self.periapsis
                e = self.eccentricity
                μ = self.body.gravitational_parameter
                h = sqrt(pe * μ * (1 + e))
            except LockError:
                try:
                    pe = self.periapsis
                    vpe = self.speed_at_periapsis
                    h = pe * vpe
                except LockError:
                    try:
                        vt = self.tangent_speed_at_epoch
                        r = self.radius_at_epoch
                        h = vt * r
                    except LockError:
                        p,q = self.perifocal_position_at_epoch
                        vp,vq = self.perifocal_velocity_at_epoch
                        h = (p * vq) - (q * vp)
        return h

    @locked_property
    def universal_anomaly_at_epoch(self):
        t0 = self.epoch
        χ0 = self.universal_anomaly_at_time(t0)
        return χ0

    @locked_property
    def mean_anomaly_at_epoch(self):
        θ0 = self.true_anomaly_at_epoch
        M0 = self.mean_anomaly_at_true_anomaly(θ0)
        return M0

    @locked_property
    def true_anomaly_at_epoch(self):
        try:
            e = self.eccentricity
            h = self.specific_angular_momentum
            r = self.radius_at_epoch
            μ = self.body.gravitational_parameter
            θ = arccos((h**2 / (μ * r) - 1) / e)
        except LockError:
            p,q = self.perifocal_position_at_epoch
            r = self.radius_at_epoch
            θ = arccos(p / r)

        try:
            vr = self.radial_speed_at_epoch
            if vr < 0:
                θ = 2*π - θ
        except LockError:
            try:
                vr = self.radial_speed_at_epoch
                if vr < 0:
                    θ = 2*π - θ
            except LockError:
                try:
                    p,q = self.perifocal_position_at_epoch
                    if q < 0:
                        θ = 2*π - θ
                except LockError:
                    '''not enough info to pick out quadrant'''
                    pass

        return θ

    @locked_property
    def position_at_epoch(self):
        raise LockError from None

    @locked_property
    def velocity_at_epoch(self):
        raise LockError from None

    @locked_property
    def tangent_speed_at_epoch(self):
        try:
            v = self.speed_at_epoch
            γ = self.flight_path_angle_at_epoch
            vt = v / sqrt(tan(γ)**2 + 1)
        except LockError:
            v = self.speed_at_epoch
            vr = self.radial_speed_at_epoch
            vt = sqrt(v**2 - vr**2)
        return vt

    @locked_property
    def radial_speed_at_epoch(self):
        try:
            v = self.speed_at_epoch
            γ = self.flight_path_angle_at_epoch
            vr = v / sqrt((1/tan(γ)**2) + 1)
        except LockError:
            x,y = self.position_at_epoch
            vx,vy = self.velocity_at_epoch
            r = self.radius_at_epoch
            vr = (vx*x + vy*y) / r
        return vr


    @locked_property
    def radius_at_epoch(self):
        try:
            if self.orbit_type is parabolic:
                raise LockError from None
            e = self.eccentricity
            a = self.semi_major_axis
            θ = self.true_anomaly_at_epoch
            r = a * (1 - e**2) / (1 + e * cos(θ))
        except LockError:
            try:
                x,y = self.position_at_epoch
                r = sqrt(x**2 + y**2)
            except LockError:
                p,q = self.perifocal_position_at_epoch
                r = sqrt(p**2 + q**2)
        return r

    @locked_property
    def apoapsis(self):
        try:
            if self.orbit_type is parabolic:
                raise LockError from None
            e = self.eccentricity
            a = self.semi_major_axis
            ap = a * (1 + e)
        except LockError:
            ap_alt = self.apoapsis_altitude
            R = self.body.equatorial_radius
            ap = ap_alt + R
        return ap

    @locked_property
    def periapsis(self):
        try:
            if self.orbit_type is parabolic:
                raise LockError from None
            e = self.eccentricity
            a = self.semi_major_axis
            pe = a * (1 - e)
        except LockError:
            try:
                pe_alt = self.periapsis_altitude
                R = self.body.equatorial_radius
                pe = pe_alt + R
            except LockError:
                if self.orbit_type is parabolic:
                    vpe = self.speed_at_periapsis
                    μ = self.body.gravitational_parameter
                    pe = 2 * μ / vpe**2
                else:
                    raise LockError from None
        return pe

    @locked_property
    def apoapsis_altitude(self):
        ap = self.apoapsis
        R = self.body.equatorial_radius
        ap_alt = ap - R
        return ap_alt

    @locked_property
    def periapsis_altitude(self):
        pe = self.periapsis
        R = self.body.equatorial_radius
        pe_alt = pe - R
        return pe_alt

    @locked_property
    def period(self):
        if self.orbit_type.isopen:
            return nan
        a = self.semi_major_axis
        μ = self.body.gravitational_parameter
        T = 2*π * sqrt(a**3 / μ)
        return T

    @locked_property
    def speed_at_epoch(self):
        try:
            r = self.radius_at_epoch
            v = self.speed_at_radius(r)
        except LockError:
            try:
                vx,vy = self.velocity_at_epoch
                v = sqrt(vx**2 + vy**2)
            except LockError:
                vp,vq = self.perifocal_velocity_at_epoch
                v = sqrt(vp**2 + vq**2)
        return v

    @locked_property
    def speed_at_periapsis(self):
        if self.orbit_type is parabolic:
            raise LockError from None
        a = self.semi_major_axis
        pe = self.periapsis
        μ = self.body.gravitational_parameter
        vpe = sqrt(μ * (2/pe - 1/a))
        return vpe

    @locked_property
    def speed_at_apoapsis(self):
        if self.orbit_type.isopen:
            return nan
        a = self.semi_major_axis
        ap = self.apoapsis
        μ = self.body.gravitational_parameter
        vap = sqrt(μ * (2/ap - 1/a))
        return vap

    @locked_property
    def radius_at_average_true_anomaly(self):
        if self.orbit_type is parabolic:
            raise LockError from None
        a = self.semi_major_axis
        e = self.eccentricity
        r = a * sqrt(1 - e**2)
        return r

    @locked_property
    def true_anomaly_at_semi_minor_axis(self):
        e = self.eccentricity
        θ = arccos(-e)
        return θ

    @locked_property
    def flight_path_angle_at_epoch(self):
        e = self.eccentricity
        θ = self.true_anomaly_at_epoch
        γ = arctan2(e * sin(θ), 1 + e * cos(θ))
        return γ

    @locked_property
    def flight_path_angle_max(self):
        # first, calculate true anomaly (θ) where dγ/dθ = 0
        θ = self.true_anomaly_at_semi_minor_axis
        γ = self.flight_path_angle_at_true_anomaly(θ)
        return γ

    @locked_property
    def true_anomaly_at_infinity(self):
        if self.orbit_type is not hyperbolic:
            return nan
        e = self.eccentricity
        η = arccos(-1/e)
        return η

    @locked_property
    def speed_at_infinity(self):
        if self.orbit_type.isclosed:
            return nan
        elif self.orbit_type is parabolic:
            return 0
        a = self.semi_major_axis
        μ = self.body.gravitational_parameter
        v = sqrt(μ / -a)
        return v

    @locked_property
    def turn_angle(self):
        if self.orbit_type is not hyperbolic:
            return nan
        e = self.eccentricity
        δ = 2 * arcsin(1/e)
        return δ

    @locked_property
    def impact_parameter(self):
        if self.orbit_type is not hyperbolic:
            return nan
        a = self.semi_major_axis
        δ = self.turn_angle
        b = -a / tan(δ/2)
        return b

    @locked_property
    def aiming_radius(self):
        if self.orbit_type is not hyperbolic:
            return nan
        e = self.eccentricity
        a = self.semi_major_axis
        Δ = -a * sqrt(e**2 - 1)
        return Δ

    @locked_property
    def mean_motion(self):
        T = self.period
        return 2*π / T

    ### tranformation vectors to Cartesian coordinates
    @locked_property
    def transform(self):
        '''P⃗,Q⃗,W⃗'''
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
        return (P,Q,W)

    @cached_property
    def lagrange_coefficients_at_epoch(self):
        t0 = self.epoch
        return lagrange_coefficients_at_time(t0)

    @cached_property
    def perifocal_position_at_epoch(self):
        e = self.eccentricity
        h = self.specific_angular_momentum
        θ = self.true_anomaly_at_epoch
        μ = self.body.gravitational_parameter
        x = h**2 / (μ * (1 + e * cos(θ)))
        p = x * cos(θ)
        q = x * sin(θ)
        return p,q

    @cached_property
    def perifocal_velocity_at_epoch(self):
        e = self.eccentricity
        h = self.specific_angular_momentum
        θ = self.true_anomaly_at_epoch
        μ = self.body.gravitational_parameter
        vp = - (μ / h) * sin(θ)
        vq = (μ / h) * (e + cos(θ))
        return vp,vq

    def true_anomaly_from_periapsis(self,r):
        e = self.eccentricity
        h = self.specific_angular_momentum
        μ = self.body.gravitational_parameter
        θ = arccos((h**2 / (μ * r) - 1) / e)
        return θ

    def true_anomaly_to_periapsis(self,r):
        θ0 = self.true_anomaly_from_periapsis(r)
        θ = 2*π - θ0
        return θ

    def true_anomaly_at_radius(self,r):
        e = self.eccentricity
        h = self.specific_angular_momentum
        μ = self.body.gravitational_parameter
        θ = arccos(((h**2 / (μ * r)) - 1) / e)
        return θ

    def speed_at_radius(self,r):
        if self.orbit_type is parabolic:
            raise LockError from None
        a = self.semi_major_axis
        μ = self.body.gravitational_parameter
        v = sqrt(2*((-μ/(2*a)) + (μ/r)))
        return v

    def tangent_speed_at_time(self,t):
        h = self.specific_angular_momentum
        r = self.radius_at_time(t)
        vt = h / r
        return vt

    def radial_speed_at_time(self,t):
        e = self.eccentricity
        h = self.specific_angular_momentum
        μ = self.body.gravitational_parameter
        θ = self.true_anomaly_at_time(t)
        vr = (μ / h) * e * sin(θ)
        return vr

    def speed_at_time(self,t):
        vr = self.radial_speed_at_time(t)
        vt = self.tangent_speed_at_time(t)
        v = sqrt(vt**2 + vr**2)
        return v

    def flight_path_angle_at_true_anomaly(self,θ):
        otype = self.orbit_type
        if otype is parabolic:
            γ = θ/2
        else:
            e = self.eccentricity
            γ = arctan2(e * sin(θ),1 + e * cos(θ))
        return γ

    def radius_at_true_anomaly(self,θ):
        e = self.eccentricity
        h = self.specific_angular_momentum
        μ = self.body.gravitational_parameter
        r = (h**2 / μ) * (1 / (1 + e * cos(θ)))
        return r

    def radius_at_time(self,t):
        θ = self.true_anomaly_at_time(t)
        r = self.radius_at_true_anomaly(θ)
        return r

    def position_at_time(self,t):
        x0,y0 = self.position_at_epoch
        vx0,vy0 = self.velocity_at_epoch
        f,g,dfdχ,dgdχ = self.lagrange_coefficients_at_time(t)
        x = f * x0 + g * vx0
        y = f * y0 + g * vy0
        return np.hstack([x,y])

    def position_at_true_anomaly(self,θ):
        h = self.specific_angular_momentum
        μ = self.body.gravitational_parameter
        θ0 = self.true_anomaly_at_epoch
        r0 = self.radius_at_epoch
        r = self.radius_at_true_anomaly(θ)
        Δθ = θ - θ0

        f = 1 - (μ * r / h**2) * (1 - cos(Δθ))
        g = (r * r0 / h) * sin(Δθ)

        x0,y0 = self.position_at_epoch
        vx0,vy0 = self.velocity_at_epoch

        x = f * x0 + g * vx0
        y = f * y0 + g * vy0
        return np.hstack([x,y])

    def velocity_at_time(self,t):
        x0,y0 = self.position_at_epoch
        vx0,vy0 = self.velocity_at_epoch
        f,g,dfdχ,dgdχ = self.lagrange_coefficients_at_time(t)
        vx = dfdχ * x0 + dgdχ * vx0
        vy = dfdχ * y0 + dgdχ * vy0
        return np.hstack([vx,vy])

    def velocity_at_true_anomaly(self,θ):
        h = self.specific_angular_momentum
        μ = self.body.gravitational_parameter
        θ0 = self.true_anomaly_at_epoch
        r0 = self.radius_at_epoch
        r = self.radius_at_true_anomaly(θ)
        Δθ = θ - θ0

        dfdθ = (μ/h) * ((1-cos(Δθ))/sin(Δθ)) \
             * ((μ/h**2) * (1 - cos(Δθ)) - 1/r0 - 1/r)
        dgdθ = 1 - (μ * r0 / h**2) * (1 - cos(Δθ))

        x0,y0 = self.position_at_epoch
        vx0,vy0 = self.velocity_at_epoch

        vx = dfdθ * x0 + dgdθ * vx0
        vy = dfdθ * y0 + dgdθ * vy0
        return np.hstack([vx,vy])

    def universal_anomaly_at_time(self,t):
        t0 = self.epoch
        r0 = self.radius_at_epoch
        vr0 = self.radial_speed_at_epoch
        a = self.semi_major_axis
        μ = self.body.gravitational_parameter

        α = 1 / a

        A = 1 - α * r0
        B = r0 * vr0 / sqrt(μ)
        C = r0
        D = -sqrt(μ)
        def f(t):
            def _f(χ,α=α,t0=t0,A=A,B=B,C=C,D=D):
                return A * χ**3 * stumpff_s(α * χ**2) \
                     + B * χ**2 * stumpff_c(α * χ**2) \
                     + C * χ \
                     + D * (t - t0)
            return _f

        E = -r0 * vr0 * α / sqrt(μ)
        F = 1 - α * r0
        G = r0 * vr0 / sqrt(μ)
        H = r0
        def dfdχ(χ,α=α,E=E,F=F,G=G,H=H):
            return E * χ**3 * stumpff_s(α * χ**2) \
                 + F * χ**2 * stumpff_c(α * χ**2) \
                 + G * χ \
                 + H

        def χinit(t,t0=t0):
            return sqrt(μ) * abs(α) * (t - t0)

        if hasattr(t,'__iter__'):
            χ = np.array([opt.newton(f(_t),χinit(_t),fprime=dfdχ) \
                          for _t in t])
        else:
            χ = opt.newton(f(t),χinit(t),fprime=dfdχ)

        return χ

    def lagrange_coefficients_at_time(self,t):
        t0 = self.epoch
        r0 = self.radius_at_epoch
        x0 = self.position_at_epoch
        v0 = self.velocity_at_epoch
        a = self.semi_major_axis
        μ = self.body.gravitational_parameter
        χ = self.universal_anomaly_at_time(t)

        α = 1 / a
        z = α * χ**2

        cz = stumpff_c(z)
        sz = stumpff_s(z)

        f = 1 - (χ**2 / r0) * cz
        g = (t - t0) - (1/sqrt(μ)) * χ**3 * sz

        x = f * x0 + g * v0
        r = sqrt(sum(x**2))

        dfdχ = (sqrt(μ) / (r * r0)) * (α * χ**3 * sz - χ)
        dgdχ = 1 - (χ**2 / r) * cz

        return f,g,dfdχ,dgdχ

    def mean_anomaly_at_true_anomaly(self,θ):
        otype = self.orbit_type
        if otype.isclosed:
            e = self.eccentricity
            E = self.eccentric_anomaly_at_true_anomaly(θ)
            M = E - ((e*sqrt(1-e**2)*sin(θ)) / (1+e*cos(θ)))
        elif otype is parabolic:
            M = (1/2) * tan(θ/2) + (1/6) * tan(θ/2)**3
        elif otype is hyperbolic:
            e = self.eccentricity
            E = self.eccentric_anomaly_at_true_anomaly(θ)
            M = e * sinh(E) - E
        return M

    def mean_anomaly_at_eccentric_anomaly(self,E):
        e = self.eccentricity
        M = E - e * sin(E)
        return M


    def mean_anomaly_at_time(self,t):
        t0 = self.epoch
        M0 = self.mean_anomaly_at_epoch
        otype = self.orbit_type
        if otype.isclosed:
            n = self.mean_motion
            M = M0 + n * (t - t0)
        elif otype is parabolic:
            h = self.specific_angular_momentum
            μ = self.body.gravitational_parameter
            M = M0 + (μ**2 / h**3) * (t - t0)
        elif otype is hyperbolic:
            e = self.eccentricity
            h = self.specific_angular_momentum
            μ = self.body.gravitational_parameter
            M = M0 + (μ**2 / h**3) * (e**2 - 1)**(3/2) * (t - t0)
        return M

    def eccentric_anomaly_at_mean_anomaly(self,M):
        otype = self.orbit_type
        if otype is elliptic:
            e = self.eccentricity
            def f(E,M,e=e):
                return E - e * sin(E) - M
            def dfdE(E,M,e=e):
                return 1 - e * cos(E)
            if hasattr(M,'__iter__'):
                E = np.array([opt.newton(f,π,args=(m,),fprime=dfdE) \
                              for m in M])
            else:
                E = opt.newton(f,π,args=(M,),fprime=dfdE)
            return E
        elif otype.isopen:
            e = self.eccentricity
            def f(F,M,e=e):
                return e * sinh(F) - F - M
            def dfdF(F,M,e=e):
                return e * cosh(F) - 1
            if hasattr(M,'__iter__'):
                F = np.array([opt.newton(f,π,args=(m,),fprime=dfdF) \
                              for m in M])
            else:
                F = opt.newton(f,π,args=(M,),fprime=dfdF)
            return F

    def eccentric_anomaly_at_time(self,t):
        M = self.mean_anomaly_at_time(t)
        E = self.eccentric_anomaly_at_mean_anomaly(M)
        return E

    def eccentric_anomaly_at_true_anomaly(self,θ):
        otype = self.orbit_type
        if otype.isclosed:
            e = self.eccentricity
            E = 2 * arctan2(sqrt(1-e)*sin(θ/2),sqrt(1+e)*cos(θ/2))
        else:
            e = self.eccentricity
            E = 2 * arctanh(sqrt((e-1)/(e+1)) * tan(θ/2))
        return E

    def true_anomaly_at_time(self,t):
        otype = self.orbit_type
        if otype is circular:
            n = self.mean_motion
            θ = n * t
        elif otype is elliptic:
            e = self.eccentricity
            E = self.eccentric_anomaly_at_time(t)
            θ = π - 2 * arctan2(sqrt(1-e) * cos(E/2),
                                sqrt(1+e) * sin(E/2))
        elif otype is parabolic:
            M = self.mean_anomaly_at_time(t)
            θ = 2 * arctan(
                (3*M + sqrt((3*M)**2 + 1))**(1/3)
              - (3*M + sqrt((3*M)**2 + 1))**(-1/3) )
        elif otype is hyperbolic:
            e = self.eccentricity
            F = self.eccentric_anomaly_at_time(t)
            θ = 2 * arctan2(sqrt(e+1) * sinh(F/2),
                            sqrt(e-1) * cosh(F/2))
        return θ

    def time_at_true_anomaly(self,θ):
        otype = self.orbit_type
        if otype.isclosed:
            T = self.period
            M = self.mean_anomaly_at_true_anomaly(θ)
            t = (M / (2*π)) * T
        else:
            h = self.specific_angular_momentum
            e = self.eccentricity
            M = self.mean_anomaly_at_true_anomaly(θ)
            μ = self.body.gravitational_parameter
            t = (h**3 / (μ**2 * (e**2 - 1)**(3/2))) * M
        return t

    @staticmethod
    def from_radius_at_true_anomaly(r1,θ1,r2,θ2,body):
        '''Sets the eccentricity (e) and specific angular momentum (h)
        from two points (radius,true anomaly): (r1,θ1), (r2,θ2)
        '''
        orbit = Orbit(body=body)
        μ = orbit.body.gravitational_parameter
        e = (r2 - r1) / (r1 * cos(θ1) - r2 * cos(θ2))
        h = sqrt(r1 * μ * (1 + e * cos(θ1)))
        orbit.eccentricity = e
        orbit.specific_angular_momentum = h
        return orbit

    @staticmethod
    def from_altitude_at_true_anomaly(z1,θ1,z2,θ2,body):
        '''Sets the eccentricity (e) and specific angular momentum (h)
        from two points (altitude,true anomaly): (z1,θ1), (z2,θ2)
        '''
        R = body.equatorial_radius
        r1 = z1 + R
        r2 = z2 + R
        orbit = Orbit.from_radius_at_true_anomaly(r1,θ1,r2,θ2,body)
        return orbit
