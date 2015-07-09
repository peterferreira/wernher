import numpy as np

π = np.pi
sqrt = np.sqrt

class CelestialBody(object):
    _ra = {
        'Kerbin' :    0        ,
        'Minmus' :  140 * π/180,
        'Dres'   :  -65 * π/180,
    }

    def __init__(self,**kwargs):
        for k,v in kwargs.items():
            setattr(self,k,v)

    #self.name                    = body.name.lower()
    #self.gravitational_parameter = body.gravitational_parameter
    #self.equatorial_radius       = body.equatorial_radius
    #self.rotational_speed        = body.rotational_speed

    @poperty
    def epoch(self):
        return 0

    @property
    def right_ascension_at_epoch(self):
        # RA of prime meridian (lon = 0) at epoch (t = 0)
        return CelestialBody._ra.get(self.name,0)

    @property
    def equatorial_surface_rotational_speed(self):
        ω = self.rotational_speed
        r = self.equatorial_radius
        return ω * r

    def right_ascension_at_time(self,t):
        '''RA of prime meridian'''
        α0 = self.right_ascension_at_epoch
        t0 = self.epoch
        ω = self.rotational_speed
        return (α0 + ω * (t - t0)) % (2*π)

    @property
    def stationary_radius(self):
        μ = self.gravitational_parameter
        ω = self.rotational_speed
        return (μ / ω**2)**(1/3)

    @property
    def stationary_altitude(self):
        R = self.equatorial_radius
        r = self.stationary_radius
        return r - R

    @property
    def stationary_speed(self):
        μ = self.gravitational_parameter
        r = self.stationary_radius
        return sqrt(μ / r)

    def escape_speed(self,r):
        μ = self.gravitational_parameter
        return sqrt(2*μ / r)
