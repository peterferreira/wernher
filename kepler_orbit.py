import os
import numpy as np
from scipy import optimize as opt
from matplotlib import pyplot

from mpl_toolkits.basemap import Basemap

π = np.pi

sqrt = np.sqrt
sin = np.sin
cos = np.cos
tan = np.tan
arcsin = np.arcsin
arccos = np.arccos
arctan = np.arctan
arctan2 = np.arctan2

def plot_basemap(bodyname,bodyradius,*args,**kwargs):
    maptype = kwargs.pop('maptype','sat')
    zoomlevel = kwargs.pop('zoomlevel',0)

    curdir = os.path.dirname(os.path.realpath(__file__))
    fpath = os.path.join(curdir,
        'map_images',bodyname.lower(),maptype,str(zoomlevel))
    ffmt = '{col}_{row}.png'

    images = []

    def cols(zoom):
        return range(2**(zoom+1))

    def rows(zoom):
        return range(2**zoom)

    for col in cols(zoomlevel):
        imrow = []
        for row in rows(zoomlevel):
            fname = os.path.join(fpath,ffmt.format(col=col,row=row))
            data = pyplot.imread(fname)
            data_u8 = (data[::-1,...,...] * np.iinfo(np.uint8).max)\
                        .astype(np.uint8,casting='unsafe')
            imrow.append(data_u8)
            del data
        images.append(np.vstack(imrow))

    image = np.hstack(images)

    fig, ax = pyplot.subplots(figsize=(15,7))

    use_basemap = False
    if use_basemap:
        bmap = Basemap(resolution=None,rsphere=bodyradius,ax=ax)
        im = bmap.imshow(image)
        major_parallels = np.linspace(-90,90,2*3+1)
        major_meridians = np.linspace(-180,180,2*4+1)
        line_opts = dict(labels=[1,0,0,1], dashes=[3,9])
        bmap.drawparallels(major_parallels, **line_opts)
        bmap.drawmeridians(major_meridians, **line_opts)
    else:
        im = ax.imshow(image,extent=[-180,180,-90,90],aspect='auto')

if __name__ == '__main__':
    plot_basemap('duna',320000,maptype='biome',zoomlevel=2)
    pyplot.show()

class CelestialBody(object):
    _ra = {
        'Minmus' :  140 * π/180,
        'Dres'   :  -65 * π/180,
    }

    def __init__(self,body):
        self.name = body.name
        if self.name not in ['Sun','Kerbol']:
            self.orbit = Orbit(body.orbit)
        else:
            self.orbit = None

        self.gravitational_parameter = body.gravitational_parameter
        self.equatorial_radius = body.equatorial_radius
        self.rotational_speed = body.rotational_speed

        self.epoch = 0

    @property
    def right_ascension_at_epoch(self):
        # RA of prime meridian (lon = 0) at epoch (t = 0)
        return CelestialBody._ra.get(self.name,0)

    @property
    def equatorial_surface_rotational_speed(self):
        ω = self.rotational_speed
        r = self.equatorial_radius
        return ω * r

    def right_ascension(self,t):
        '''RA of prime meridian'''
        α0 = self.right_ascension_at_epoch
        t0 = self.epoch
        ω = self.rotational_speed
        return (α0 + ω * (t - t0)) % (2*π)

class Orbit(object):

    def __init__(self, orbit):

        ### Orbital Elements (fundamental contants of the orbit)
        self.eccentricity                = orbit.eccentricity                # e
        self.semi_major_axis             = orbit.semi_major_axis             # a
        self.inclination                 = orbit.inclination                 # i
        self.longitude_of_ascending_node = orbit.longitude_of_ascending_node # Ω
        self.argument_of_periapsis       = orbit.argument_of_periapsis       # ω
        self.mean_anomaly_at_epoch       = orbit.mean_anomaly_at_epoch       # M0

        ### Reference time for mean anomaly
        self.epoch                       = orbit.epoch                       # t0

        ### Orbiting this CelestialBody
        self.body = CelestialBody(orbit.body)

    @staticmethod
    def from_cartesian(position,velocity):
        '''Orbit from r⃗, v⃗'''
        pass

    ### derived constants of the orbit
    @property
    def semi_minor_axis(self):
        '''b'''
        a = self.semi_major_axis
        e = self.eccentricity
        return a * sqrt(1 - e**2)

    @property
    def apoapsis(self):
        '''ap: distance from center of body'''
        e = self.eccentricity
        a = self.semi_major_axis
        return a * (1 + e)

    @property
    def periapsis(self):
        '''pe: distance from center of body'''
        e = self.eccentricity
        a = self.semi_major_axis
        return a * (1 - e)

    @property
    def longitude_of_periapsis(self):
        '''ϖ'''
        Ω = self.longitude_of_ascending_node
        ω = self.argument_of_periapsis
        return (Ω + ω) % (2*π)

    @property
    def mean_motion(self):
        '''n'''
        μ = self.body.gravitational_parameter
        a = self.semi_major_axis
        return sqrt(μ / a**3)

    @property
    def period(self):
        '''T'''
        n = self.mean_motion
        return 2*π / n

    @property
    def apoapsis_altitude(self):
        '''ap: distance from eq radius (surface) of body'''
        ap = self.apoapsis
        R = self.body.equatorial_radius
        return ap - R

    @property
    def periapsis_altitude(self):
        '''pe: distance from eq radius (surface) of body'''
        pe = self.periapsis
        R = self.body.equatorial_radius
        return pe - R

    ### tranformation vectors to Cartesian coordinates
    @property
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
        return P,Q,W

    ### predictions of the orbit
    def mean_anomaly(self,t=None):
        M0 = self.mean_anomaly_at_epoch
        if t is None:
            return M0
        else:
            t0 = self.epoch
            n = self.mean_motion
            return M0 + n * (t - t0)

    def true_anomaly(self,t=None):
        e = self.eccentricity
        E = self.eccentric_anomaly(t)
        return π - 2 * arctan2(sqrt(1-e) * cos(0.5*E),
                               sqrt(1+e) * sin(0.5*E))

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
        ν = self.true_anomaly(t)
        return a * (1 - e**2) / (1 + e * cos(ν))

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

    def eccentric_anomaly(self,t=None):
        '''
        This only works for eccentricity < (1 - ε) where ε is small
        Above this, one should use the universal variable formulation
        '''
        e = self.eccentricity
        M = self.mean_anomaly(t)
        def fn(E,e,M):
            return E - e * sin(E) - M
        if hasattr(t,'__iter__'):
            return np.array([opt.newton(fn, 0.5, args=(e,m)) for m in M])
        else:
            return opt.newton(fn, 0.5, args=(e,M))

    _='''
    @property
    def s(self,t=None):
        '' '
        universal formulation variable which satisfies:
        ds / dt = 1 / r
        '' '

        def _stumpff_c1(x):
            sqrtx = sqrt(x)
            return sin(sqrtx) / sqrtx
        def _stumpff_c2(x):
            return

        t = t if t is not None else self.epoch
        a = self.semi_major_axis
        μ = self.body.gravitational_parameter
        r = self.radius()
        drdt =


        def fn(s,a,μ,r0,drdt,t,t0):
            return r0 * s * stumpff_c1(α * s**2) \
                 + r0 * dr0dt * s**2 * stumpff_c2(α * s**2) \
                 + μ * s**3 * stumpff_c3(α * s**2) \
                 + t0 - t
    '''


    def position(self,t=None):
        '''celestial position'''
        a = self.semi_major_axis
        b = self.semi_minor_axis
        e = self.eccentricity
        E = self.eccentric_anomaly(t)
        P,Q,W = self.transform
        if hasattr(t,'__iter__'):
            E = E.reshape((1,len(E)))
            P = P.reshape((3,1))
            Q = Q.reshape((3,1))
            W = W.reshape((3,1))
        return a * (cos(E) - e) * P + b * sin(E) * Q

    def velocity(self,t=None):
        '''celestial velocity'''
        n = self.mean_motion
        a = self.semi_major_axis
        b = self.semi_minor_axis
        e = self.eccentricity
        E = self.eccentric_anomaly(t)
        P,Q,W = self.transform
        dEdt = n / (1 - e * cos(E))
        if hasattr(t,'__iter__'):
            E = E.reshape((1,len(E)))
            P = P.reshape((3,1))
            Q = Q.reshape((3,1))
            W = W.reshape((3,1))
        return - a * sin(E) * dEdt * P + b * cos(E) * dEdt * Q

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
