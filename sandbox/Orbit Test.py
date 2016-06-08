# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run -i '../Common.ipynb'

# <codecell>

%load_ext autoreload
%autoreload 2
sys.path.append('home/goetz/local/lib64/python')
%aimport ksp

# <codecell>

class Bunch(dict):
    def __contains__(self, k):
        try:
            return dict.__contains__(self, k) or hasattr(self, k)
        except:
            return False
    def __getattr__(self, k):
        try:
            return object.__getattribute__(self, k)
        except AttributeError:
            try:
                return self[k]
            except KeyError:
                raise AttributeError(k)
    def __setattr__(self, k, v):
        try:
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                self[k] = v
            except:
                raise AttributeError(k)
        else:
            object.__setattr__(self, k, v)
    def __delattr__(self, k):
        try:
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                del self[k]
            except KeyError:
                raise AttributeError(k)
        else:
            object.__delattr__(self, k)

# <codecell>

π = np.pi
deg = 180/π

b = Bunch(
    name = 'Kerbin',
    gravitational_parameter = 3.5316e12,
    equatorial_radius = 600000,
    rotational_speed = 2*π / 21600)

oe = Bunch(
    eccentricity = 0.5,
    inclination = π/3,
    longitude_of_ascending_node = π+0.1,
    argument_of_periapsis = π-2,
    semi_major_axis = 3420030,
    mean_anomaly_at_epoch = π+0.1,
    epoch = 1000,
    body = b)

o = ksp.Orbit(oe)
eo = ksp.EllipticOrbit(oe)

assert o.true_anomaly_at_epoch == eo.true_anomaly_at_epoch
assert o.specific_angular_momentum == eo.specific_angular_momentum

oe.eccentricity = 1.4

ho = ksp.HyperbolicOrbit(oe)

assert o.true_anomaly_at_epoch == ho.true_anomaly_at_epoch
assert o.specific_angular_momentum == ho.specific_angular_momentum

