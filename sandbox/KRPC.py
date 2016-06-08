# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import krpc
import datetime
import time
%run -i '../Common.ipynb'
#conn = krpc.connect(name='laptop0', address='192.168.1.9')
conn = krpc.connect(name='laptop0', address='goetz.homenet.org')
ksc = conn.space_center

# <codecell>

ksc.active_vessel.orbit.apoapsis

# <codecell>

%load_ext autoreload
%autoreload 2
sys.path.append('home/goetz/local/lib64/python')
%aimport ksp

# <codecell>

π = np.pi
deg = π/180

def test_close(name,x,y,rtol=1e-5,atol=1e-5):
    if not np.isclose(x,y,rtol=rtol,atol=atol):
        print(name+': {} != {}'.format(x,y))

# <codecell>

def orbit_from_krpc_orbit(ksc,obt):
    body = ksp.CelestialBody(
        name = obt.body.name,
        gravitational_parameter = obt.body.gravitational_parameter,
        equatorial_radius = obt.body.equatorial_radius,
        rotational_speed = obt.body.rotational_speed)
    return ksp.Orbit(
        t0 = ksc.ut,
        i  = obt.inclination,
        Ω  = obt.longitude_of_ascending_node,
        ω  = obt.argument_of_periapsis,
        e  = obt.eccentricity,
        a  = obt.semi_major_axis,
        M0 = obt.mean_anomaly_at_epoch,
        body = body)

# <codecell>

vessel = ksc.active_vessel
obt = vessel.orbit
obtrf = obt.body.reference_frame
f = vessel.flight(obtrf)
o = orbit_from_krpc_orbit(ksc,obt)

print('celestial body:',obt.body.name)
print('orbit type:',o.orbit_type)
test_close('    ap',o.apoapsis, obt.apoapsis, 1e-3)
test_close('    pe',o.periapsis,obt.periapsis, 1e-2)
test_close('ap alt',o.apoapsis_altitude, obt.apoapsis_altitude, 1e-3)
test_close('pe alt',o.periapsis_altitude, obt.periapsis_altitude, 1e-3)
test_close('period',o.period, obt.period, 1e-3)
test_close('     b',o.semi_minor_axis, obt.semi_minor_axis, 1e-3)
test_close('     M',o.mean_anomaly_at_epoch, obt.mean_anomaly, 1e-3, 1e-3)
test_close('     r',o.radius_at_epoch, obt.radius, 1e-3)

# KRPC's orbital speed is incorrect
#test_close(' speed',o.speed_at_epoch, obt.speed)

print('true anomaly:',o.true_anomaly_at_epoch/deg)
print('lon of pe:',(o.longitude_of_periapsis + 2*π)/deg)
#print('orbital speed:',o.speed(t))

test_close('lat',o.latitude_at_epoch/deg,f.latitude, rtol=1e-2, atol=1e-3)
test_close('lon',o.longitude_at_epoch/deg,f.longitude, rtol=1e-3)
print('Δlon:',f.longitude-o.longitude_at_epoch/deg)

for t in np.linspace(o.epoch,o.epoch+5*60*60,5):
    print('t',t)
    print('   ',o.position_at_time(t))

# <codecell>

h∞ = np.inf

# <codecell>

deg = 180/π
vessel = ksc.active_vessel
obt = vessel.orbit
obtrf = obt.body.reference_frame
obtorf = obt.body.orbital_reference_frame
obtnrrf = obt.body.non_rotating_reference_frame
t = ksc.ut
o = KeplerOrbit(obt)
f = vessel.flight(obtrf)

def longitude(kepler_orbit, t=None):
    x,y,z = o.position(t)
    return arctan2(y,x) % (2*π)

def surface_longitude(kepler_orbit, t=None):
    t = t if t is not None else self.epoch
    l = self.longitude(t)
    l0 = self.body.longitude(t)
    return ((l - l0) % (2*π)) - π/2

print('true surface lon:',f.longitude)
print('surface lon:', o.surface_longitude() * deg)

# <codecell>

deg = π/180
km=1000
vessel = ksc.active_vessel
obt = vessel.orbit
o = orbit_from_krpc_orbit(ksc,obt)
print('a:',o.semi_major_axis/km)
print('μ: {:g}'.format(o.body.gravitational_parameter))
print('period:',o.period/(60*60*6))

# <codecell>

for i in range(10):
    vessel = ksc.active_vessel
    obt = vessel.orbit
    o = orbit_from_krpc_orbit(ksc,obt)
    t = o.epoch
    tpe = o.time_to_periapsis_at_epoch
    print(o.periapsis)
    print(o.radius_at_time(t+tpe))

# <codecell>

deg = π/180
km=1000
vessel = ksc.active_vessel
obt = vessel.orbit
o = orbit_from_krpc_orbit(ksc,obt)

t = o.epoch
tpe = o.time_to_periapsis_at_epoch
npoints = 10
tmin = t + tpe - 10*60
tmax = tmin + 20*60

#tmin = t - 0.25*o.period
#tmax = t + 1.25*o.period

tt = np.linspace(tmin,tmax,npoints)
lat = o.latitude_at_time(tt) / deg
lon = o.longitude_at_time(tt) / deg
r = o.radius_at_time(tt)
rmin,rmax = r.min(),r.max()
r = (r - rmin) / (rmax - rmin)

print('time to periapsis:',tpe/(60*60*6))
print('r:',o.radius_at_time(tt)/km)

mview = ksp.MapView(o.body)

fig,ax = pyplot.subplots(figsize=(16,8))
mview.plot_basemap(ax)
_=ksp.MapView.plot_track(ax,lat,lon,r)
_=ksp.MapView.plot_marker(ax,
    o.latitude_at_epoch/deg,
    o.longitude_at_epoch/deg)
_=ksp.MapView.plot_marker(ax,
    o.latitude_at_time(t+tpe)/deg,
    o.longitude_at_time(t+tpe)/deg,
    color='yellow')
#fig.savefig('/home/goetz/dres_ground_track.png')

print(o)

# <codecell>

t0 = ksc.ut
time.sleep(60)
t1 = ksc.ut
print(t0,t1)
print(t1-t0)

# <codecell>

import wernher
from matplotlib import pyplot
import numpy as np

π = np.pi
deg = π/180
km = 1000

body = wernher.CelestialBody(
    name = 'kerbin',
    gravitational_parameter = 3.5316e12,
    equatorial_radius = 600*km,
    rotational_speed = 2*π/21600)
orbit = wernher.Orbit(
    t0 = 0,
    i  = 30*deg,
    Ω  = 0*deg,
    ω  = 15*deg,
    pe_alt  = 100*km,
    ap_alt  = 200*km,
    M0 = -45*deg,
    body = body)

# ground track consists of 200 points
npoints = 200

# start in the past by 1/4 of the orbital period
tmin = orbit.epoch - 0.25*orbit.period

# plot 1.5 periods of ground track
tmax = tmin + 1.5*orbit.period

# array of times - evenly spaced
tt = np.linspace(tmin,tmax,npoints)

# array of lattitude and longitudes, converted to degrees
lat = orbit.latitude_at_time(tt) / deg
lon = orbit.longitude_at_time(tt) / deg

# calculate radius and normalize to the range [0,1]
r = orbit.radius_at_time(tt)
rmin,rmax = r.min(),r.max()
r = (r - rmin) / (rmax - rmin)

# create figure and add map view, track and position marker
fig,ax = pyplot.subplots()
mview = wernher.MapView(orbit.body)
mview.zoomlevel = 1
mview.plot_basemap(ax)
tk = wernher.MapView.plot_track(ax,lat,lon,r)
mk = wernher.MapView.plot_marker(ax,
    orbit.latitude_at_epoch/deg,
    orbit.longitude_at_epoch/deg)

# show plot in new window
fig.savefig('kerbin_ground_track.png')
pyplot.show()

# <codecell>

0.2 * 10*6*60*60

# <codecell>

tt = []
lat = []
lon = []
data = []

vessel = ksc.active_vessel
obt = vessel.orbit
obtrf = obt.body.reference_frame

for step in range(100):
    t = ksc.ut
    o = KeplerOrbit(obt)
    f = vessel.flight(obtrf)
    
    tt.append(t)
    lat.append(f.latitude)
    lon.append(f.longitude)
    data.append(o)
    
    time.sleep(1)

# <codecell>

print([d.radius() for d in data])

# <codecell>

i = 32.223*π/180
ω = 57.836*π/180
λ1 = -142.483*π/180
Ω = 111.892*π/180
ν = 25.794*π/180
l = ω + ν
λ2 = (np.arctan(np.cos(i)*np.tan(l)) + λ1)
Δλ = λ2 - λ1
print(λ2 * 180/π, Δλ * 180/π)

