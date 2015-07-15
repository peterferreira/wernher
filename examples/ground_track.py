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

# calculate radius (to be used to set color of track
r = orbit.radius_at_time(tt)

# create figure and add map view, track and position marker
fig,ax = pyplot.subplots()
mview = wernher.MapView(orbit.body)
mview.zoomlevel = 1
mview.plot_basemap(ax)

# plot ground track
tk = wernher.MapView.plot_track(ax,lat,lon,r)

# place marker for the vessel location
mk = wernher.MapView.plot_marker(ax,
    orbit.latitude_at_epoch/deg,
    orbit.longitude_at_epoch/deg,
    color='red')

# next periapsis marker
mk = wernher.MapView.plot_marker(ax,
    orbit.latitude_at_periapsis()/deg,
    orbit.longitude_at_periapsis()/deg,
    marker='v', color='cyan')

# next apoapsis marker
mk = wernher.MapView.plot_marker(ax,
    orbit.latitude_at_apoapsis()/deg,
    orbit.longitude_at_apoapsis()/deg,
    marker='^', color='magenta')

# show plot in new window
pyplot.show()
