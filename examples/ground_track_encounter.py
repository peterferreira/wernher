import wernher
from matplotlib import pyplot, cm
import numpy as np


'''
This is a hyperbolic orbit and we are only interested in
the ground track around periapsis. First, we get the
speed at periapsis and the radius of the celestial body
and use these to choose an appropriate time interval
to plot.
'''

π = np.pi
deg = π/180
km = 1000

orbit = wernher.Orbit(
    i  = 0.14149227768205455,
    Ω  = 3.977254620789031,
    ω  = 0.22395653996553322,
    e  = 1.0779572208620696,
    a  = -8855744.039847286,
    M0 = -6.0078569863130475,
    t0 = 34056398.642449796,
    body = wernher.CelestialBody(
        name = 'Kerbin',
        equatorial_radius = 600000.0,
        gravitational_parameter = 3531600035840.0,
        rotational_speed = 0.0002908894093707204,
    ),
)

R = orbit.body.equatorial_radius
vpe = orbit.speed_at_periapsis
Δt = 10*π*R / vpe

t0 = orbit.epoch
tpe = orbit.time_to_periapsis_at_epoch

# ground track consists of 200 points
npoints = 200

# start in the past by 1/4 of the orbital period
tmin = t0 + tpe - 0.5*Δt

# plot 1.5 periods of ground track
tmax = tmin + Δt

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
tk = wernher.MapView.plot_track(ax,lat,lon,r)#,cmap=cm.cubehelix_r)

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

# show plot in new window
pyplot.show()
