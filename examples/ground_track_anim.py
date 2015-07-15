import wernher
from matplotlib import pyplot, animation
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


fig,ax = pyplot.subplots()

def initfig():
    # create figure and add map view, track and position marker
    mview = wernher.MapView(orbit.body)
    mview.zoomlevel = 1
    im = mview.plot_basemap(ax)
    return im,

def updatefig(tick):
    if tick == 0:
        return []

    t = orbit.epoch + tick*60

    # ground track consists of 200 points
    npoints = 200

    # start in the past by 1/4 of the orbital period
    tmin = t - 0.25*orbit.period

    # plot 1.5 periods of ground track
    tmax = tmin + 1.5*orbit.period

    # array of times - evenly spaced
    tt = np.linspace(tmin,tmax,npoints)

    # array of lattitude and longitudes, converted to degrees
    lat = orbit.latitude_at_time(tt) / deg
    lon = orbit.longitude_at_time(tt) / deg

    # calculate radius (to be used to set color of track
    r = orbit.radius_at_time(tt)

    # plot ground track
    ret = wernher.MapView.plot_track(ax,lat,lon,r)

    # place marker for the vessel location
    ret.append(wernher.MapView.plot_marker(ax,
        orbit.latitude_at_time(t)/deg,
        orbit.longitude_at_time(t)/deg,
        color='red'))

    # next periapsis marker
    ret.append(wernher.MapView.plot_marker(ax,
        orbit.latitude_at_periapsis_after_time(t)/deg,
        orbit.longitude_at_periapsis_after_time(t)/deg,
        marker='v', color='cyan'))

    # next apoapsis marker
    ret.append(wernher.MapView.plot_marker(ax,
        orbit.latitude_at_apoapsis_after_time(t)/deg,
        orbit.longitude_at_apoapsis_after_time(t)/deg,
        marker='^', color='magenta'))

    return ret



ani = animation.FuncAnimation(
    fig,
    updatefig,
    interval=10,
    blit=True,
    init_func=initfig)

# show plot in new window
pyplot.show()
