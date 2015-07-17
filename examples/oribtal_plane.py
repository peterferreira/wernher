import numpy as np
from matplotlib import pyplot
import wernher

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
    pe_alt  = 800*km,
    ap_alt  = 2300*km,
    M0 = -45*deg,
    body = body)


θθ = np.linspace(0,2*π,300)
rr = orbit.radius_at_true_anomaly(θθ)
R = orbit.body.equatorial_radius

fig,ax = pyplot.subplots()
ax.set_aspect('equal')

# planet
ax.fill_between(R*np.cos(θθ), R*np.sin(θθ), lw=0)

# orbit
ax.plot(rr*np.cos(θθ), rr*np.sin(θθ), lw=2)

ax.annotate(str(orbit.pe/km)+' km',
    xy=( orbit.pe,0,),
    xycoords = 'data',
    xytext=(-10,0),
    textcoords = 'offset points',
    ha="right", va="center",
    bbox=dict(boxstyle="round", color='white', linewidth=0, alpha=0.7),
    )
ax.annotate(str(orbit.ap/km)+' km',
    xy=(-orbit.ap,0),
    xycoords = 'data',
    xytext=(10,0),
    textcoords = 'offset points',
    ha="left", va="center",
    bbox=dict(boxstyle="round", color='white', linewidth=0, alpha=0.7),
    )

pyplot.show()
