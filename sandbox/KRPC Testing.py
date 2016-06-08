# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run -i 'KRPC.ipynb'

# <codecell>

conn = krpc.connect(name='laptop', address='192.168.1.9')

ksc = conn.space_center
vessel = ksc.active_vessel
obt = vessel.orbit
ap = vessel.auto_pilot
con = vessel.control

vrf = vessel.reference_frame
srfrf = vessel.surface_reference_frame
vobtrf = vessel.orbital_reference_frame
obtrf = obt.body.reference_frame
obtorf = obt.body.orbital_reference_frame
obtnrrf = obt.body.non_rotating_reference_frame

flight = lambda rf: vessel.flight(rf)

# <codecell>

t = ksc.ut
o = KeplerOrbit(obt)

f = flight(obtorf)
print(obt.time_to_apoapsis, obt.time_to_periapsis)
print(f.longitude)
print(o.Ω * 180/π)
print(o.ν * 180/π)

# <codecell>

speed = conn.add_stream(getattr, flight(srfrf), 'speed')
altitude = conn.add_stream(getattr, flight(obtrf), 'mean_altitude')
apoapsis = conn.add_stream(getattr, obt, 'apoapsis_altitude')

# <codecell>

con.throttle = 0.6

ap.set_rotation(90, 90, roll=90)

time.sleep(1)

con.activate_next_stage()

while flight(obtrf).speed < 100.:
    time.sleep(0.1)
    
ap.set_rotation(80, 90, roll=90)

while flight(obtrf).mean_altitude < 5000.:
    time.sleep(0.1)

ap.disengage()
ap.sas = True
ap.sas_mode = ksc.SASMode.prograde

while obt.apoapsis_altitude < 80000:
    time.sleep(0.1)
    
ap.sas_mode = ksc.SASMode.stability_assist
ap.sas = False

while abs(obt.eccentricity) > 0.1:
    obt.apoapsis
    ap.set_direction(, 90, roll=90)
    

    
ap.disengage()
con.throttle = 0.

# <codecell>

ksc.SASMode.prograde

# <codecell>

speed.remove()
altitude.remove()
apoapsis.remove()

# <codecell>

def prelaunch(conn):
    ksc = conn.space_center
    vessel = ksc.active_vessel
    obtbody_rf = vessel.orbit.body.reference_frame
    
    flight = vessel.flight
    ap = vessel.auto_pilot
    cont = vessel.control
    
    vessel
    


ut = conn.add_stream(getattr, ksc, 'ut')
mean_altitude = conn.add_stream(getattr, flight(), 'mean_altitude')
#position = conn.add_stream(vessel.position, obtbody_rf)

timestamp = []
altitude = []

t0 = ut()
alt = mean_altitude()
while alt < 80000:
    t1 = ut()
    alt = mean_altitude()
    if abs(t1 - t0) > 0.001:
        timestamp.append(t1)
        altitude.append(alt)
        t0 = t1
        time.sleep(1./25.)

# <codecell>

print(ut())

# <codecell>

pyplot.plot(timestamp,altitude)

# <codecell>


print(vessel.name)
print(vessel.met)
print(vessel.mass)
print(vessel.position(vessel.orbit.body.reference_frame))

# <codecell>

def latlon(vessel):
    x,y,z = vessel.position(vessel.orbit.body.reference_frame)
    r = np.sqrt(x*x + y*y + z*z)
    lat = 90. - np.arccos(y / r) * 180. / np.pi
    lon = np.arctan2(z, x) * 180. / np.pi
    return lat,lon

# <codecell>

data = []

# <codecell>

image = pyplot.imread('/home/goetz/kerbin.jpg')
fig, ax = pyplot.subplots(figsize=(15,7))
im = ax.imshow(image)
ax.set_autoscale_on(False)

xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()

lat,lon = latlon(vessel)
xmap = ((lon + 180.) / 360.) * (xmax - xmin) + xmin
ymap = ((lat + 90.) / 180.) * (ymax - ymin) + ymin

pt = ax.plot(xmap,ymap, marker='o', color='cyan')

