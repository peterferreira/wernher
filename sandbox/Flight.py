# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run -i '../Common.ipynb'
import krpc
import wernher

# <codecell>

con = krpc.connect(name='laptop0', address='192.168.1.2')
ksc = con.space_center
vessel = ksc.active_vessel

# <codecell>

vessel.flight(vessel.orbit.body.reference_frame).mean_altitude
vessel.flight(vessel.reference_frame).roll

# <codecell>


def dot_product(x, y):
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]


def magnitude(x):
    return math.sqrt(x[0]**2 + x[1]**2 + x[2]**2)

def angle_between_vectors(x, y):
    """ Compute the angle between vector x and y """
    dp = dot_product(x, y)
    if dp == 0:
        return 0
    xm = magnitude(x)
    ym = magnitude(y)
    return math.acos(dp / (xm*ym)) * (180. / math.pi)

def vessel_pitch(vessel):
    vessel_direction = vessel.direction(vessel.surface_reference_frame)

    # Get the direction of the vessel in the horizon plane
    horizon_direction = (0, vessel_direction[1], vessel_direction[2])

    # Compute the pitch - the angle between the vessels direction and the direction in the horizon plane
    pitch = angle_between_vectors(vessel_direction, horizon_direction)
    if vessel_direction[0] < 0:
        pitch = -pitch
    return pitch

# <codecell>

cont_alt = wernher.Controller(set_point=5000,kp=1/3,t0=ksc.ut)
cont_alt.min = -15
cont_alt.max =  15
cont_alt.ziegler_nichols(ku=1/3,tu=6,control_type='no_overshoot')

cont_pitch = wernher.Controller(set_point=5,kp=1/30,t0=ksc.ut)
cont_pitch.min = -1
cont_pitch.max =  1
cont_pitch.ziegler_nichols(ku=1/25,tu=1,control_type='no_overshoot')

while True:
    t = ksc.ut
    flight = vessel.flight(vessel.orbit.body.reference_frame)
    alt = flight.mean_altitude
    pitch = vessel_pitch(vessel)
    cont_pitch.set_point = cont_alt(alt,t)
    vessel.control.pitch = cont_pitch(pitch,t)
    time.sleep(0.1)

# <codecell>

con_alt = wernher.Controller(set_point=5000)
con_alt.ziegler_nichols(ku=1/2000,tu=33)
con_alt.min = -1
con_alt.max =  1

#con_pitch = wernher.Controller(set_point=0,kp=15)
#con_pitch.min = -1
#con_pitch.max =  1

while True:
    t = ksc.ut
    flight = vessel.flight(vessel.orbit.body.reference_frame)
    alt = flight.mean_altitude
    #pitch = vessel_pitch(vessel)
    #con_pitch.set_point = con_alt(alt,t)
    vessel.control.pitch = con_alt(alt,t)
    time.sleep(0.1)

