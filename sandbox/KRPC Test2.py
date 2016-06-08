# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run -i '../Common.ipynb'
import krpc

# <codecell>

conn = krpc.connect(name='laptop2', address='192.168.1.2')

ksc = conn.space_center
vessel = ksc.active_vessel
obt = vessel.orbit
ap = vessel.auto_pilot
con = vessel.control

# <codecell>

def torque(vessel):
    surfs = vessel.parts.with_module('ModuleControlSurface')
    pitch,yaw,roll = 0,0,0
    for surf in surfs:
        modl = [x.fields for x in surf.modules if x.name == 'ModuleControlSurface'][0]
        print(surf.name)
        p = surf.position(vessel.reference_frame)
        d = surf.direction(vessel.reference_frame)
        r = surf.rotation(vessel.reference_frame)
        print('   ',p)
        print('   ',d)
        print('   ',r)
        if modl['Pitch']:
            pass
        if modl['Yaw']:
            pass
        if modl['Roll']:
            pass
torque(vessel)

# <codecell>

surf = vessel.parts.with_module('ModuleControlSurface')[0]

# <codecell>

dir(surf.modules[0])

# <codecell>

vrf = vessel.reference_frame
srfrf = vessel.surface_reference_frame
vobtrf = vessel.orbital_reference_frame
obtrf = obt.body.reference_frame
obtorf = obt.body.orbital_reference_frame
obtnrrf = obt.body.non_rotating_reference_frame

flight = lambda rf: vessel.flight(rf)

# <codecell>

compare(ksc)

# <codecell>

t = ksc.ut
o = KeplerOrbit(obt)
print(o)

f = flight(obtorf)
print(f.longitude, f.speed)

f = flight(obtnrrf)
print(f.longitude, f.speed)

f = flight(obtrf)
print(f.longitude, f.speed)

