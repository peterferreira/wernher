import numpy as np
import tables
from telemetry import Telemetry
from matplotlib import pyplot


h5file = tables.open_file('telemetry.h5','r')
table = h5file.root.telemetry.readout

t = np.linspace(0,60,60*10)
a = 20*10

fig,ax = pyplot.subplots(5,1)
for i,l in zip([0,1],['altitude','speed']):
    ax[i].plot(t[:a],[x[l] for x in table.iterrows()][:a])
    ax[i].set_title(l)

for i,l,c in zip([2,3,4],['pitch','heading','roll'],['con_yaw','con_pitch','con_roll']):
    ax[i].plot(t[:a],[x[l] for x in table.iterrows()][:a])
    axr = ax[i].twinx()
    axr.plot(t[:a],[x[c] for x in table.iterrows()][:a], color='g')
    ax[i].set_title(l+' / '+c)

fig,ax = pyplot.subplots(2)

pitch = np.array([x['pitch'] for x in table.iterrows()][:a])
dpitch = np.abs(pitch[1:] - pitch[:-1])
dt = t[:a][1:] - t[:a][:-1]
ddpitch = dpitch[1:] - dpitch[:-1]
ddt = dt[1:] - dt[:-1]

ax[0].plot(t[:a][:-1],dpitch/dt)
ax[1].plot(t[:a][:-2],ddpitch/ddt)

pyplot.show()


#t, dpitch/dt
#3.57, 55
#5.6, 389
#
#ddpitch/ddt = 165 deg / s^2
#        = 2.88 rad / s^2
