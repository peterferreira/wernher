import numpy as np
import tables
from telemetry import Telemetry
from matplotlib import pyplot


h5file = tables.open_file('science1b_telemetry.h5','r')
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

pyplot.show()
