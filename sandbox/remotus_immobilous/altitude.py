import time
import math
import krpc

linkup = krpc.connect('192.168.1.2', name='altitude readout')
ksc = linkup.space_center
ves = ksc.active_vessel

alt = linkup.add_stream(getattr, ves.flight(), 'mean_altitude')

def round_distance(dist):
    fmt = '{dist:.0f} {pre}m'
    prefix = ['','k','M','G']
    d = int(math.log(dist+0.1, 10) / 3)
    i = d % len(prefix)
    return fmt.format(dist=dist/(10**(3*d)), pre=prefix[i])


while True:
    print(round_distance(alt()))
    time.sleep(1)
