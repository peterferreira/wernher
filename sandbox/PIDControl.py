# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run -i '../Common.ipynb'
from wernher import Controller

# <codecell>

π = np.pi
arctan = np.arctan

class Device(object):
    def __init__(self):
        # position, velocity, acceleration
        self.x = 0
        self.v = 0
        
        # mass, drag, spring-constant, spring equilibrium point
        self.m = 1
        self.b = 0
        self.k = 0
        self.x0 = 0
        
        # gain (maximum control), sensitivity
        self.g = 1
        self.s = 1
        
    def control_force(self,c):
        '''represents physical limitation of the control'''
        
        # bring instance variables into local scope
        g = self.g
        s = self.s
        
        return g * (2/π) * arctan((s/g) * c)
        
    def force(self,c):
        '''control force plus the external forces (drag, etc)'''
        
        # bring instance variables into local scope
        x = self.x
        v = self.v
        b = self.b
        k = self.k
        x0 = self.x0
        
        F = self.control_force(c)
        
        return F - b*v - k*(x - x0)
        
    def __call__(self,c,Δt):
        '''set the control to c for Δt time
        
        return the position and velocity after Δt has passed
        
        This device has a characteristic "mass" which gives
        it a certain momentum that the control must fight against
            
            ΣF = Δp / Δt
            ΣF = m * Δv / Δt
            v = v + (ΣF/m)*Δt
            x = x + v*Δt
        '''
        # bring instance variables into local scope
        x = self.x
        v = self.v
        m = self.m
        
        # sum of forces
        ΣF = self.force(c)
        
        # change in velocity and position
        x = x + v*Δt
        v = v + (ΣF/m)*Δt
        
        # save parameters to class instance
        self.x = x
        self.v = v
        
        return x,v

# <codecell>

def plot_behavior(xset,kp,ki,kd,m,b,k,g=1,s=1, total_time=10, npoints=100):
    
    dev = Device()
    dev.m = m
    dev.b = b
    dev.k = k
    dev.g = g
    dev.s = s
    
    cont = Controller()
    cont.kp = kp
    cont.ki = ki
    cont.kd = kd
    
    tt = np.linspace(0,total_time, npoints)
    
    xx0 = np.zeros(tt.shape)
    xx0[5:] = xset
    
    xx = [0]
    vv = [0]
    aa = [0]
    
    t0 = tt[0]
    for t,x0 in zip(tt[1:],xx0[1:]):
        cont.set_point = x0
        a = cont(dev.x,t)
        x,v = dev(a,t-t0)
        aa.append(a)
        xx.append(x)
        vv.append(v)
        t0 = t
    
    fig,ax = pyplot.subplots(1,2,figsize=(10,3))
    ptset, = ax[1].plot(tt, xx0, color='blue',label='set point')
    ptcont, = ax[1].plot(tt, aa, color='red',label='control')
    ptval, = ax[0].plot(tt, xx, color='green',label='actual')
    ptvel, = ax[1].plot(tt, vv, color='magenta',label='velocity')
    ax[0].margins(0,0.1)
    ax[1].margins(0,0.1)
    l = ax[1].legend((ptset,ptcont,ptval,ptvel), ('set point','control','actual','velocity'),
              loc='best')
    l.set_zorder(1)
    
    return tt, xx0, aa, xx, vv

# <codecell>

plot_behavior(xset=1, kp=0, ki=0, kd=0, m=1, b=10, k=20, total_time=3, npoints=100)

# <codecell>

plot_behavior(xset=1, kp=300, ki=0, kd=0, m=1, b=10, k=20, g=300, s=1, total_time=2, npoints=100)
plot_behavior(xset=1, kp=300, ki=0, kd=0, m=1, b=10, k=20, g=1000, s=1, total_time=2, npoints=100)
plot_behavior(xset=1, kp=300, ki=0, kd=0, m=1, b=10, k=20, g=1000, s=2, total_time=2, npoints=100)

# <codecell>

plot_behavior(xset=1, kp=30, ki=70, kd=0, m=1, b=10, k=20, g=30, s=3, total_time=4, npoints=200)

# <codecell>

plot_behavior(xset=1, kp=350, ki=300, kd=50, m=1, b=10, k=20, g=350, s=1, total_time=4, npoints=400)

# <codecell>

kp = 350
tt,xx0,aa,xx,vv = plot_behavior(xset=1, kp=kp, ki=0, kd=0, m=1, b=10, k=20, g=1000, s=2, total_time=5, npoints=600)

Fs = 1/(tt[1] - tt[0]) # sampling rate
n = len(xx[len(xx)//2:]) # length of the signal
k = np.arange(n)
T = n/Fs
frq = k/T # two sides frequency range
frq = frq[range(n//2)] # one side frequency range

Y = np.fft.fft(xx[len(xx)//2:])/n # fft computing and normalization
Y = Y[range(n//2)]

fig, ax = plt.subplots()
ax.plot(frq,abs(Y),'r') # plotting the spectrum
ax.set_xlabel('Freq (Hz)')
ax.set_ylabel('|Y(freq)|')
ax.set_ylim(0,1.1*max(abs(Y[2:])))

cont = Controller(set_point=1)
cont.ziegler_nichols(ku=kp,tu=1/frq[2:][np.argmax(Y[2:])],control_type='no_overshoot')

print(cont.kp,cont.ki,cont.kd)

_=plot_behavior(xset=1, kp=cont.kp, ki=cont.ki, kd=cont.kd, m=1, b=10, k=20, g=1000, s=2, total_time=5, npoints=600)

# <codecell>

g,s = 300,1
x = np.linspace(-5*g,5*g,100)
fig,ax = pyplot.subplots(figsize=(10,10))
pt = ax.plot(x,g * (2/π) * np.arctan((s/g) * x))
ax.set_aspect('equal')
ax.grid(True)

