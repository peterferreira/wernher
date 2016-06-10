import numpy as np

inf = np.inf
isclose = np.isclose

class Controller(object):
    '''Single Axis PID Controller'''
    def __init__(self,set_point=0,limits=(-inf,inf),kp=1,ki=0,kd=0,t0=0):
        # set point and control constants
        self.set_point = set_point
        self.cmin,self.cmax = limits
        self.kp = kp
        self.ki = ki
        self.kd = kd

        # time of previous call, running integral and
        # proportional term of previous call
        self.t0 = t0
        self.kiI  = 0
        self.P0 = 0

        # response value of previous call
        self.c = 0

    def __call__(self,x,t):
        # return previous value if no time has passed
        if isclose(t - self.t0, 0):
            return self.c

        # bring instance variables into local scope
        xset = self.set_point
        cmin = self.cmin
        cmax = self.cmax
        kp = self.kp
        ki = self.ki
        kd = self.kd

        # if parameters are all zero or None, return set point
        if not any([kp,ki,kd]):
            self.t0 = t
            return xset

        # bring instance variables into local scope
        t0 = self.t0
        kiI  = self.kiI
        P0 = self.P0

        # calculate PID terms
        Δt = t - t0
        P = xset - x
        ΔP = P - P0
        D = ΔP / Δt

        if self.ki > 0:
            if not (cmin < kp*P < cmax):
                kiI = 0
            else:
                kiI += ki * P * Δt
                kiI = min(max(kiI,cmin),cmax)

        c = kp*P + kiI + kd*D

        # clip output to specified limits
        c = min(max(c,cmin),cmax)

        # save parameters to class instance
        self.t0 = t
        self.kiI = kiI
        self.P0 = P
        self.c  = c

        return c

    @property
    def ti(self):
        '''integral time'''
        return self.kp / self.ki
    @ti.setter
    def ti(self,ti):
        self.ki = self.kp / ti

    @property
    def td(self):
        '''derivative time'''
        return self.kd / self.kp
    @td.setter
    def td(self,td):
        self.kd = self.kp * td

    @property
    def ku(self):
        '''ultimate gain, assuming classic ziegler-nichols pid scheme'''
        return (1/.6)*self.kp
    @ku.setter
    def ku(self,ku):
        self.kp = .6*ku

    @property
    def tu(self):
        '''period of oscillation at ultimate gain'''
        return 2*self.kp/self.ki
    @tu.setter
    def tu(self,tu):
        self.ki = 2*self.kp/tu
        self.kd = self.kp*tu/8

    def ziegler_nichols(self,ku,tu,control_type='pid'):
        '''
            ku = ultimate gain
            tu = period of oscillation at ultimate gain
        '''
        converter = dict(
            p = lambda ku,tu: (.5*ku, 0, 0),
            pi = lambda ku,tu: (.45*ku, 1.2*(.45*ku)/tu, 0),
            pd = lambda ku,tu: (.8*ku, 0, (.8*ku)*tu/8),
            pid = lambda ku,tu: (.6*ku, 2*(.6*ku)/tu, (.6*ku)*tu/8),
            pessen = lambda ku,tu: (.7*ku, 2.5*(.7*ku)/tu, 3*(.7*ku)*tu/20),
            some_overshoot = lambda ku,tu: (.33*ku, 2*(.33*ku)/tu, (.33*ku)*tu/3),
            no_overshoot = lambda ku,tu: (.2*ku, 2*(.2*ku)/tu, (.2*ku)*tu/3)
        )
        self.kp,self.ki,self.kd = converter[control_type.lower()](ku,tu)
