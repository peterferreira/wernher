import numpy as np

inf = np.inf
isclose = np.isclose

class Controller(object):
    '''Single Axis PID Controller
    goals:
        self-tuning
        edge-case switching to PD or PI
    '''

    def __init__(self,set_point=0,kp=1,ki=0,kd=0,t0=0):
        # set point and control constants
        self.set_point = set_point
        self.kp = kp
        self.ki = ki
        self.kd = kd

        # time of previous call, running integral and
        # proportional term of previous call
        self.t0 = t0
        self.I  = 0
        self.P0 = 0

        # limits
        self.min = -inf
        self.max =  inf

        # response value of previous call
        self.c = 0

    def __call__(self,x,t):
        # return previous value if no time has passed
        if isclose(t - self.t0, 0):
            return self.c

        # bring instance variables into local scope
        xset = self.set_point
        kp = self.kp
        ki = self.ki
        kd = self.kd

        # if parameters are all zero or None, return set point
        if not all([kp,ki,kd]):
            self.t0 = t
            return xset

        # bring instance variables into local scope
        t0 = self.t0
        I  = self.I
        P0 = self.P0
        ti = self.ti

        # calculate PID terms
        Δt = t - t0
        P = xset - x
        ΔP = P - P0
        D = ΔP / Δt

        # freeze integral for a small time on
        # a large disturbance
        if abs(kp*ΔP) > 0.5*(self.max - self.min):
            self._t0_freeze_I = t
        else:
            try:
                if (t - self._t0_freeze_I) > ti:
                    del self._t0_freeze_I
                    I += P * Δt
            except AttributeError:
                I += P * Δt

        # clip proportional gain and turn off
        # integral term if P*km is out of the
        # control range
        if not (self.min < kp*P < self.max):
            P = min(max(P, self.min/kp), self.max/kp)
            I = 0
        else:
            I = min(max(I, self.min/ki), self.max/ki)

        c = kp*P + ki*I + kd*D

        # clip output to specified limits
        c = min(max(c, self.min), self.max)

        # save parameters to class instance
        self.t0 = t
        self.I  = I
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

    def ziegler_nichols(ku,tu,control_type='pid'):
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
