import numpy as np

inf = np.inf
isclose = np.isclose

class Controller(object):
    '''Single Axis PID Controller'''
    def __init__(self, kp=1, ki=0, kd=0,
                 set_point=0, deadband=0.001, cmin=-1, cmax=1,
                 t0=0):
        # set point and control constants
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.set_point = set_point
        self.deadband = deadband
        self.cmin = cmin
        self.cmax = cmax

        # time of previous call, running integral and
        # proportional term of previous call
        self.t0 = t0
        self.kiI = 0
        self.P0 = 0

        # response value of previous call
        self.c = 0

    def __call__(self,x,t):
        # if parameters are all zero or None, return zero
        if not any([self.kp,self.ki,self.kd]):
            self.t0 = t
            return 0

        Δt = t - self.t0

        # return previous value if no time has passed
        if isclose(Δt, 0):
            return self.c

        P = self.set_point - x

        if abs(P) < self.deadband:
            kpP = 0
        else:
            kpP = self.kp * P

        if Δt > self.td:
            kdD = 0
        else:
            ΔP = P - self.P0
            D = ΔP / Δt
            kdD = self.kd*D

        if Δt > self.ti:
            kiI = 0
        else:
            if self.ki > 0:
                if not (self.cmin < kpP < self.cmax):
                    kiI = 0
                else:
                    kiI = self.kiI + self.ki * P * Δt
                    kiI = min(max(kiI,self.cmin),self.cmax)

        self.c = kpP + kiI + kdD

        # clip output to specified limits
        self.c = min(max(self.c,self.cmin),self.cmax)

        # save parameters to class instance
        self.t0 = t
        self.kiI = kiI
        self.P0 = P

        return self.c

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
