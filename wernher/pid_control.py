import numpy as np

class Controller(object):
    '''Single Axis PID Controller
    goals:
        self-tuning
        edge-case switching to PD or PI
    '''

    def __init__(self,set_point=0,kp=1,ki=0,kd=0,t0=0):
        # set point
        self.set_point = set_point

        # control constants
        self.kp = kp
        self.ki = ki
        self.kd = kd

        # time and proportional term of previous call, integral
        self.t0 = t0
        self.I  = 0
        self.P0 = 0

        # response value of previous call
        self.c = 0

    def __call__(self,x,t):
        # just return previous value if no time has passed
        if np.isclose(t,self.t0):
            return self.c

        # bring instance variables into local scope
        xset = self.set_point
        kp = self.kp
        ki = self.ki
        kd = self.kd

        if all([np.isclose(k,0) for k in [kp,ki,kd]]):
            self.t0 = t
            return xset

        # bring instance variables into local scope
        t0 = self.t0
        I  = self.I
        P0 = self.P0

        # calculate PID control
        Δt = t - t0
        P = xset - x
        I += P * Δt
        D = (P - P0) / Δt

        c = kp*P + ki*I + kd*D

        # save parameters to class instance
        self.t0 = t
        self.I  = I
        self.P0 = P
        self.c  = c

        return c
