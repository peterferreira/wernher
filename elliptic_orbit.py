from .locked_property import locked_property, LockError

class EllipticOrbit(object):

    def __init__(self,**kwargs):
        for k,v in kwargs.items():
            setattr(self,k,v)

    ### Specific Angular Momentum (h) vs Semi-Major Axis (a)
    @locked_property
    def specific_angular_momentum(self):
        e = self.eccentricity
        a = self.semi_major_axis
        μ = self.body.gravitational_parameter
        h = sqrt(a * μ * (1 - e**2))
        return h

    @locked_property
    def semi_major_axis(self):
        h = self.specific_angular_momentum
        e = self.eccentricity
        μ = self.body.gravitational_parameter
        a = sqrt(h * μ * (1 - e**2))
        return a

    ### Mean Anomaly at Epoch (M0) vs True Anomaly at Epoch (θ0)
    @locked_property
    def mean_anomaly_at_epoch(self):
        e = self.eccentricity
        θ0 = self.true_anomaly_at_epoch
        M0 = 2 * arctan2(sqrt(1+e)*cos(θ0/2),sqrt(1-e)*sin(θ0/2)) \
                - ((e*sqrt(1-e**2)*sin(θ0)) / (1+e*cos(θ0)))
        return M0

    @locked_property
    def true_anomaly_at_epoch(self):
        θ0 = self.true_anomaly()
        return θ0

    @locked_property
    def mean_motion(self):
        e = self.eccentricity
        a = self.semi_major_axis
        μ = self.body.gravitational_parameter
        n = sqrt(μ / a**3)
        return n

    def mean_anomaly(self,t=None):
        M0 = self.mean_anomaly_at_epoch
        if t is None:
            return M0
        else:
            t0 = self.epoch
            n = self.mean_motion
            return M0 + n * (t - t0)

    def eccentric_anomaly(self,t=None):
        e = self.eccentricity
        M = self.mean_anomaly(t)
        def fn(E,e,M):
            return E - e * sin(E) - M
        if hasattr(t,'__iter__'):
            return np.array([opt.newton(fn, π, args=(e,m)) for m in M])
        else:
            return opt.newton(fn, π, args=(e,M))

    def true_anomaly(self,t=None):
        e = self.eccentricity
        E = self.eccentric_anomaly(t)
        return π - 2 * arctan2(sqrt(1-e) * cos(E/2),
                               sqrt(1+e) * sin(E/2))
