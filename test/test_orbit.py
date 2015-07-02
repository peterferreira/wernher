import unittest

import numpy as np
from ksp.orbit2 import Orbit, OrbitType


arctan2 = np.arctan2
arccos = np.arccos
sqrt = np.sqrt
cos = np.cos

π = np.pi
deg = π/180
km = 1000
h = 60*60

class Bunch(dict):
    def __contains__(self, k):
        try:
            return dict.__contains__(self, k) or hasattr(self, k)
        except:
            return False
    def __getattr__(self, k):
        try:
            return object.__getattribute__(self, k)
        except AttributeError:
            try:
                return self[k]
            except KeyError:
                raise AttributeError(k)
    def __setattr__(self, k, v):
        try:
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                self[k] = v
            except:
                raise AttributeError(k)
        else:
            object.__setattr__(self, k, v)
    def __delattr__(self, k):
        try:
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                del self[k]
            except KeyError:
                raise AttributeError(k)
        else:
            object.__delattr__(self, k)


class TestOrbit(unittest.TestCase):

    def isclose(self,x,y):
        msg = '{} != {}'.format(x,y)
        self.assertTrue(np.isclose(x,y,rtol=1e-4,atol=1e-4),msg)

    def test_curtis_ex2_5(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            rotational_speed = 72.9217e-6,
            )
        o = Orbit(body=b)
        self.isclose(o.body.stationary_radius,42164*km)
        self.isclose(o.body.stationary_altitude,35786*km)
        self.isclose(o.body.stationary_speed,3.0747*km)



    def test_curtis_ex2_7(self):

        # page 91 (Orbital Mech for Eng, 3rd ed)
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        oe = Bunch(
            body = b,
            periapsis_altitude = 400*km,
            apoapsis_altitude = 4000*km,
            )

        o = Orbit(**oe)
        self.isclose(o.eccentricity,0.2098)

        o = Orbit(**oe)
        self.isclose(o.specific_angular_momentum,57172*km**2)

        o = Orbit(**oe)
        self.isclose(o.speed_at_periapsis,8.435*km)

        o = Orbit(**oe)
        self.isclose(o.speed_at_apoapsis,5.509*km)

        o = Orbit(**oe)
        self.isclose(o.semi_major_axis,8578*km)

        o = Orbit(**oe)
        self.isclose(o.period,2.1963*h)

        o = Orbit(**oe)
        self.isclose(o.radius_at_average_true_anomaly,8387*km)

        o = Orbit(**oe)
        rθ = o.radius_at_average_true_anomaly
        θ = o.true_anomaly_from_periapsis(r=rθ)
        self.isclose(θ/deg,96.09)
        self.isclose(o.speed_at_radius(rθ),6.970*km)
        self.isclose(o.flight_path_angle_at_true_anomaly(θ)/deg,12.047)
        self.isclose(o.true_anomaly_at_semi_minor_axis/deg,102.113)
        self.isclose(o.flight_path_angle_max/deg,12.113)

    def test_curtis_ex2_7_aliases(self):

        # page 91 (Orbital Mech for Eng, 3rd ed)
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        oe = Bunch(
            body = b,
            pe_alt = 400*km,
            ap_alt = 4000*km,
            )

        o = Orbit(**oe)
        self.isclose(o.e,0.2098)

        o = Orbit(**oe)
        self.isclose(o.h,57172*km**2)

        o = Orbit(**oe)
        self.isclose(o.a,8578*km)

    def test_curtis_ex2_8(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )

        z1 = 1545*km
        θ1 = 126*deg
        z2 = 852*km
        θ2 = 58*deg

        o = Orbit.from_altitude_at_true_anomaly(z1,θ1,z2,θ2,b)

        self.isclose(o.e,0.08164)
        self.isclose(o.h,54830*km**2)
        self.isclose(o.pe,6974*km)
        self.isclose(o.pe_alt,595.5*km)
        self.isclose(o.a,7593*km)
        self.isclose(o.T,1.8291*h)


    def test_curtis_ex2_9(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        o = Orbit(body=b)
        o.eccentricity = 1
        o.periapsis = 7000*km
        self.isclose(o.specific_angular_momentum,74700*km**2)

        r1 = 8000*km
        r2 = 16000*km

        θ1 = o.true_anomaly_at_radius(r1)
        θ2 = o.true_anomaly_at_radius(r2)

        self.isclose(θ1/deg,41.41)
        self.isclose(θ2/deg,97.18)

        Δθ = θ2 - θ1
        d = sqrt(r1**2 + r2**2 - 2*r1*r2*cos(Δθ))
        self.isclose(d,13266.5*km)

    def test_curtis_ex2_10(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        o = Orbit(body=b)
        o.radius_at_epoch = 14600*km
        o.speed_at_epoch = 8.6*km
        o.flight_path_angle_at_epoch = 50*deg

        self.assertTrue(o.orbit_type is OrbitType.hyperbolic)

        self.isclose(o.tangent_speed_at_epoch,5.528*km)
        self.isclose(o.specific_angular_momentum,80708*km**2)
        self.isclose(o.radial_speed_at_epoch,6.588*km)
        self.isclose(o.eccentricity,1.3393)
        self.isclose(o.true_anomaly_at_epoch,84.889*deg)

        self.isclose(o.periapsis,6986*km)
        self.isclose(o.semi_major_axis,-20590*km)
        self.isclose(o.speed_at_infinity,sqrt(19.36)*km)
        self.isclose(o.turn_angle,96.60*deg)
        self.isclose(o.aiming_radius,18344*km)

    def test_curtis_ex2_11(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        o = Orbit(body=b)
        o.e = 0.3
        o.h = 60000*km**2
        o.θ0 = 120*deg

        p,q = o.perifocal_position_at_epoch
        self.isclose(p,-5312.7*km)
        self.isclose(q,9201.9*km)

        vp,vq = o.perifocal_velocity_at_epoch
        self.isclose(vp,-5.7533*km)
        self.isclose(vq,-1.3287*km)

    def test_curtis_ex2_12(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        o = Orbit(body=b)
        o.perifocal_position_at_epoch = (7000*km, 9000*km)
        o.perifocal_velocity_at_epoch = (-5*km, 7*km)

        self.isclose(o.h,94000*km**2)
        self.isclose(o.θ0,52.125*deg)
        self.isclose(o.e,1.538)

        self.assertTrue(o.orbit_type is OrbitType.hyperbolic)

    def test_curtis_ex2_13(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        o = Orbit(body=b)
        o.position_at_epoch = (8182.4*km, -6865.9*km)
        o.velocity_at_epoch = (0.47572*km, 8.8116*km)

        self.isclose(o.radius_at_epoch,10681*km) # 10861 in book (error)
        self.isclose(o.speed_at_epoch,8.8244*km)
        self.isclose(o.radial_speed_at_epoch,-5.2996*km)
        self.isclose(o.h,75366*km**2)

        θ0 = o.true_anomaly_at_epoch
        self.isclose(θ0,288.44*deg)

        r = o.radius_at_true_anomaly(θ0 + 120*deg)
        self.isclose(r,8378.8*km)

    def test_curtis_ex2_14(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        o = Orbit(body=b)
        o.position_at_epoch = (8182.4*km, -6865.9*km)
        o.velocity_at_epoch = (0.47572*km, 8.8116*km)

        θ0 = o.true_anomaly_at_epoch

        x,y = o.position_at_true_anomaly(θ0 + 120*deg)
        vx,vy = o.velocity_at_true_anomaly(θ0 + 120*deg)
        self.isclose(x,1454.9*km)
        self.isclose(y,8251.6*km)
        self.isclose(vx,-8.1323*km)
        self.isclose(vy,5.6785*km)

    def test_curtis_ex3_1(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        o = Orbit(body=b)
        o.pe = 9600*km
        o.ap = 21000*km
        o.epoch = 0
        o.θ0 = 0

        self.isclose(o.e,0.37255)
        self.isclose(o.h,72472*km**2)
        self.isclose(o.T,18834)

        E = o.eccentric_anomaly_at_true_anomaly(120*deg)
        self.isclose(E,1.7281)
        M = o.mean_anomaly_at_eccentric_anomaly(E)
        self.isclose(M,1.3601)
        t = o.time_at_true_anomaly(120*deg)
        self.isclose(t,4077)


    def test_curtis_ex3_2(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3,
            )
        o = Orbit(body=b)
        o.pe = 9600*km
        o.ap = 21000*km
        o.epoch = 0
        o.θ0 = 0

        t = 3*60*60
        M = o.mean_anomaly_at_time(t)
        self.isclose(M,3.6029)
        E = o.eccentric_anomaly_at_time(t)
        self.isclose(E,3.4794)
        θ = o.true_anomaly_at_time(t)
        self.isclose(θ,193.18*deg)

    def test_curtis_ex3_3(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3)
        o = Orbit(body=b,
            pe_alt = 500*km,
            ap_alt = 5000*km)

        self.isclose(o.e,0.24649)
        self.isclose(o.h,58458*km**2)
        self.isclose(o.a,9128*km)
        self.isclose(o.T,8679.1)

        A = o.e
        B = -(1-o.e**2) * o.a / o.body.equatorial_radius
        C = -1
        D = arctan2(B,A)
        E = arccos((C/A) * cos(arctan2(B,A)))

        θb = D + E
        θc = D - E
        θa = 2*π - θb
        θd = 2*π - θc

        self.isclose(θb,57.423*deg)
        self.isclose(θc,-216.64*deg)

        self.isclose(o.eccentric_anomaly_at_true_anomaly(θb),0.80521)
        self.isclose(o.mean_anomaly_at_true_anomaly(θb),0.62749)

        ta = o.time_at_true_anomaly(θa)
        tb = o.time_at_true_anomaly(θb)
        tc = o.time_at_true_anomaly(θc)
        td = o.time_at_true_anomaly(θd)

        self.isclose(tb,866.77)

        self.isclose(o.T - (ta-tb), 1733.5)
        self.isclose(td-tc, 2715.5)


    def test_curtis_ex3_4(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3)
        o = Orbit(body=b,
            orbit_type = OrbitType.parabolic,
            vpe = 10*km,
            θ0 = 0,
            t0 = 0)

        self.isclose(o.e,1)
        self.isclose(o.pe,7972*km)
        h = o.h
        self.isclose(h,79720*km**2)
        self.isclose(o.mean_anomaly_at_time(6*60*60),6.7737)

        θ = o.true_anomaly_at_time(6*60*60)
        self.isclose(θ,144.75*deg)

        # possible (rounding?) error in book: 86899*km
        self.isclose(o.radius_at_time(6*60*60),86976*km)


    def test_curtis_ex3_5(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3)
        o = Orbit(body=b,
            vpe = 15*km,
            pe_alt = 300*km,
            θ0 = 0,
            t0 = 0)

        θ = 100*deg

        self.isclose(o.h,100170*km**2)
        self.isclose(o.e,2.7696)
        self.isclose(o.true_anomaly_at_infinity,111.17*deg)
        self.isclose(o.radius_at_true_anomaly(θ),48497*km)
        self.isclose(o.eccentric_anomaly_at_true_anomaly(θ),2.2927)
        self.isclose(o.time_at_true_anomaly(θ),4141.4)

        t0 = o.time_at_true_anomaly(θ)
        t = t0 + 3*60*60

        self.isclose(o.mean_anomaly_at_time(t),40.690)
        self.isclose(o.eccentric_anomaly_at_time(t),3.4631)
        self.isclose(o.true_anomaly_at_time(t),107.78*deg)
        self.isclose(o.radius_at_time(t),163180*km)
        self.isclose(o.tangent_speed_at_time(t),0.61386*km)
        self.isclose(o.radial_speed_at_time(t),10.494*km)
        self.isclose(o.speed_at_time(t),10.512*km)
        self.isclose(o.speed_at_infinity,10.277*km)




if __name__ == '__main__':
    unittest.main()
