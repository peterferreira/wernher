import unittest

import numpy as np
from ksp import Orbit, OrbitType, CelestialBody, plot_orbit_3d
from matplotlib import pyplot

arctan2 = np.arctan2
arccos = np.arccos
sqrt = np.sqrt
cos = np.cos

π = np.pi
deg = π/180
km = 1000

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

    def isclose(self,x,y,rtol=1e-4,atol=1e-4):
        msg = '{} != {}'.format(x,y)
        self.assertTrue(np.isclose(x,y,rtol=rtol,atol=atol),msg)

    def test_curtis_ex2_5(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                rotational_speed = 72.9217e-6,
                )
            o = Orbit(body=b)
            return o
        o = setup_orbit()
        self.isclose(o.body.stationary_radius,42164*km)
        o = setup_orbit()
        self.isclose(o.body.stationary_altitude,35786*km)
        o = setup_orbit()
        self.isclose(o.body.stationary_speed,3.0747*km)



    def test_curtis_ex2_7(self):

        # page 91 (Orbital Mech for Eng, 3rd ed)
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            oe = Bunch(
                body = b,
                periapsis_altitude = 400*km,
                apoapsis_altitude = 4000*km,
                )
            return Orbit(**oe)

        h = 60*60

        o = setup_orbit()
        self.isclose(o.eccentricity,0.2098)
        o = setup_orbit()
        self.isclose(o.specific_angular_momentum,57172*km**2)
        o = setup_orbit()
        self.isclose(o.speed_at_periapsis,8.435*km)
        o = setup_orbit()
        self.isclose(o.speed_at_apoapsis,5.509*km)
        o = setup_orbit()
        self.isclose(o.semi_major_axis,8578*km)
        o = setup_orbit()
        self.isclose(o.period,2.1963*h)
        o = setup_orbit()
        self.isclose(o.radius_at_average_true_anomaly,8387*km)
        o = setup_orbit()
        rθ = o.radius_at_average_true_anomaly
        θ = o.true_anomaly_from_periapsis(r=rθ)
        self.isclose(θ/deg,96.09)
        self.isclose(o.speed_at_radius(rθ),6.970*km)
        self.isclose(o.flight_path_angle_at_true_anomaly(θ)/deg,12.047)
        o = setup_orbit()
        self.isclose(o.true_anomaly_at_semi_minor_axis/deg,102.113)
        o = setup_orbit()
        self.isclose(o.flight_path_angle_max/deg,12.113)

    def test_curtis_ex2_7_aliases(self):

        # page 91 (Orbital Mech for Eng, 3rd ed)
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            oe = Bunch(
                body = b,
                pe_alt = 400*km,
                ap_alt = 4000*km,
                )
            return Orbit(**oe)

        o = setup_orbit()
        self.isclose(o.e,0.2098)
        o = setup_orbit()
        self.isclose(o.h,57172*km**2)
        o = setup_orbit()
        self.isclose(o.a,8578*km)

    def test_curtis_ex2_8(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )

            z1 = 1545*km
            θ1 = 126*deg
            z2 = 852*km
            θ2 = 58*deg
            return Orbit.from_altitude_at_true_anomaly(z1,θ1,z2,θ2,b)

        o = setup_orbit()
        self.isclose(o.e,0.08164)
        o = setup_orbit()
        self.isclose(o.h,54830*km**2)
        o = setup_orbit()
        self.isclose(o.pe,6974*km)
        o = setup_orbit()
        self.isclose(o.pe_alt,595.5*km)
        o = setup_orbit()
        self.isclose(o.a,7593*km)
        o = setup_orbit()
        self.isclose(o.T,1.8291*60*60)


    def test_curtis_ex2_9(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            o = Orbit(body=b)
            o.eccentricity = 1
            o.periapsis = 7000*km
            return o

        o = setup_orbit()
        self.isclose(o.specific_angular_momentum,74700*km**2)

        o = setup_orbit()
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
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            o = Orbit(body=b)
            o.radius_at_epoch = 14600*km
            o.speed_at_epoch = 8.6*km
            o.flight_path_angle_at_epoch = 50*deg
            return o

        o = setup_orbit()
        self.assertTrue(o.orbit_type is OrbitType.hyperbolic)

        o = setup_orbit()
        self.isclose(o.tangent_speed_at_epoch,5.528*km)
        o = setup_orbit()
        self.isclose(o.specific_angular_momentum,80708*km**2)
        o = setup_orbit()
        self.isclose(o.radial_speed_at_epoch,6.588*km)
        o = setup_orbit()
        self.isclose(o.eccentricity,1.3393)
        o = setup_orbit()
        self.isclose(o.true_anomaly_at_epoch,84.889*deg)

        o = setup_orbit()
        self.isclose(o.periapsis,6986*km)
        o = setup_orbit()
        self.isclose(o.semi_major_axis,-20590*km)
        o = setup_orbit()
        self.isclose(o.speed_at_infinity,sqrt(19.36)*km)
        o = setup_orbit()
        self.isclose(o.turn_angle,96.60*deg)
        o = setup_orbit()
        self.isclose(o.aiming_radius,18344*km)

    def test_curtis_ex2_11(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            o = Orbit(body=b)
            o.e = 0.3
            o.h = 60000*km**2
            o.θ0 = 120*deg
            return o

        o = setup_orbit()
        p,q = o.perifocal_position_at_epoch
        self.isclose(p,-5312.7*km)
        self.isclose(q,9201.9*km)

        o = setup_orbit()
        vp,vq = o.perifocal_velocity_at_epoch
        self.isclose(vp,-5.7533*km)
        self.isclose(vq,-1.3287*km)

    def test_curtis_ex2_12(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            o = Orbit(body=b)
            o.position_at_epoch = (7000*km, 9000*km, 0)
            o.velocity_at_epoch = (-5*km, 7*km, 0)
            return o

        o = setup_orbit()
        self.isclose(o.e, 1.1077)
        o = setup_orbit()
        self.isclose(o.h,94000*km**2)
        o = setup_orbit()
        self.isclose(o.pe,10517.5*km)
        o = setup_orbit()
        self.isclose(o.radius_at_epoch,11401.8*km)
        o = setup_orbit()
        self.isclose(o.θ0,31.5224*deg)
        o = setup_orbit()
        self.isclose(o.speed_at_epoch,8.6023*km)
        o = setup_orbit()
        self.isclose(o.h,94000*km**2)
        o = setup_orbit()
        self.isclose(o.e, 1.1077)

        o = setup_orbit()
        self.assertTrue(o.orbit_type is OrbitType.hyperbolic)

    def test_curtis_ex2_13(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            o = Orbit(body=b)
            o.position_at_epoch = (8182.4*km, -6865.9*km, 0)
            o.velocity_at_epoch = (0.47572*km, 8.8116*km, 0)
            return o

        o = setup_orbit()
        self.isclose(o.radius_at_epoch,10681*km) # 10861 in book (error)
        o = setup_orbit()
        self.isclose(o.speed_at_epoch,8.8244*km)
        o = setup_orbit()
        self.isclose(o.radial_speed_at_epoch,-5.2996*km)
        o = setup_orbit()
        self.isclose(o.h,75366*km**2)

        o = setup_orbit()
        θ0 = o.true_anomaly_at_epoch
        self.isclose(θ0,288.44*deg)

        r = o.radius_at_true_anomaly(θ0 + 120*deg)
        self.isclose(r,8378.8*km)

    def test_curtis_ex2_14(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            o = Orbit(body=b)
            o.position_at_epoch = (8182.4*km, -6865.9*km, 0)
            o.velocity_at_epoch = (0.47572*km, 8.8116*km, 0)
            return o

        o = setup_orbit()
        θ0 = o.true_anomaly_at_epoch

        x,y,z = o.position_at_true_anomaly(θ0 + 120*deg)
        self.isclose(x,1454.9*km)
        self.isclose(y,8251.6*km)
        self.isclose(z,0)

        o = setup_orbit()
        θ0 = o.true_anomaly_at_epoch

        vx,vy,vz = o.velocity_at_true_anomaly(θ0 + 120*deg)
        self.isclose(vx,-8.1323*km)
        self.isclose(vy,5.6785*km)
        self.isclose(vz,0)

    def test_curtis_ex3_1(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            o = Orbit(body=b)
            o.pe = 9600*km
            o.ap = 21000*km
            o.epoch = 0
            o.θ0 = 0
            return o

        o = setup_orbit()
        self.isclose(o.e,0.37255)
        o = setup_orbit()
        self.isclose(o.h,72472*km**2)
        o = setup_orbit()
        self.isclose(o.T,18834)

        o = setup_orbit()
        E = o.eccentric_anomaly_at_true_anomaly(120*deg)
        self.isclose(E,1.7281)
        M = o.mean_anomaly_at_eccentric_anomaly(E)
        self.isclose(M,1.3601)

        o = setup_orbit()
        t = o.time_at_true_anomaly(120*deg)
        self.isclose(t,4077)


    def test_curtis_ex3_2(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3,
                )
            o = Orbit(body=b)
            o.pe = 9600*km
            o.ap = 21000*km
            o.epoch = 0
            o.θ0 = 0
            return o

        t = 3*60*60

        o = setup_orbit()
        M = o.mean_anomaly_at_time(t)
        self.isclose(M,3.6029)
        o = setup_orbit()
        E = o.eccentric_anomaly_at_time(t)
        self.isclose(E,3.4794)
        o = setup_orbit()
        θ = o.true_anomaly_at_time(t)
        self.isclose(θ,193.18*deg)

    def test_curtis_ex3_3(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)
            o = Orbit(body=b,
                pe_alt = 500*km,
                ap_alt = 5000*km)
            return o

        o = setup_orbit()
        self.isclose(o.e,0.24649)
        o = setup_orbit()
        self.isclose(o.h,58458*km**2)
        o = setup_orbit()
        self.isclose(o.a,9128*km)
        o = setup_orbit()
        self.isclose(o.T,8679.1)

        o = setup_orbit()
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

        o.epoch = 0
        ta = o.time_at_true_anomaly(θa)
        tb = o.time_at_true_anomaly(θb)
        tc = o.time_at_true_anomaly(θc)
        td = o.time_at_true_anomaly(θd)

        self.isclose(tb,866.77)

        self.isclose(o.T - (ta-tb), 1733.5)
        self.isclose(td-tc, 2715.5)


    def test_curtis_ex3_4(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)
            return Orbit(body=b,
                orbit_type = OrbitType.parabolic,
                vpe = 10*km,
                θ0 = 0,
                t0 = 0,
                i = 0,
                Ω = 0,
                ω = 0)


        o = setup_orbit()
        self.isclose(o.e,1)
        o = setup_orbit()
        self.isclose(o.pe,7972*km)
        o = setup_orbit()
        h = o.h
        self.isclose(h,79720*km**2)
        o = setup_orbit()
        self.isclose(o.mean_anomaly_at_time(6*60*60),6.7737)

        o = setup_orbit()
        θ = o.true_anomaly_at_time(6*60*60)
        self.isclose(θ,144.75*deg)

        o = setup_orbit()
        # possible (rounding?) error in book: 86899*km
        self.isclose(o.radius_at_time(6*60*60),86976*km)


    def test_curtis_ex3_5(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)
            return Orbit(body=b,
                vpe = 15*km,
                pe_alt = 300*km,
                θ0 = 0,
                t0 = 0,
                i = 0,
                Ω = 0,
                ω = 0)

        θ = 100*deg

        o = setup_orbit()
        self.isclose(o.h,100170*km**2)
        o = setup_orbit()
        self.isclose(o.e,2.7696)
        o = setup_orbit()
        self.isclose(o.true_anomaly_at_infinity,111.17*deg)
        o = setup_orbit()
        self.isclose(o.radius_at_true_anomaly(θ),48497*km)
        o = setup_orbit()
        self.isclose(o.eccentric_anomaly_at_true_anomaly(θ),2.2927)
        o = setup_orbit()
        self.isclose(o.time_at_true_anomaly(θ),4141.4)

        o = setup_orbit()
        t0 = o.time_at_true_anomaly(θ)
        t = t0 + 3*60*60

        o = setup_orbit()
        self.isclose(o.mean_anomaly_at_time(t),40.690)
        o = setup_orbit()
        self.isclose(o.eccentric_anomaly_at_time(t),3.4631)
        o = setup_orbit()
        self.isclose(o.true_anomaly_at_time(t),107.78*deg)
        o = setup_orbit()
        self.isclose(o.radius_at_time(t),163180*km)
        o = setup_orbit()
        self.isclose(o.tangent_speed_at_time(t),0.61386*km)
        o = setup_orbit()
        self.isclose(o.radial_speed_at_time(t),10.494*km)
        o = setup_orbit()
        self.isclose(o.speed_at_time(t),10.512*km)
        o = setup_orbit()
        self.isclose(o.speed_at_infinity,10.277*km)


    def test_curtis_ex3_6(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)

            o = Orbit(body=b)
            o.epoch = 0
            o.true_anomaly_at_epoch = 30*deg
            o.radius_at_epoch = 10000*km
            o.speed_at_epoch = 10*km
            return o

        o = setup_orbit()
        self.isclose(o.e,1.4682)
        o = setup_orbit()
        self.isclose(o.h,95154*km**2)
        o = setup_orbit()
        self.isclose(o.radial_speed_at_epoch,3.0752*km)
        o = setup_orbit()
        self.isclose(o.eccentric_anomaly_at_time(0),0.23448)
        o = setup_orbit()
        self.isclose(o.a,-19655*km)
        o = setup_orbit()
        self.isclose(o.universal_anomaly_at_time(3600),128.51*km**0.5)



    def test_curtis_ex3_7(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)

            o = Orbit(body=b)
            o.epoch = 0
            o.position_at_epoch = (7000*km,-12124*km,0)
            o.velocity_at_epoch = (2.6679*km,4.6210*km,0)
            return o

        t = 60*60

        o = setup_orbit()
        self.isclose(o.radius_at_epoch, 14000*km)
        o = setup_orbit()
        self.isclose(o.speed_at_epoch, 5.3359*km)
        o = setup_orbit()
        self.isclose(o.radial_speed_at_epoch, -2.6679*km)
        o = setup_orbit()
        self.isclose(o.semi_major_axis, 13999*km)
        o = setup_orbit()
        self.isclose(1/o.semi_major_axis, 7.1429e-5 / km)

        o = setup_orbit()
        self.isclose(o.universal_anomaly_at_time(t), 253.53*km**0.5)

        o = setup_orbit()
        f,g,dfdχ,dgdχ = o.lagrange_coefficients_at_time(t)

        self.isclose(f,-0.54123)
        self.isclose(g,184.13)
        self.isclose(dfdχ,-0.00055298)
        self.isclose(dgdχ,-1.6593)

        o = setup_orbit()
        r = o.radius_at_time(t)
        self.isclose(r,8113.9*km)

        o = setup_orbit()
        x = o.position_at_time(t)
        self.isclose(x[0],-3297.8*km)
        self.isclose(x[1],7413.9*km)
        self.isclose(x[2],0)

        o = setup_orbit()
        v = o.velocity_at_time(t)
        self.isclose(v[0],-8.2977*km)
        self.isclose(v[1],-0.96404*km)
        self.isclose(v[2],0)

    def test_curtis_ex3_16(self):
        def setup_orbit():
            b = CelestialBody(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)
            o = Orbit(body=b)
            o.epoch = 0
            o.periapsis = 6600*km
            o.speed_at_periapsis = 1.2 * b.escape_speed(6600*km)
            return o

        o = setup_orbit()
        t0 = o.time_at_true_anomaly(-90*deg)
        t1 = o.time_at_true_anomaly(90*deg)
        self.isclose(t1-t0, 0.9992*60*60)

        o.θ0 = 0
        r = o.radius_at_time(24*60*60)
        self.isclose(r, 656610*km)

    def test_curtis_ex3_19(self):
        def setup_orbit():
            b = CelestialBody(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)
            o = Orbit(body=b)
            o.epoch = 0
            o.radius_at_epoch = 7200*km
            o.radial_speed_at_epoch = 1*km
            o.semi_major_axis = 10000*km
            o.i = 0
            o.Ω = 0
            o.ω = 0
            return o

        #o = setup_orbit()
        #F0 = o.eccentric_anomaly_at_epoch
        #F = o.eccentric_anomaly_at_time(o.epoch + 60*60)
        #print(np.sqrt(-o.a) * (F - F0))

        o = setup_orbit()
        χ = o.universal_anomaly_at_time(o.epoch + 60*60)
        #print(χ)




    def test_curtis_ex4_1(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)
            o = Orbit(body=b)

            o.epoch = 0
            o.position_at_epoch = (-5368*km, -1784*km, 3691*km)
            return o

        o = setup_orbit()
        self.isclose(o.declination_at_epoch, 33.12*deg)
        o = setup_orbit()
        self.isclose(o.right_ascension_at_epoch, 198.4*deg)


    def test_curtis_ex4_2(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)
            o = Orbit(body=b)

            o.epoch = 0
            o.position_at_epoch = (1600*km, 5310*km, 3800*km)
            o.velocity_at_epoch = (-7.350*km, 0.4600*km, 2.470*km)
            return o

        t = 3200

        o = setup_orbit()
        x = o.position_at_time(t)
        self.isclose(x[0],1091.3*km)
        self.isclose(x[1],-5199.4*km)
        self.isclose(x[2],-4480.6*km)

        o = setup_orbit()
        v = o.velocity_at_time(t)
        self.isclose(v[0],7.2284*km)
        self.isclose(v[1],1.9997*km)
        self.isclose(v[2],-0.46296*km)


    def test_curtis_ex4_3(self):

        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)
            o = Orbit(body=b)
            o.epoch = 0
            o.position_at_epoch = -6045*km, -3490*km, 2500*km
            o.velocity_at_epoch = -3.457*km, 6.618*km, 2.533*km
            return o

        o = setup_orbit()
        self.isclose(o.radius_at_epoch, 7414*km)
        o = setup_orbit()
        self.isclose(o.speed_at_epoch, 7.884*km)
        o = setup_orbit()
        self.isclose(o.radial_speed_at_epoch, 0.5575*km)
        o = setup_orbit()
        hx,hy,hz = o.specific_angular_momentum_vector
        self.isclose(hx,-25385*km**2)
        self.isclose(hy,6670*km**2)
        self.isclose(hz,-52070*km**2)
        o = setup_orbit()
        self.isclose(o.specific_angular_momentum,58310*km**2)
        o = setup_orbit()
        self.isclose(o.inclination, 153.25*deg)
        o = setup_orbit()
        nx,ny = o.node_line_vector
        self.isclose(nx,-6669.5*km**2)
        self.isclose(ny,-25385*km**2)
        o = setup_orbit()
        self.isclose(o.node_line,26247*km**2)
        o = setup_orbit()
        self.isclose(o.right_ascension, 255.3*deg)
        o = setup_orbit()
        ex,ey,ez = o.eccentricity_vector
        self.isclose(ex,-0.09160)
        self.isclose(ey,-0.1422)
        self.isclose(ez,0.02644)
        o = setup_orbit()
        self.isclose(o.eccentricity,0.1712)
        o = setup_orbit()
        self.isclose(o.argument_of_periapsis,20.07*deg)
        o = setup_orbit()
        self.isclose(o.true_anomaly_at_epoch,28.45*deg)
        o = setup_orbit()
        self.isclose(o.periapsis,7284*km)
        o = setup_orbit()
        self.isclose(o.apoapsis,10293*km)
        o = setup_orbit()
        self.isclose(o.semi_major_axis,8788*km)
        o = setup_orbit()
        self.isclose(o.period,2.2775*60*60)


    def test_curtis_ex4_7(self):
        def setup_orbit():
            b = Bunch(
                equatorial_radius = 6378*km,
                gravitational_parameter = 398600*km**3)
            o = Orbit(body=b,
                epoch = 0,
                h = 80000*km**2,
                e = 1.4,
                i = 30*deg,
                Ω = 40*deg,
                ω = 60*deg,
                θ0 = 30*deg)
            return o

        o = setup_orbit()
        x,y = o.perifocal_position_at_epoch
        self.isclose(x,6285*km)
        self.isclose(y,3628.6*km)

        o = setup_orbit()
        vx,vy = o.perifocal_velocity_at_epoch
        self.isclose(vx,-2.4913*km)
        self.isclose(vy,11.290*km)

        o = setup_orbit()
        x,y,z = o.position_at_epoch
        self.isclose(x,-4040*km)
        self.isclose(y,4815*km)
        self.isclose(z,3628.6*km)

        o = setup_orbit()
        vx,vy,vz = o.velocity_at_epoch
        self.isclose(vx,-10.386*km)
        self.isclose(vy,-4.772*km)
        self.isclose(vz,1.744*km)

    def test_pe_speed(self):
        b = Bunch(
            equatorial_radius = 6378*km,
            gravitational_parameter = 398600*km**3)

        ee = np.linspace(0,2,50)
        apo = []
        vpe = []
        for e in ee:
            o = Orbit(body=b,pe=7000*km,e=e)
            apo.append(o.apoapsis_altitude)
            vpe.append(o.speed_at_periapsis)
        apo = np.array(apo)/(km*km)
        vpe = np.array(vpe)/km

        '''
        from matplotlib import pyplot
        fig = pyplot.figure()
        ax = fig.add_subplot(1,1,1)
        axt = ax.twiny()
        ax.plot(ee,vpe,lw=3,color='blue')
        axt.plot(apo,vpe,lw=3,color='red')
        ax.set_xlabel('eccentricity',color='blue')
        axt.set_xlabel('apoapsis altitude (Mm)',color='red')
        ax.set_ylabel('speed at periapsis (Km/s)')
        ax.tick_params(axis='x', colors='blue')
        axt.tick_params(axis='x', colors='red')
        ax.spines['top'].set_color('red')
        ax.spines['bottom'].set_color('blue')
        pyplot.show()
        '''

    def test_elliptic_orbits_1(self):
        def setup_orbit():
            return Orbit(
                i  = 0.14,
                Ω  = 3.9,
                ω  = 0.2,
                e  = 0.1,
                a  = 8850000,
                M0 = 5,
                t0 = 34000000,
                body = CelestialBody(
                    name = 'kerbin',
                    equatorial_radius = 600000.0,
                    gravitational_parameter = 3531600035840.0,
                    rotational_speed = 0.0002908894093707204,
                ),
            )

        o = setup_orbit()
        #print('x0',o.radius_at_epoch,o.position_at_epoch)
        #print('v0',o.speed_at_epoch,o.velocity_at_epoch)
        #print('vr0',o.radial_speed_at_epoch)
        #print('γ',o.flight_path_angle_at_epoch)

        tt = np.linspace(o.epoch,o.epoch+1*6*60*60,50)
        χχexpect = np.sqrt(o.a) \
            * (o.eccentric_anomaly_at_time(tt) \
            - o.eccentric_anomaly_at_epoch)
        χχ = o.universal_anomaly_at_time(tt)

        #from matplotlib import pyplot
        #pyplot.plot(tt,χχexpect, color='blue')
        #pyplot.plot(tt[[0,-1]],χχexpect[[0,-1]], color='lightblue')
        #pyplot.plot(tt,χχ, color='red')
        #pyplot.plot(tt[[0,-1]],χχ[[0,-1]], color='pink')
        #pyplot.show()

        #plot_orbit_3d(o)
        #pyplot.show()



    def test_hyperbolic_orbits_1(self):
        o = Orbit(
            i  = 0.14170287439640022,
            Ω  = 3.978394514307273,
            ω  = 0.22281479333802098,
            e  = 1.077961355413604,
            a  = -8856935.204227254,
            M0 = -6.003136275330101,
            t0 = 34056451.522458464,
            body = CelestialBody(
                name = 'kerbin',
                equatorial_radius = 600000.0,
                gravitational_parameter = 3531600035840.0,
                rotational_speed = 0.0002908894093707204,
            ),
        )

        t = o.epoch
        tpe = o.time_to_periapsis_at_epoch
        npoints = 10
        tmin = t + tpe - 30
        tmax = t + tpe + 30

        tt = np.linspace(tmin,tmax,npoints)
        lat = o.latitude_at_time(tt) / deg
        lon = o.longitude_at_time(tt) / deg
        r = o.radius_at_time(tt)

        #print(list(zip(r,lat,lon)))

        #for out,r,lat,lon in zip(expected,r,lat,lon):
        #    out_r,out_lat,out_lon = out
        #    self.isclose(out_r,r)
        #    self.isclose(out_lat,lat)
        #    self.isclose(out_lon,lon)


    def test_hyperbolic_orbits_2(self):
        # same orbit pulled from KSP at different times
        # lots of jitter, so tolerance for radius at time
        # is quite large (1%)
        o1 = Orbit(
            i  = 0.14149227768205455,
            Ω  = 3.977254620789031,
            ω  = 0.22395653996553322,
            e  = 1.0779572208620696,
            a  = -8855744.039847286,
            M0 = -6.0078569863130475,
            t0 = 34056398.642449796,
            body = CelestialBody(
                name = 'Kerbin',
                equatorial_radius = 600000.0,
                gravitational_parameter = 3531600035840.0,
                rotational_speed = 0.0002908894093707204,
            ),
        )
        o2 = Orbit(
            i  = 0.14311811324451928,
            Ω  = 3.9859850089566558,
            ω  = 0.21660875587091438,
            e  = 1.07736727693924,
            a  = -8880527.593801336,
            M0 = -0.012499818155590287,
            t0 = 34140525.998688884,
            body = CelestialBody(
                name = 'Kerbin',
                equatorial_radius = 600000.0,
                gravitational_parameter = 3531600035840.0,
                rotational_speed = 0.0002908894093707204,
            ),
        )

        tpe1 = o1.time_to_periapsis_at_epoch
        tpe2 = o2.time_to_periapsis_at_epoch
        self.isclose(o1.epoch+tpe1, o2.epoch+tpe2)

        χexpect1 = np.sqrt(-o1.a) \
            * (o1.eccentric_anomaly_at_time(o1.epoch+tpe1) \
            - o1.eccentric_anomaly_at_epoch)
        χexpect2 = np.sqrt(-o2.a) \
            * (o2.eccentric_anomaly_at_time(o2.epoch+tpe2) \
            - o2.eccentric_anomaly_at_epoch)

        tt = np.linspace(o1.epoch,o1.epoch+1*6*60*60,50)
        χχexpect = np.sqrt(-o1.a) \
            * (o1.eccentric_anomaly_at_time(tt) \
            - o1.eccentric_anomaly_at_epoch)
        χχ = o1.universal_anomaly_at_time(tt)

        t0 = o1.epoch
        r0 = o1.radius_at_epoch
        vr0 = o1.radial_speed_at_epoch
        a = o1.semi_major_axis
        x0 = o1.position_at_epoch/km
        v0 = o1.velocity_at_epoch
        #print('t0',t0)
        #print('r0',r0)
        #print('vr0',vr0)
        #print('a',a)
        #print('x0',x0)
        #print('v0',v0)

        #from matplotlib import pyplot
        #pyplot.plot(tt,χχexpect, color='blue')
        #pyplot.plot(tt[[0,-1]],χχexpect[[0,-1]], color='lightblue')
        #pyplot.plot(tt,χχ, color='red')
        #pyplot.plot(tt[[0,-1]],χχ[[0,-1]], color='pink')
        #pyplot.show()

        #print('θpe1',o1.true_anomaly_at_time(o1.epoch+tpe1))
        #print('θpe2',o2.true_anomaly_at_time(o2.epoch+tpe2))

        #print('rpeθ1',o1.radius_at_true_anomaly(
        #o1.true_anomaly_at_time(o1.epoch+tpe1)))
        #print('rpeθ2',o2.radius_at_true_anomaly(
        #o2.true_anomaly_at_time(o2.epoch+tpe2)))

        f1,g1,df1,dg1 = o1.lagrange_coefficients_at_time(o1.epoch+tpe1)
        f2,g2,df2,dg2 = o2.lagrange_coefficients_at_time(o2.epoch+tpe2)
        x1 = o1.position_at_lagrange_coefficients(f1,g1)
        x2 = o2.position_at_lagrange_coefficients(f2,g2)
        #print('rpe1(χexp)',np.sqrt(sum(x1**2)))
        #print('rpe2(χexp)',np.sqrt(sum(x2**2)))

        self.isclose(o2.periapsis, o2.radius_at_time(o2.epoch+tpe2))
        self.isclose(o1.periapsis, o1.radius_at_time(o1.epoch+tpe1))

        self.isclose(o1.periapsis, o2.periapsis,rtol=1e-2)
        self.isclose(o1.radius_at_time(o1.epoch+tpe1),
                     o2.radius_at_time(o2.epoch+tpe2),
                     rtol=0.01)

if __name__ == '__main__':
    unittest.main()
