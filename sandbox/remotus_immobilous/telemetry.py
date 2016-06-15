import tables
import time
import krpc


class Telemetry(tables.IsDescription):
    time            = tables.Float64Col()
    altitude        = tables.Float64Col()
    vertical_speed  = tables.Float64Col()
    speed           = tables.Float64Col()
    pitch           = tables.Float64Col()
    heading         = tables.Float64Col()
    roll            = tables.Float64Col()
    con_pitch       = tables.Float64Col()
    con_yaw         = tables.Float64Col()
    con_roll        = tables.Float64Col()

def record(linkup):
    ksc = linkup.space_center
    vessel = ksc.active_vessel

    body = vessel.orbit.body
    surface_flight = vessel.flight(body.reference_frame)
    con = vessel.control

    altitude = linkup.add_stream(getattr, surface_flight, 'mean_altitude')
    vertical_speed = linkup.add_stream(getattr, surface_flight, 'vertical_speed')
    speed = linkup.add_stream(getattr, surface_flight, 'speed')
    pitch = linkup.add_stream(getattr, surface_flight, 'pitch')
    heading = linkup.add_stream(getattr, surface_flight, 'heading')
    roll = linkup.add_stream(getattr, surface_flight, 'roll')
    pitch = linkup.add_stream(getattr, surface_flight, 'pitch')
    heading = linkup.add_stream(getattr, surface_flight, 'heading')
    roll = linkup.add_stream(getattr, surface_flight, 'roll')
    con_pitch = linkup.add_stream(getattr, con, 'pitch')
    con_yaw = linkup.add_stream(getattr, con, 'yaw')
    con_roll = linkup.add_stream(getattr, con, 'roll')

    try:
        h5file = tables.open_file('telemetry.h5', mode='w', title='data')
        group = h5file.create_group('/', 'telemetry', 'Flight information')
        table = h5file.create_table(group, 'readout', Telemetry, 'Telemetry information')
        tel = table.row
        while True:
            tel['time'] = ksc.ut
            tel['altitude'] = altitude()
            tel['vertical_speed'] = vertical_speed()
            tel['speed'] = speed()
            tel['pitch'] = pitch()
            tel['heading'] = heading()
            tel['roll'] = roll()
            tel['con_pitch'] = con_pitch()
            tel['con_yaw'] = con_yaw()
            tel['con_roll'] = con_roll()
            tel.append()
            count -= 1
            time.sleep(0.01)
    finally:
        table.flush()
        h5file.close()

if __name__ == '__main__':
    linkup = krpc.connect('192.168.1.2', name='telemetry')
    record(linkup)
