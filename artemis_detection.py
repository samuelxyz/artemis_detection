import csv
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.ticker as tck
import numpy as np

# what orbits do i want?
# let's say jacobi constant at least 3 and pd at least 11. note lower J means more energy (intuitively J has a -v^2 term)
# note that stability index =1 is marginally stable and >1 is unstable at least in cr3bp. in ephemeris force models the ones tightest around moon and lpoint seem more consistent
# good paper here https://www.researchgate.net/publication/339028637_Stationkeeping_and_Transfer_Trajectory_Design_for_Spacecraft_in_Cislunar_Space

def load_eph_file(fname: str):
    res = []
    with open(fname, mode='r') as eph_file:
        csv_reader = csv.reader(eph_file)
        reading_eph = False
        for row in csv_reader:
            if not row: continue
            if reading_eph:
                if row[0] == '$$EOE': # end of ephemeris
                    break
                else:
                    # put in the results
                    jd, date_str, x, y, z, v_x, v_y, v_z, lt, range, range_rate, extra = row
                    entry = (
                        float(jd),
                        date_str,
                        (float(x), float(y), float(z)),
                        (float(v_x), float(v_y), float(v_z)),
                        float(lt),
                        float(range),
                        float(range_rate)
                    )
                    res.append(entry)
            else:
                if row[0] == '$$SOE': # start of ephemeris
                    reading_eph = True
                    continue
    
    return np.asarray(res, dtype=[
        ('jd', float),
        ('date_str', 'U31'),
        ('pos', float, (3,)),
        ('vel', float, (3,)),
        ('light', float),
        ('range', float),
        ('range_rate', float),
    ])

def at_time(eph, julian_date, key):
    """Get entry from the given ephemeris at the given time. 
    Uses linear interpolation between adjacent ephemeris entries for now."""
    i = np.searchsorted(eph['jd'], julian_date)
    before = eph[i-1]
    after = eph[i]
    x = (julian_date - before['jd'])/(after['jd'] - before['jd'])
    return before[key]*(1-x) + after[key]*x
    

# orion_ICRF = load_eph_file("horizons_results_orion_barycenter_1day_csv.txt")
orion_ICRF = load_eph_file("horizons_results_orion_barycenter_1hour_csv.txt")
earth_ICRF = load_eph_file("horizons_results_earth_barycenter_1hour_csv.txt")
moon_ICRF = load_eph_file("horizons_results_moon_barycenter_1hour_csv.txt")

def to_rotating_frame(eph: np.ndarray):
    """Convert the ephemeris to a coordinate system aligned along the earth-moon axis 
    and moon orbit plane"""
    res = eph.copy()
    # For now we assume eph has the same time samples as moon_ICRF
    rot = np.zeros((3, 3)) 
    for (elem, moon) in zip(res, moon_ICRF):
        # x axis is toward the moon position
        rot[:, 0] = moon['pos']/np.linalg.norm(moon['pos'])
        # z axis is along the moon angular momentum vector
        zdir = np.linalg.cross(moon['pos'], moon['vel'])
        rot[:, 2] = zdir/np.linalg.norm(zdir)
        rot[:, 1] = np.linalg.cross(rot[:, 2], rot[:, 0])
        elem['pos'] = np.linalg.solve(rot, elem['pos'])
        elem['vel'] = np.linalg.solve(rot, elem['vel'])
    return res

def to_ICRF(eph: np.ndarray):
    """The inverse of to_rotating_frame(eph)"""
    res = eph.copy()
    # For now we assume eph has the same time samples as moon_ICRF
    rot = np.zeros((3, 3)) 
    for (elem, moon) in zip(res, moon_ICRF):
        # x axis is toward the moon position
        rot[:, 0] = moon['pos']/np.linalg.norm(moon['pos'])
        # z axis is along the moon angular momentum vector
        zdir = np.linalg.cross(moon['pos'], moon['vel'])
        rot[:, 2] = zdir/np.linalg.norm(zdir)
        rot[:, 1] = np.linalg.cross(rot[:, 2], rot[:, 0])
        elem['pos'] = rot @ elem['pos']
        elem['vel'] = rot @ elem['vel']
    return res

orion_ROT = to_rotating_frame(orion_ICRF)
earth_ROT = to_rotating_frame(earth_ICRF)
moon_ROT = to_rotating_frame(moon_ICRF)

def load_orbit(fname: str):
    res = []
    with open(fname) as orbit_file:
        csv_reader = csv.reader(orbit_file)
        for row in csv_reader:
            time, x, y, z, v_x, v_y, v_z = row
            if time == 'Time (TU)': continue # skip header row
            res.append((
                time,
                (x, y, z),
                (v_x, v_y, v_z),
            ))
    orbit = np.asarray(res, dtype=[
        ('time', float),
        ('pos', float, (3,)),
        ('vel', float, (3,)),
    ])
    
    # dimensionalize
    # Actually let's dimensionalize with the current earth-moon distance 
    # instead of a circular orbit constant
    # LU = 389703 # length unit in km
    TU = 382981 # time unit in seconds
    orbit['time'] *= TU
    # orbit['pos'] *= LU
    # orbit['vel'] *= LU/TU

    orbit_ROT = moon_ROT.copy() # easiest way to create new ephemeris arrary
    # some items will be wrong like light time and range etc
    # but oh well just dont use those

    for (earth, moon, orb_entry) in zip(earth_ICRF, moon_ICRF, orbit_ROT):
        # for now just look up the closest entry from orbit_ROT
        # i could implement interpolation later or something similar
        mission_time = (earth['jd'] - earth_ICRF[0]['jd']) * 86400
        period = orbit[-1]['time']
        orbit_i = np.searchsorted(orbit['time'], mission_time % period)
        LU = earth['range'] + moon['range']
        orb_entry['pos'] = orbit['pos'][orbit_i] * LU
        orb_entry['vel'] = orbit['vel'][orbit_i] * LU/TU
        # TODO other parameters like range, range rate, light, etc could be corrected
    
    return orbit_ROT
# orbit_ROT = load_orbit('cr3bp_L1_northernhalo_3.17387.csv')
# orbit_ROT = load_orbit('cr3bp_L1_northernhalo_3.07028.csv')
orbit_ROT = load_orbit('cr3bp_L2_southernhalo_3.12603.csv')
orbit_ICRF = to_ICRF(orbit_ROT)

anim_time_rate = 3 # skip some frames for a fast-forward effect

def animate(ephs, labels, colors, title=''):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.set(xlabel='X (km)', ylabel='Y (km)', zlabel='Z (km)')
    for ax_ in (ax.xaxis, ax.yaxis, ax.zaxis):
        ax_.set_major_formatter(tck.EngFormatter())

    for eph, label, color in zip(ephs, labels, colors):
        ax.plot(*eph['pos'].T, label=label, color=color)

    data = tuple(zip(*(eph[0]['pos'] for eph in ephs)))
    points = ax.scatter(*data, c=colors, depthshade=False)
    # time_text = ax.text(-4e5, -4e5, 0, '')

    def update_anim(i):
        i = i*anim_time_rate
        data = tuple(zip(*(eph[i]['pos'] for eph in ephs)))
        # points.set_offsets(data)
        points._offsets3d = data # bruh
        points.stale = True
        # time_text.set_text(f'Time: {ephs[0][i]["date_str"]}')
        ax.set_title(f'{title}\nTime: {ephs[0][i]["date_str"]}')
        # return (points, time_text)
        return (points,)

    ax.set_aspect('equal')
    ax.legend()

    return anim.FuncAnimation(fig, update_anim, int(len(orion_ICRF['jd'])/anim_time_rate), interval=10, blit=False)

anim_ = animate((orion_ICRF, earth_ICRF, moon_ICRF, orbit_ICRF), 
                    ('Artemis 1', 'Earth', 'Moon', 'L1 Halo Observer'), 
                    ('C0', 'green', 'gray', 'red'), 
                    'Artemis 1 Mission (ICRF Barycenter Frame)')
# anim_ = animate((orion_ROT, earth_ROT, moon_ROT, orbit_ROT), 
#                    ('Artemis 1', 'Earth', 'Moon', 'L1 Halo Observer'), 
#                    ('C0', 'green', 'gray', 'red'), 
#                    'Artemis 1 Mission (Earth-Moon Rotating Barycenter Frame)')

if False:
    writer = anim.PillowWriter(fps=20,
                                    metadata=dict(artist='Samuel Tan'),
                                    bitrate=1800)
    anim_.save('ICRF.gif', writer=writer)

def simulate_surveillance(imaging_interval, fov, sensing_range):
    """Simulate a random pointing mission where `orbit_ICRF` attempts to image `orion_ICRF`.
    
    # Parameters:
    imaging_interval: float
        the number of seconds between imaging attempts
    fov: float
        the imager field of view in square degrees
    sensing range: float
        sensing range in km

    # Returns:
    hits: float
        the expected number of successful sensing attempts"""
    
    hits = 0

    for jd in np.arange(earth_ICRF[0]['jd'], earth_ICRF[-1]['jd'], imaging_interval/86400, dtype=float):
        sensor_pos = at_time(orbit_ICRF, jd, 'pos')
        orion_pos = at_time(orion_ICRF, jd, 'pos') - sensor_pos
        moon_pos = at_time(moon_ICRF, jd, 'pos') - sensor_pos
        moon_exclusion_zone = np.radians(0)

        orion_dist = np.linalg.norm(orion_pos)
        moon_dist = np.linalg.norm(moon_pos)

        if orion_dist > sensing_range: 
            continue
        if np.dot(moon_pos, orion_pos)/orion_dist/moon_dist > np.cos(moon_exclusion_zone):
            continue # orion too close to the moon, assume you cant look there
            # Could do a similar one for the sun in the future

        # We could do an actual random orientation selection and see if it matches up
        # But for now we could also just do
        whole_sky = 41253 # square degrees
        hit_chance = fov / whole_sky
        hits += hit_chance

    return hits

print(f'Hits: {simulate_surveillance(60, 5*5, 150e3)}')

plt.show()