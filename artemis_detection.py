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

def load_orbit(fname: str):
    res = []
    with open(fname) as orbit_file:
        csv_reader = csv.reader(orbit_file)
        for row in csv_reader:
            time, x, y, z, v_x, v_y, v_z = row
            if time == 'Time (TU)': continue
            res.append({
                time,
                (x, y, z),
                (v_x, v_y, v_z),
            })
    orbit = np.asarray(res, dtype=[
        ('time', float),
        ('pos', float, (3,)),
        ('vel', float, (3,)),
    ])
    
    # dimensionalize
    LU = 389703 # length unit in km
    TU = 382981 # time unit in seconds
    orbit['time'] *= TU
    orbit['pos'] *= LU
    orbit['vel'] *= LU/TU
    
    return orbit

def at_time(eph, julian_date):
    '''Get entry from the given ephemeris at the given time'''
    # Maybe at some future point we can do a fancier interpolation or propagation
    # For now let's just do nearest point lookup
    # Actually do we even need this function? nvm let's leave it for now
    pass

def plot_ephem(eph, **kwargs):
    # posns = [entry["pos"] for entry in eph]
    # posns_zipped = tuple(zip(*posns))
    # ax.plot(*posns_zipped, **kwargs)
    ax.plot(*eph['pos'].T, **kwargs)

# orion_ICRF = load_eph_file("horizons_results_orion_barycenter_1day_csv.txt")
orion_ICRF = load_eph_file("horizons_results_orion_barycenter_1hour_csv.txt")
earth_ICRF = load_eph_file("horizons_results_earth_barycenter_1hour_csv.txt")
moon_ICRF = load_eph_file("horizons_results_moon_barycenter_1hour_csv.txt")

def to_rotating_frame(eph: np.ndarray):
    '''Convert the ephemeris to a coordinate system aligned along the earth-moon axis 
    and moon orbit plane'''
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

orion_ROT = to_rotating_frame(orion_ICRF)
earth_ROT = to_rotating_frame(earth_ICRF)
moon_ROT = to_rotating_frame(moon_ICRF)

anim_time_rate = 3 # skip some frames for a fast-forward effect

# Animate ICRF

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set(xlabel='X (km)', ylabel='Y (km)', zlabel='Z (km)')
for ax_ in (ax.xaxis, ax.yaxis, ax.zaxis):
    ax_.set_major_formatter(tck.EngFormatter())
plot_ephem(orion_ICRF, label='Artemis 1')
plot_ephem(earth_ICRF, color='green', label='Earth')
plot_ephem(moon_ICRF, color='gray', label='Moon')

orion_pos = orion_ICRF[0]["pos"]
earth_pos = earth_ICRF[0]["pos"]
moon_pos = moon_ICRF[0]["pos"]
data = tuple(zip(orion_pos, earth_pos, moon_pos))
points = ax.scatter(*data, c=('blue', 'green', 'gray'), depthshade=False)
time_text = ax.text(-4e5, -4e5, 0, '')

def update_anim(i):
    i = i*anim_time_rate
    orion_pos = orion_ICRF[i]["pos"]
    earth_pos = earth_ICRF[i]["pos"]
    moon_pos = moon_ICRF[i]["pos"]
    data = tuple(zip(orion_pos, earth_pos, moon_pos))
    # points.set_offsets(data)
    points._offsets3d = data # bruh
    points.stale = True
    time_text.set_text(f'Time: {orion_ICRF[i]["date_str"]}')
    return (points, time_text)

ax.set_aspect('equal')
ax.legend()

anim_ = anim.FuncAnimation(fig, update_anim, int(len(orion_ICRF['jd'])/anim_time_rate), interval=10, blit=False)

# Animate ROT

fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
ax2.set(xlabel='X (km)', ylabel='Y (km)', zlabel='Z (km)')
for ax_ in (ax2.xaxis, ax2.yaxis, ax2.zaxis):
    ax_.set_major_formatter(tck.EngFormatter())
ax2.plot(*orion_ROT['pos'].T, label='Artemis 1')
ax2.plot(*earth_ROT['pos'].T, color='green', label='Earth')
ax2.plot(*moon_ROT['pos'].T, color='gray', label='Moon')

orion_pos2 = orion_ROT[0]["pos"]
earth_pos2 = earth_ROT[0]["pos"]
moon_pos2 = moon_ROT[0]["pos"]
data2 = tuple(zip(orion_pos2, earth_pos2, moon_pos2))
points2 = ax2.scatter(*data, c=('blue', 'green', 'gray'), depthshade=False)
time_text2 = ax2.text(-4e5, -4e5, -4e4, '')

def update_anim2(i):
    i = i*anim_time_rate
    orion_pos = orion_ROT[i]["pos"]
    earth_pos = earth_ROT[i]["pos"]
    moon_pos = moon_ROT[i]["pos"]
    data2 = tuple(zip(orion_pos, earth_pos, moon_pos))
    # points.set_offsets(data)
    points2._offsets3d = data2 # bruh
    points2.stale = True
    time_text2.set_text(f'Time: {orion_ICRF[i]["date_str"]}')
    return (points2, time_text2)

ax2.set_aspect('equal')
ax2.legend()

anim2 = anim.FuncAnimation(fig2, update_anim2, int(len(orion_ROT['jd'])/anim_time_rate), interval=10, blit=False)

plt.show()