import csv
import matplotlib.pyplot as plt
import matplotlib.animation as anim

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
                    entry = {
                        "JD": float(jd),
                        "date_str": date_str,
                        "pos": (float(x), float(y), float(z)),
                        "vel": (float(v_x), float(v_y), float(v_z)),
                        "light": float(lt),
                        "range": float(range),
                        "range_rate": float(range_rate),
                        "extra": extra
                    }
                    res.append(entry)
            else:
                if row[0] == '$$SOE': # start of ephemeris
                    reading_eph = True
                    continue
    
    return res

def at_time(eph, julian_date):
    '''Get entry from the given ephemeris at the given time'''
    # Maybe at some future point we can do a fancier interpolation or propagation
    # For now let's just do nearest point lookup
    # Actually do we even need this function? nvm let's leave it for now
    pass

def plot_ephem(eph, **kwargs):
    posns = [entry["pos"] for entry in eph]
    posns_zipped = tuple(zip(*posns))
    ax.plot(*posns_zipped, **kwargs)

# orion_ICRF = load_eph_file("horizons_results_orion_barycenter_1day_csv.txt")
orion_ICRF = load_eph_file("horizons_results_orion_barycenter_1hour_csv.txt")
earth_ICRF = load_eph_file("horizons_results_earth_barycenter_1hour_csv.txt")
moon_ICRF = load_eph_file("horizons_results_moon_barycenter_1hour_csv.txt")

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set(xlabel='X', ylabel='Y', zlabel='Z')
ax.xaxis.set_major_formatter('{x:.1e}')
ax.yaxis.set_major_formatter('{x:.1e}')
ax.zaxis.set_major_formatter('{x:.1e}')
plot_ephem(orion_ICRF, label='Artemis 1')
plot_ephem(earth_ICRF, color='green', label='Earth')
plot_ephem(moon_ICRF, color='gray', label='Moon')

orion_pos = orion_ICRF[0]["pos"]
earth_pos = earth_ICRF[0]["pos"]
moon_pos = moon_ICRF[0]["pos"]
data = tuple(zip(orion_pos, earth_pos, moon_pos))
points = ax.scatter(*data, c=('blue', 'green', 'gray'))
time_text = ax.text(-4e5, -4e5, -4e4, '')

time_rate = 3 # skip some frames for a fast-forward effect

def update_anim(i):
    i = i*time_rate
    orion_pos = orion_ICRF[i]["pos"]
    earth_pos = earth_ICRF[i]["pos"]
    moon_pos = moon_ICRF[i]["pos"]
    data = tuple(zip(orion_pos, earth_pos, moon_pos))
    # points.set_offsets(data)
    points._offsets3d = data # bruh
    points.stale = True
    time_text.set_text(f'Time: {orion_ICRF[i]["date_str"]}')
    return (points, time_text)

ax.legend()

anim_ = anim.FuncAnimation(fig, update_anim, int(len(orion_ICRF)/time_rate), interval=10, blit=False)

plt.show()