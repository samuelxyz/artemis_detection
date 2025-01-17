import csv
import matplotlib.pyplot as plt

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

orion_ICRF = load_eph_file("horizons_results_orion_barycenter_1day_csv.txt")

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set(xlabel='X', ylabel='Y', zlabel='Z')
ax.xaxis.set_major_formatter('{x:.1e}')
ax.yaxis.set_major_formatter('{x:.1e}')
ax.zaxis.set_major_formatter('{x:.1e}')
posns = [entry["pos"] for entry in orion_ICRF]
posns_zipped = tuple(zip(*posns))
ax.plot(*posns_zipped)

plt.show()