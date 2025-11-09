from astroquery.jplhorizons import Horizons
import spiceypy as spice
import pandas as pd

def get_spice_planet_vectors(planet_id: int, date: float):
    et = (date - 2451545.0) * 86400.0

    target = str(planet_id)
    observer = "0"  # barycentrum
    ref_frame = "ECLIPJ2000"
    abcorr = "NONE"

    state, lt = spice.spkezr(target, et, ref_frame, abcorr, observer)

    df = pd.DataFrame([{
        "x": state[0],
        "y": state[1],
        "z": state[2],
        "vx": state[3],
        "vy": state[4],
        "vz": state[5]
    }])
    return df

def get_planet_vectors(planet_id, date):
    obj = Horizons(id=planet_id, location ='@0', epochs=[date])
    vectors = obj.vectors()
    df = vectors[['datetime_str', 'x', 'y', 'z', 'vx', 'vy', 'vz']].to_pandas()
    return df