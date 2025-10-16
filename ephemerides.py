from astroquery.jplhorizons import Horizons
import pandas as pd

def get_planet_vectors(planet_id, date):
    obj = Horizons(id=planet_id, location ='@0', epochs=[date])
    vectors = obj.vectors()
    df = vectors[['datetime_str', 'x', 'y', 'z', 'vx', 'vy', 'vz']].to_pandas()
    return df