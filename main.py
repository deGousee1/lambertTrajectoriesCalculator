import astroquery
import scipy
import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import EarthLocation
from astroquery.jplhorizons import Horizons
from scipy.optimize import newton
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)

from utils import get_julian_date, get_planet_id

print("Library initialization completed successfully")

AU = const.au.value
DAY = 86400


from ephemerides import get_planet_vectors
date = "2026-12-01"
date_julian = get_julian_date(date)

planetname = input("First planet name: ")
planet1name = planetname
planet1id = get_planet_id(planetname)

planetname = input("Second planet name: ")
planet2name = planetname
planet2id = get_planet_id(planetname)

first_v = get_planet_vectors(planet1id, date_julian)
second_v = get_planet_vectors(planet2id, date_julian)

r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * AU
r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * AU

v_first = np.array(first_v[["vx", "vy", "vz"]].iloc[0]) * AU/DAY
v_second = np.array(second_v[["vx", "vy", "vz"]].iloc[0]) * AU / DAY

print(first_v)
print(second_v)
print(planet1name, "position in meters:", r_first)
print(planet2name, "position in meters:", r_second)
print(planet1name, "velocity in m/s:", v_first)
print(planet2name, "velocity in m/s:", v_second)
print("UTC date:", date)
print("Julian date:", date_julian)