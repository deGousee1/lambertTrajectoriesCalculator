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
date = input("Date of departure (yyyy-mm-dd): ")
date_julian = get_julian_date(date)

planetName = input("First planet name: ")
planet1name = planetName
planet1id = get_planet_id(planetName)

planetName = input("Second planet name: ")
planet2name = planetName
planet2id = get_planet_id(planetName)

first_v = get_planet_vectors(planet1id, date_julian)
second_v = get_planet_vectors(planet2id, date_julian)

r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * AU
r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * AU

v_first = np.array(first_v[["vx", "vy", "vz"]].iloc[0]) * AU/DAY
v_second = np.array(second_v[["vx", "vy", "vz"]].iloc[0]) * AU / DAY

total_v_first = np.linalg.norm(v_first)
total_v_second = np.linalg.norm(v_second)

r1_norm=np.linalg.norm(r_first)
r2_norm=np.linalg.norm(r_second)
angleEquation = np.dot(r_first, r_second) / (r1_norm * r2_norm)
r1r2angle=np.arccos(angleEquation)

print("UTC date:", date, "Julian date:", date_julian)
print(first_v)
print(planet1name, "position in meters:", r_first)
print(planet1name, "velocity in m/s:", v_first)
print(planet1name, "total velocity in m/s:", total_v_first)
print(second_v)
print(planet2name, "position in meters:", r_second)
print(planet2name, "velocity in m/s:", v_second)
print(planet2name, "total velocity in m/s:", total_v_second)
print("Angle between", planet1name, "and", planetName, "relative to the Sun:", np.round(np.degrees(r1r2angle), 2), "°")