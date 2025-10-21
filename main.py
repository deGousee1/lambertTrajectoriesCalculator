import astroquery
import scipy
import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import EarthLocation
from astroquery.jplhorizons import Horizons
from scipy.optimize import newton
import pandas as pd
from lambert import get_ToF_estimate, get_Corrected_ToF_estimate, get_LambertV, get_vInfinity, get_orbSpeed, \
    get_Peri_Speed

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
from utils import get_julian_date, get_planet_id, get_Clear_ToF_Time, debug_vectors

print("Library initialization completed successfully")

AU = const.au.value
DAY = 86400

from ephemerides import get_planet_vectors
date = input("Date of departure (yyyy-mm-dd): ")
#date = "2025-06-15"
date_julian = get_julian_date(date)

planetName = input("First planet name: ")
#planetName = "Earth"
planet1name = planetName
planet1id = get_planet_id(planetName)

planetName = input("Second planet name: ")
#planetName = "Jupiter"
planet2name = planetName
planet2id = get_planet_id(planetName)

departOrbitHeight = float(input("Departure Orbit Height (km): "))*1000
arrivalOrbitHeight = float(input("Arrival Orbit Height (km): "))*1000

#departOrbitHeight = 200000
#arrivalOrbitHeight = 2000000

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
r1r2angle = np.arccos(np.dot(r_first, r_second) / (r1_norm * r2_norm))

secondsToF = get_ToF_estimate(planet2name, r1_norm)
ToFdays = secondsToF/86400
JulianArrivalETA = ToFdays + date_julian
correctedToF = get_Corrected_ToF_estimate(date_julian, JulianArrivalETA, planet1id, planet2id)
#correctedToF = secondsToF eliminuje drugi stopień obliczeń
daysToF, hoursToF, minutesToF, secsToF = get_Clear_ToF_Time(correctedToF)
correctedToFdays = correctedToF/86400
JulianArrivalCorrected = correctedToFdays + date_julian
v1, v2 = get_LambertV(JulianArrivalCorrected, correctedToFdays, date_julian, planet1id, planet2id, correctedToF)
v1_norm = np.linalg.norm(v1)
v2_norm = np.linalg.norm(v2)

v_arrivalSecond = get_planet_vectors(planet2id, JulianArrivalCorrected)
v_arrivalSecond = np.array(v_arrivalSecond[["vx", "vy", "vz"]].iloc[0]) * AU / DAY

planetVector = v_first
shipVector = v1
vInfinityDepart = get_vInfinity(planetVector, shipVector)
planetVector = v_arrivalSecond
shipVector = v2
vInfinityArrival = get_vInfinity(planetVector, shipVector)

planetName = planet1name
orbitHeight = float(departOrbitHeight)
departOrbitSpeed = get_orbSpeed(orbitHeight, planetName)
planetName = planet2name
orbitHeight = float(arrivalOrbitHeight)
arrivalOrbitSpeed = get_orbSpeed(orbitHeight, planetName)

planetName = planet1name
orbitHeight = float(departOrbitHeight)
vInfinity = vInfinityDepart
departPeriSpeed = get_Peri_Speed(orbitHeight, planetName, vInfinity)
planetName = planet2name
orbitHeight = float(arrivalOrbitHeight)
vInfinity = vInfinityArrival
arrivalPeriSpeed = get_Peri_Speed(orbitHeight, planetName, vInfinity)

departDeltaV = (departPeriSpeed - departOrbitSpeed)
arrivalDeltaV = (arrivalPeriSpeed - arrivalOrbitSpeed)
debug_vectors(v_first, v1, r_first, "departure")
debug_vectors(v_arrivalSecond, v2, r_second, "arrival")

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
print("Estimated time of flight:", daysToF, "days", hoursToF, "hours", minutesToF, "minutes", round(secsToF), "seconds")
print("V1:", v1, "V2:", v2)
print("V1 norm:", v1_norm, "V2 norm:", v2_norm)
print("Departure burn vector (m/s):", departDeltaV, ". Arrival capture burn vector (m/s):", arrivalDeltaV)
print("Delta V needed for transfer from", planet1name, "orbit at height of", departOrbitHeight, "meters:", np.round(np.linalg.norm(departDeltaV), 1), "m/s")
print("Delta V needed for capture at", planet2name, "orbit at", arrivalOrbitHeight, "meters:", np.round(np.linalg.norm(arrivalDeltaV), 1) ,"m/s")
#print("Departure peri speed:", departOrbitSpeed, "Arrival peri speed:", arrivalOrbitSpeed)
#print("V infinity departure:", vInfinityDepart, "V infinity arrival:", vInfinityArrival)
#print(v1_norm - (total_v_first + departOrbitSpeed), "m/s transfer burn")
#print(v2_norm - (total_v_second + arrivalOrbitSpeed), "m/s capture burn")