import astroquery
import scipy
import time
from tqdm import tqdm
import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation
from astroquery.jplhorizons import Horizons
from scipy.optimize import newton
import pandas as pd
import matplotlib.pyplot as plt
from lambert import get_ToF_estimate, get_Corrected_ToF_estimate, get_LambertV, get_vInfinity, get_orbSpeed, \
    get_Peri_Speed, get_Optimal_Launch_Angle, get_Delta_V

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
from utils import get_julian_date, get_planet_id, get_Clear_ToF_Time, julian_to_utc

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

scanRange=40
scanStep=4
start_jd_array = np.arange(date_julian-scanRange, date_julian+scanRange, scanStep)
tof_days_array = np.arange(correctedToFdays-scanRange, correctedToFdays+scanRange, scanStep)
utc_dates = [Time(jd, format='jd').to_datetime() for jd in start_jd_array]

rangeA: int=int(((scanRange/scanStep)*2)**2)
progressBar = tqdm(total=int(rangeA), desc="Drawing the porkchop plot")

deltaV_matrix = np.zeros((len(tof_days_array), len(start_jd_array)))
for i in tqdm(range(0, rangeA), desc="Drawing the porkchop plot"):
    for i, tof in enumerate(tof_days_array):
        for j, jd_start_val in enumerate(start_jd_array):
            jd_arrival = jd_start_val + tof
            v1, v2 = get_LambertV(
                JulianArrivalCorrected=jd_arrival,
                correctedToFdays=tof,
                date_julian=jd_start_val,
                planet1id=planet1id,
                planet2id=planet2id,
                correctedToF=tof * 86400
            )
            #Obliczenia deltaV
            jd_arrival = jd_start_val + tof
            arrivalDeltaV, departDeltaV = get_Delta_V(planet2id
                                                           , v1=v1
                                                           , v2=v2
                                                           , planet1name=planet1name
                                                           , planet2name=planet2name
                                                           , departOrbitHeight=departOrbitHeight
                                                           , arrivalOrbitHeight=arrivalOrbitHeight
                                                           , v_first=v_first
                                                           , JulianArrivalCorrected=jd_arrival)
            deltaV = departDeltaV# + arrivalDeltaV
            deltaV_matrix[i, j] = deltaV
            progressBar.update(1)


plt.figure(figsize=(10,6))
X, Y = np.meshgrid(utc_dates, tof_days_array)
deltaV_matrix_masked = np.ma.masked_greater(deltaV_matrix, 20000)
plt.contourf(X, Y, deltaV_matrix_masked, levels=50, cmap='viridis')
plt.colorbar(label='Delta-V [m/s]')
plt.xlabel('Data startu')
plt.ylabel('Czas lotu [dni]')
plt.title('Porkchop plot')
plt.show()
min_idx = np.unravel_index(np.argmin(deltaV_matrix), deltaV_matrix.shape)
i_min, j_min = min_idx
best_tof = tof_days_array[i_min]
jd = float(start_jd_array[j_min])
utcBestLaunch = julian_to_utc(jd)
best_deltaV = deltaV_matrix[i_min, j_min]
print("Best time of flight:", best_tof, "Best launch date:",utcBestLaunch, "Best deltaV possible:", best_deltaV)

scanRange=5
scanStep=1
start_jd_array = np.arange(jd-scanRange, jd+scanRange, scanStep)
tof_days_array = np.arange(best_tof-scanRange, best_tof+scanRange, scanStep)
utc_dates = [Time(jd, format='jd').to_datetime() for jd in start_jd_array]

progressBar2 = tqdm(total=int(scanRange / scanStep), desc="Drawing the porkchop plot")

deltaV_matrix = np.zeros((len(tof_days_array), len(start_jd_array)))
for i in tqdm(range(0, (int(scanRange/scanStep)*2)**2), desc="Drawing the porkchop plot"):
    for i, tof in enumerate(tof_days_array):
        for j, jd_start_val in enumerate(start_jd_array):
            jd_arrival = jd_start_val + tof
            v1, v2 = get_LambertV(
                JulianArrivalCorrected=jd_arrival,
                correctedToFdays=tof,
                date_julian=jd_start_val,
                planet1id=planet1id,
                planet2id=planet2id,
                correctedToF=tof * 86400
            )
            #Obliczenia deltaV
            jd_arrival = jd_start_val + tof
            arrivalDeltaV, departDeltaV = get_Delta_V(planet2id
                                                           , v1=v1
                                                           , v2=v2
                                                           , planet1name=planet1name
                                                           , planet2name=planet2name
                                                           , departOrbitHeight=departOrbitHeight
                                                           , arrivalOrbitHeight=arrivalOrbitHeight
                                                           , v_first=v_first
                                                           , JulianArrivalCorrected=jd_arrival)
            deltaV = departDeltaV# + arrivalDeltaV
            deltaV_matrix[i, j] = deltaV
            progressBar2.update(1)

utc_labels = [d.strftime("%m-%d") for d in utc_dates]
plt.figure(figsize=(10,6))
X, Y = np.meshgrid(utc_dates, tof_days_array)
deltaV_matrix_masked = np.ma.masked_greater(deltaV_matrix, 20000)
plt.contourf(X, Y, deltaV_matrix_masked, levels=50, cmap='viridis')
plt.colorbar(label='Delta-V [m/s]')
plt.xlabel('Data startu')
plt.ylabel('Czas lotu [dni]')
plt.title('Porkchop plot')
plt.show()
min_idx = np.unravel_index(np.argmin(deltaV_matrix), deltaV_matrix.shape)
i_min, j_min = min_idx
best_tof = tof_days_array[i_min]
jd = float(start_jd_array[j_min])
utcBestLaunch = julian_to_utc(jd)
best_deltaV = deltaV_matrix[i_min, j_min]
print("Best time of flight:", best_tof, "Best launch date:",utcBestLaunch, "Best deltaV possible:", best_deltaV)

JulianArrivalBest: float = float(best_tof) + jd

v1, v2 = get_LambertV(JulianArrivalCorrected=JulianArrivalBest, correctedToFdays=best_tof, date_julian=jd, planet1id=planet1id, planet2id=planet2id, correctedToF=best_tof * 86400)
v1_norm = np.linalg.norm(v1)
v2_norm = np.linalg.norm(v2)

optimalAngle = abs(get_Optimal_Launch_Angle(planet2name=planet2name, correctedToFdays=best_tof))

v_arrivalSecond = get_planet_vectors(planet_id=planet2id, date=JulianArrivalBest)
v_arrivalSecond = np.array(v_arrivalSecond[["vx", "vy", "vz"]].iloc[0]) * AU / DAY

first_v = get_planet_vectors(planet1id, jd)

r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * AU

v_first = np.array(first_v[["vx", "vy", "vz"]].iloc[0]) * AU/DAY

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

arrivalDate = julian_to_utc(jd)

angleArr = np.degrees(np.arccos(np.dot(v2, v_arrivalSecond) / (np.linalg.norm(v2)*np.linalg.norm(v_arrivalSecond))))

print("Angle between", planet1name, "and", planetName, "relative to the Sun at departure:", np.round(np.degrees(r1r2angle), 2), "°")
print("Optimal launch angle:", np.round((optimalAngle), 2))
print("UTC departure date:", utcBestLaunch, "Julian departure date:", date_julian)
#print(first_v)
#print(planet1name, "position in meters:", r_first)
#print(planet1name, "velocity in m/s:", v_first)
#print(planet1name, "total velocity in m/s:", total_v_first)
#print(second_v)
#print(planet2name, "position in meters:", r_second)
#print(planet2name, "velocity in m/s:", v_second)
#print(planet2name, "total velocity in m/s:", total_v_second)
print("Time of flight:", daysToF, "days", hoursToF, "hours", minutesToF, "minutes", round(secsToF), "seconds")
print("UTC arrival date:", arrivalDate, "Julian arrival date:", np.round(JulianArrivalBest, 1))
#print("V1:", v1, "V2:", v2)
#print("V1 norm:", v1_norm, "V2 norm:", v2_norm)
#("Departure burn vector (m/s):", departDeltaV, ". Arrival capture burn vector (m/s):", arrivalDeltaV)
print("Delta V needed for transfer from", planet1name, "orbit at height of", departOrbitHeight/1000, "km:", np.round(np.linalg.norm(departDeltaV), 1), "m/s")
print("Delta V needed for capture at", planet2name, "orbit at", arrivalOrbitHeight/1000, "km:", np.round(np.linalg.norm(arrivalDeltaV), 1) ,"m/s")
print("Angle at arrival:", angleArr, "°")
if angleArr>20:
    print("WARNING! Catastrophically inefficient trajectory! The launch date is most likely set far from the Hohmann transfer window!")
elif angleArr>15:
    print("Warning! Inefficient trajectory! The launch date may not be set during the perfect transfer window.")
elif angleArr>5:
    print("Warning! Slightly inefficient trajectory! The launch date may not be set during the perfect transfer window.")
else:
    print("Good Hohmann transfer found!")