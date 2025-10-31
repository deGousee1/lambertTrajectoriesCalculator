import sys
import numpy as np
#from astropy import constants as const
from astropy.time import Time
import spiceypy as spice
import pandas as pd
import matplotlib.pyplot as plt

from db import get_planet_orbPeriod
from lambert import get_ToF_estimate, get_Corrected_ToF_estimate, get_LambertV, get_vInfinity, get_orbSpeed, \
    get_Peri_Speed, get_Optimal_Launch_Angle, get_Delta_V

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
from utils import get_julian_date, get_planet_id, get_Clear_ToF_Time, julian_to_utc
#AU = const.au.value
DAY = 86400

spice.furnsh("kernels/de440.bsp")
spice.furnsh("kernels/naif0012.tls")
spice.furnsh("kernels/pck00010.tpc")

from ephemerides import get_spice_planet_vectors
date = input("Date of departure (yyyy-mm-dd): ")
date_julian = get_julian_date(date)

planetName = input("First planet name: ")
planet1name = planetName
planet1id = get_planet_id(planetName)

planetName = input("Second planet name: ")
planet2name = planetName
planet2id = get_planet_id(planetName)

departOrbitHeight = float(input("Departure Orbit Height (km): "))*1000
arrivalOrbitHeight = float(input("Arrival Orbit Height (km): "))*1000


first_v = get_spice_planet_vectors(planet1id, date_julian) #Spice implemented
second_v = get_spice_planet_vectors(planet2id, date_julian) #Spice implemented

r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000
r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000

v_first = np.array(first_v[["vx", "vy", "vz"]].iloc[0]) * 1000
v_second = np.array(second_v[["vx", "vy", "vz"]].iloc[0]) * 1000

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

p1T = get_planet_orbPeriod(planet1name)
p2T = get_planet_orbPeriod(planet2name)
if p1T < p2T:
    outward=True
else:
    outward=False

optimalAngle = abs(get_Optimal_Launch_Angle(planet2name, correctedToFdays, outward))
#print("IDEALNY KĄT AAAAAA:", optimalAngle)
Tsyn=1/abs((1/get_planet_orbPeriod(planet1name) - 1/get_planet_orbPeriod(planet2name)))

scanRange=Tsyn*0.5
scanStep=round(Tsyn*0.01)
start_jd_array = np.arange(date_julian-scanRange, date_julian+scanRange, scanStep)
iterationGoal = len(start_jd_array)
utc_dates = [Time(jd, format='jd').to_datetime() for jd in start_jd_array]
angleLoopCounter = 0
angleDifference=0
angle_matrix = np.zeros(len(start_jd_array))
for i, jd_start_val in enumerate(start_jd_array):
    jd_start_val: float = jd_start_val
    first_v = get_spice_planet_vectors(planet1id, jd_start_val)  # Spice implemented
    second_v = get_spice_planet_vectors(planet2id, jd_start_val)  # Spice implemented

    r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000
    r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
    r1_norm = np.linalg.norm(r_first)
    r2_norm = np.linalg.norm(r_second)

    optimalAngle = abs(get_Optimal_Launch_Angle(planet2name, correctedToFdays, outward))
    realAngle = np.degrees(np.arccos(np.dot(r_first, r_second) / (r1_norm * r2_norm)))
    angleDifference = abs(optimalAngle - realAngle)
    angle_matrix[i] = angleDifference
    angleLoopCounter += 1
    sys.stdout.write(
        f"\rTransfer window search progress: {angleLoopCounter} of {round(iterationGoal)}"
    )
    sys.stdout.flush()
print()
plt.plot(utc_dates, angle_matrix, marker='o')
plt.xlabel('Maneuver start date', fontsize=11)
plt.ylabel('Angle difference')
plt.title('Optimal transfer angle deviation plot')
plt.show()
min_idx = np.argmin(angle_matrix)
min_angle = angle_matrix[min_idx]

angleDifference=0
jd_start_val = date_julian-scanRange
jd_start_val: float = jd_start_val
currentBestAngle = 181
worstAngle = 0
bestAngleDateJ = 0
currentBestAngleDateJ = 0
worstAngleDateJ = 0
windowFound = False
optimalAngle = abs(get_Optimal_Launch_Angle(planet2name, correctedToFdays, outward))

for jdDate in start_jd_array: # Finding maximum angle difference
    jdDate: float = jdDate
    first_v = get_spice_planet_vectors(planet1id, jdDate)  # Spice implemented
    second_v = get_spice_planet_vectors(planet2id, jdDate)  # Spice implemented

    r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000
    r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
    v_firstA = np.array(first_v[["vx", "vy", "vz"]].iloc[0]) * 1000
    v_secondA = np.array(second_v[["vx", "vy", "vz"]].iloc[0]) * 1000
    r1_norm = np.linalg.norm(r_first)
    r2_norm = np.linalg.norm(r_second)

    realAngle = np.degrees(np.arccos(np.dot(r_first, r_second) / (r1_norm * r2_norm)))
    angleDifference = abs(optimalAngle - realAngle)
    if optimalAngle < 100 and optimalAngle > 80:
        if angleDifference > worstAngle:
            velocityVectorAngle= np.arccos(np.dot(v_firstA, v_secondA) / (np.linalg.norm(v_firstA)*np.linalg.norm(v_secondA)))
            utcScanDateDebug = julian_to_utc(jdDate)
            utcAngleDateDebug = julian_to_utc(worstAngleDateJ)
            #print("Velocity Vector Angle:", np.degrees(velocityVectorAngle), utcScanDateDebug)
            #print("Worst angle:", worstAngle, worstAngleDateJ)
            #windowFound = True
            if np.degrees(velocityVectorAngle) > 100:
                if angleDifference > worstAngle:
                    worstAngle = angleDifference
                    worstAngleDateJ = jdDate
    else:
        if angleDifference > worstAngle:
            worstAngle = angleDifference
            worstAngleDateJ = jdDate
    angleLoopCounter += 1
    sys.stdout.write(
        f"\rTransfer window search progress: {angleLoopCounter} of {round(iterationGoal)}"
    )
    sys.stdout.flush()
utcWorstAngleDate = julian_to_utc(worstAngleDateJ)
print("Worst angle date:", utcWorstAngleDate)
print("Worst angle :", worstAngle)

start_jd_array = np.arange(worstAngleDateJ, worstAngleDateJ+scanRange*2, scanStep)

for jdDate in start_jd_array:
    if windowFound == False:
        jdDate: float = jdDate
        first_v = get_spice_planet_vectors(planet1id, jdDate)  # Spice implemented
        second_v = get_spice_planet_vectors(planet2id, jdDate)  # Spice implemented

        r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000
        r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
        r1_norm = np.linalg.norm(r_first)
        r2_norm = np.linalg.norm(r_second)
        optimalAngle = abs(get_Optimal_Launch_Angle(planet2name, correctedToFdays, outward))
        realAngle = np.degrees(np.arccos(np.dot(r_first, r_second) / (r1_norm * r2_norm)))
        angleDifference = abs(optimalAngle - realAngle)
        if angleDifference < currentBestAngle:
            currentBestAngle = angleDifference
            currentBestAngleDateJ = jdDate
            if outward == True:
                if currentBestAngle < 5:
                    bestAngle = currentBestAngle
                    bestAngleDateJ = jdDate
                    windowFound = True
    angleLoopCounter += 1
    sys.stdout.write(
        f"\rTransfer window search progress: {angleLoopCounter} of {round(iterationGoal)}"
    )
    sys.stdout.flush()
print()
bestAngleDateJ = currentBestAngleDateJ
bestAngle = currentBestAngle
if outward == True:
    bestAngleDateJ = currentBestAngleDateJ
    bestAngle = currentBestAngle

utcTransferWindow = julian_to_utc(bestAngleDateJ)
print("Best transfer window date Julian:", bestAngleDateJ)
print("Best transfer window date:", utcTransferWindow)
print("Best transfer window angle:", bestAngle)

#v1, v2 = get_LambertV(JulianArrivalCorrected, date_julian, planet1id, planet2id, correctedToF)
#arrivalDeltaV, departDeltaV = get_Delta_V(planet2id
   #                                                    , v1
  #                                                     , v2
  #                                                     , planet1name
  #                                                     , planet2name
  #                                                     , departOrbitHeight
 #                                                      , arrivalOrbitHeight
 #                                                      , v_first
#                                                       , JulianArrivalCorrected)
#maxPorkValue=departOrbitHeight+2000

scanRange=Tsyn*0.1
if scanRange < 100:
    scanRange = 100
scanStep=round(scanRange*0.05)
if scanStep < 1:
    scanStep = 1
print("Scan range:", scanRange, "Scan step:", scanStep)
scanStepToF = 0.01*correctedToFdays
start_jd_array = np.arange(bestAngleDateJ-scanRange, bestAngleDateJ+scanRange, scanStep)
tof_days_array = np.arange(correctedToFdays*0.5, correctedToFdays+scanRange, scanStepToF)
iterationGoal = len(start_jd_array)*len(tof_days_array)
utc_dates = [Time(jd, format='jd').to_datetime() for jd in start_jd_array]
firstLoopCounter = 0
deltaV_matrix = np.zeros((len(tof_days_array), len(start_jd_array)))
k=1
for i, tof in enumerate(tof_days_array):
    for j, jd_start_val in enumerate(start_jd_array):
        tof:float = tof
        jd_arrival = jd_start_val + tof
        v1, v2 = get_LambertV(
            JulianArrivalCorrected=jd_arrival,
            date_julian=jd_start_val,
            planet1id=planet1id,
            planet2id=planet2id,
            correctedToF=tof * 86400,
            k=k
        )
        if (v1 == 100).all():
            v1, v2 = get_LambertV(
                JulianArrivalCorrected=jd_arrival,
                date_julian=jd_start_val,
                planet1id=planet1id,
                planet2id=planet2id,
                correctedToF=tof * 86400,
                k=100
            )
            k=100
        #Obliczenia deltaV
        jd_arrival = jd_start_val + tof
        departDeltaV, arrivalDeltaV = get_Delta_V(planet2id=planet2id,
                                                  planet1id=planet1id
                                                       , v1=v1
                                                       , v2=v2
                                                       , planet1name=planet1name
                                                       , planet2name=planet2name
                                                       , departOrbitHeight=departOrbitHeight
                                                       , arrivalOrbitHeight=arrivalOrbitHeight
                                                       , JulianArrivalCorrected=jd_arrival
                                                       , jd=jd_start_val)
        deltaV = departDeltaV + arrivalDeltaV
        deltaV_matrix[i, j] = deltaV
        firstLoopCounter += 1
        if k==100:
            sys.stdout.write(
                f"\rPorkchop iteration progress: {firstLoopCounter} of {round(iterationGoal)}. Error detected! Lambert will now run in slower backup mode." #Można dodać odniesienie do instrukcji gdzie będzie powiedziane co w tym wypadku zrobić
            )
            sys.stdout.flush()
        else:
            sys.stdout.write(
                f"\rPorkchop iteration progress: {firstLoopCounter} of {round(iterationGoal)}"
            )
            sys.stdout.flush()

print()
print("First rough sieve porkchop graph done!")

plt.figure(figsize=(10,6))
X, Y = np.meshgrid(utc_dates, tof_days_array)

min_idx = np.unravel_index(np.argmin(deltaV_matrix), deltaV_matrix.shape)
i_min, j_min = min_idx
best_tof = tof_days_array[i_min]
jd = float(start_jd_array[j_min])
utcBestLaunch = julian_to_utc(jd)
best_deltaV = deltaV_matrix[i_min, j_min]
print("Best time of flight:", best_tof, "Best launch date:",utcBestLaunch, "Best deltaV possible:", best_deltaV)

deltaV_matrix_masked = np.ma.masked_greater(deltaV_matrix, best_deltaV*2)
plt.contourf(X, Y, deltaV_matrix_masked, levels=50, cmap='viridis')
plt.colorbar(label='Delta-V [m/s]')
plt.xlabel('Data startu')
plt.ylabel('Czas lotu [dni]')
plt.title('Porkchop plot')
plt.show()
#min_idx = np.unravel_index(np.argmin(deltaV_matrix), deltaV_matrix.shape)
#i_min, j_min = min_idx
#best_tof = tof_days_array[i_min]
#jd = float(start_jd_array[j_min])
#utcBestLaunch = julian_to_utc(jd)
#best_deltaV = deltaV_matrix[i_min, j_min]
#print("Best time of flight:", best_tof, "Best launch date:",utcBestLaunch, "Best deltaV possible:", best_deltaV)

first_v2 = get_spice_planet_vectors(planet1id, jd) #Spice implemented
r_first2 = np.array(first_v2[["x", "y", "z"]].iloc[0]) * 1000
v_first2 = np.array(first_v2[["vx", "vy", "vz"]].iloc[0]) * 1000

scanRange=20#round(scanRange*0.2)
scanStep=1
start_jd_array = np.arange(jd-scanRange, jd+scanRange, scanStep)
tof_days_array = np.arange(best_tof-scanRange, best_tof+scanRange, scanStep)
iterationGoal = len(start_jd_array)*len(tof_days_array)
utc_dates = [Time(jd, format='jd').to_datetime() for jd in start_jd_array]
secondLoopCounter = 0
deltaV_matrix = np.zeros((len(tof_days_array), len(start_jd_array)))
k=1
for i, tof in enumerate(tof_days_array):
    for j, jd_start_val in enumerate(start_jd_array):
        jd_arrival = jd_start_val + tof
        v1, v2 = get_LambertV(
            JulianArrivalCorrected=jd_arrival,
            date_julian=jd_start_val,
            planet1id=planet1id,
            planet2id=planet2id,
            correctedToF=tof * 86400,
            k=k
        )
        if (v1 == 100).all():
            v1, v2 = get_LambertV(
                JulianArrivalCorrected=jd_arrival,
                date_julian=jd_start_val,
                planet1id=planet1id,
                planet2id=planet2id,
                correctedToF=tof * 86400,
                k=100
            )
            k = 100

        #Obliczenia deltaV
        jd_arrival = jd_start_val + tof
        departDeltaV, arrivalDeltaV = get_Delta_V(planet2id=planet2id
                                                  , planet1id=planet1id
                                                       , v1=v1
                                                       , v2=v2
                                                       , planet1name=planet1name
                                                       , planet2name=planet2name
                                                       , departOrbitHeight=departOrbitHeight
                                                       , arrivalOrbitHeight=arrivalOrbitHeight
                                                       , JulianArrivalCorrected=jd_arrival
                                                       , jd=jd_start_val)
        deltaV = departDeltaV + arrivalDeltaV
        deltaV_matrix[i, j] = deltaV
        secondLoopCounter += 1
        if k==100:
            sys.stdout.write(
                f"\rPorkchop iteration progress: {secondLoopCounter} of {round(iterationGoal)}. Error detected! Lambert will now run in slower backup mode." #Można dodać odniesienie do instrukcji gdzie będzie powiedziane co w tym wypadku zrobić
            )
            sys.stdout.flush()
        else:
            sys.stdout.write(
                f"\rPorkchop iteration progress: {secondLoopCounter} of {round(iterationGoal)}"
            )
            sys.stdout.flush()

print()
print("Second fine sieve porkchop graph done!")

utc_labels = [d.strftime("%m-%d") for d in utc_dates]
plt.figure(figsize=(10,6))
X, Y = np.meshgrid(utc_dates, tof_days_array)
deltaV_matrix_masked = np.ma.masked_greater(deltaV_matrix, deltaV*2)
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

first_v3 = get_spice_planet_vectors(planet1id, jd) #Spice implemented
r_first3 = np.array(first_v3[["x", "y", "z"]].iloc[0]) * 1000
v_first3 = np.array(first_v3[["vx", "vy", "vz"]].iloc[0]) * 1000

second_v = get_spice_planet_vectors(planet2id, jd) #Spice implemented
r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
v_second = np.array(second_v[["vx", "vy", "vz"]].iloc[0]) * 1000

JulianArrivalBest: float = float(best_tof) + jd
best_tof = float(best_tof)
print("Best ToF:", best_tof)
print("Julian Arrival Best:", JulianArrivalBest, "Jd (should be departure date)", jd, "ToF:", best_tof*86400)
k=1
v1, v2 = get_LambertV(JulianArrivalCorrected=JulianArrivalBest, date_julian=jd, planet1id=planet1id, planet2id=planet2id, correctedToF=best_tof * 86400, k=k)
if (v1 == 100).all():
    v1, v2 = get_LambertV(
        JulianArrivalCorrected=JulianArrivalBest,
        date_julian=jd,
        planet1id=planet1id,
        planet2id=planet2id,
        correctedToF=best_tof * 86400,
        k=100
    )
    sys.stdout.write(
        f"\rError detected! Lambert will now run in slower backup mode."
        # Można dodać odniesienie do instrukcji gdzie będzie powiedziane co w tym wypadku zrobić
    )
    sys.stdout.flush()
    print()
v1_norm = np.linalg.norm(v1)
v2_norm = np.linalg.norm(v2)

departDeltaV, arrivalDeltaV = get_Delta_V(planet2id=planet2id
                                                       , planet1id=planet1id
                                                       , v1=v1
                                                       , v2=v2
                                                       , planet1name=planet1name
                                                       , planet2name=planet2name
                                                       , departOrbitHeight=departOrbitHeight
                                                       , arrivalOrbitHeight=arrivalOrbitHeight
                                                       , JulianArrivalCorrected=JulianArrivalBest
                                                       , jd=jd)
deltaV = departDeltaV# + arrivalDeltaV

#optimalAngle = get_Optimal_Launch_Angle(planet2name=planet2name, correctedToFdays=best_tof, outward=outward)

v_arrivalSecond = get_spice_planet_vectors(planet_id=planet2id, date=JulianArrivalBest) #Spice implemented
v_arrivalSecond = np.array(v_arrivalSecond[["vx", "vy", "vz"]].iloc[0]) * 1000

first_v = get_spice_planet_vectors(planet1id, jd) #Spice implemented

r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000

v_first = np.array(first_v[["vx", "vy", "vz"]].iloc[0]) * 1000

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

#departDeltaV = (departPeriSpeed - departOrbitSpeed)
#arrivalDeltaV = (arrivalPeriSpeed - arrivalOrbitSpeed)

arrivalDate = julian_to_utc(JulianArrivalBest)

angleArr = np.degrees(np.arccos(np.dot(v2, v_arrivalSecond) / (np.linalg.norm(v2)*np.linalg.norm(v_arrivalSecond))))

daysToF, hoursToF, minutesToF, secsToF = get_Clear_ToF_Time(correctedToF=best_tof*86400)

r1r2angle = np.arccos(np.dot(r_first3, r_second) / (r1_norm * r2_norm))

print("Angle between", planet1name, "and", planetName, "relative to the Sun at departure:", np.round(np.degrees(r1r2angle), 2), "°")
print("Optimal launch angle:", np.round((optimalAngle), 2))
print("UTC departure date:", utcBestLaunch, "Julian departure date:", jd)
#print(first_v)
#print(planet1name, "position in meters:", r_first)
#print(planet1name, "velocity in m/s:", v_first)
#print(planet1name, "total velocity in m/s:", total_v_first)
#print(second_v)
#print(planet2name, "position in meters:", r_second)
#print(planet2name, "velocity in m/s:", v_second)
#print(planet2name, "total velocity in m/s:", total_v_second)
print("Time of flight:", daysToF, "days", hoursToF, "hours", minutesToF, "minutes", round(secsToF), "seconds")
print("UTC arrival date:", arrivalDate, "Julian arrival date:", JulianArrivalBest)
#print("V1:", v1, "V2:", v2)
#print("V1 norm:", v1_norm, "V2 norm:", v2_norm)
#("Departure burn vector (m/s):", departDeltaV, ". Arrival capture burn vector (m/s):", arrivalDeltaV)
print("Delta V needed for transfer from", planet1name, "orbit at height of", departOrbitHeight/1000, "km:", np.round(departDeltaV, 1), "m/s")
print("Delta V needed for capture at", planet2name, "orbit at", arrivalOrbitHeight/1000, "km:", np.round(arrivalDeltaV, 1) ,"m/s")
print("Angle at arrival:", angleArr, "°")
if angleArr>25:
    print("WARNING! Catastrophically inefficient trajectory! The launch date is most likely set far from the Hohmann transfer window!")
elif angleArr>15:
    print("Warning! Inefficient trajectory! The launch date may not be set during the perfect transfer window.")
elif angleArr>5:
    print("Warning! Slightly inefficient trajectory! The launch date may not be set during the perfect transfer window.")
else:
    print("Good Hohmann transfer found!")