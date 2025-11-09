import numpy as np
import spiceypy as spice
from db import get_planet_orbPeriod
from lambert import get_ToF_estimate, get_Corrected_ToF_estimate, get_LambertV, get_Optimal_Launch_Angle, get_Delta_V
from plots import transfer_Angle_Scan, porkchop_plot
from utils import get_Clear_ToF_Time, julian_to_utc, ask_for_Entry_Data
from ephemerides import get_spice_planet_vectors
DAY = 86400
#Wczytanie efemerydów
spice.furnsh("kernels/de440.bsp")
spice.furnsh("kernels/naif0012.tls")
spice.furnsh("kernels/pck00010.tpc")
#Zdobycie danych wejściowych
date_julian, planet1name, planet1id, planet2name, planet2id, departOrbitHeight, arrivalOrbitHeight = ask_for_Entry_Data()

first_v = get_spice_planet_vectors(planet1id, date_julian)
second_v = get_spice_planet_vectors(planet2id, date_julian)

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

#określenie czy lecimy z planety wewnętrznej na zewnętrzą czy na odwrót
p1T = get_planet_orbPeriod(planet1name)
p2T = get_planet_orbPeriod(planet2name)
if p1T < p2T:
    outward=True
else:
    outward=False

#obliczenie idealnego heliocentrycznego kąta między planetami do transferu
optimalAngle = abs(get_Optimal_Launch_Angle(planet2name, correctedToFdays, outward))
#obliczenie czasu synodycznego planet
Tsyn=1/abs((1/get_planet_orbPeriod(planet1name) - 1/get_planet_orbPeriod(planet2name)))

#znalezienie przybliżonych okien startowych i idealnych katów do transferu
bestAngleDateJ, bestAngle, utcTransferWindow, worstAngle, utcWorstAngleDate = transfer_Angle_Scan(Tsyn, date_julian, planet1id, planet2id, planet2name, correctedToFdays, outward)
print("Best transfer window date Julian:", bestAngleDateJ)
print("Best transfer window date:", utcTransferWindow)
print("Best transfer window angle:", bestAngle)
print("Worst angle date:", utcWorstAngleDate)
print("Worst angle :", worstAngle)

#określenie zakresu skanu i dokładności dla pierwszego porkchop plota
scanRange=Tsyn*0.1
if scanRange < 100:
    scanRange = 100
scanStep=round(scanRange*0.05)
if scanStep < 1:
    scanStep = 1
print("Scan range:", scanRange, "Scan step:", scanStep)
scanStepToF = 0.01*correctedToFdays
#pierwszy szeroki porkchop plot
porkchopNumber = 1
jd, best_tof, best_deltaV = porkchop_plot(scanRange, scanStep, scanStepToF,bestAngleDateJ, correctedToFdays, planet1id, planet2id, planet1name, planet2name, departOrbitHeight, arrivalOrbitHeight, porkchopNumber)

first_v2 = get_spice_planet_vectors(planet1id, jd) #Spice implemented
r_first2 = np.array(first_v2[["x", "y", "z"]].iloc[0]) * 1000
v_first2 = np.array(first_v2[["vx", "vy", "vz"]].iloc[0]) * 1000

#określenie zakresu skanu i dokładności dla drugiego wykresu
scanRange=20#round(scanRange*0.2)
scanStep=1
#drugi dokładniejszy porkchop plot
porkchopNumber = 2
jd, best_tof, best_deltaV = porkchop_plot(scanRange, scanStep, scanStep, jd, best_tof, planet1id, planet2id, planet1name, planet2name, departOrbitHeight, arrivalOrbitHeight, porkchopNumber)
utcBestLaunch = julian_to_utc(jd)
print("Best time of flight:", best_tof, "Best launch date:",utcBestLaunch, "Best deltaV possible:", best_deltaV)

first_v3 = get_spice_planet_vectors(planet1id, jd)
r_first3 = np.array(first_v3[["x", "y", "z"]].iloc[0]) * 1000
v_first3 = np.array(first_v3[["vx", "vy", "vz"]].iloc[0]) * 1000

second_v = get_spice_planet_vectors(planet2id, jd)
r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
v_second = np.array(second_v[["vx", "vy", "vz"]].iloc[0]) * 1000

JulianArrivalBest: float = float(best_tof) + jd
best_tof = float(best_tof)
print("Best ToF:", best_tof)
print("Julian Arrival Best:", JulianArrivalBest, "Jd (should be departure date)", jd, "ToF:", best_tof*86400)
v1, v2 = get_LambertV(JulianArrivalCorrected=JulianArrivalBest, date_julian=jd, planet1id=planet1id, planet2id=planet2id, correctedToF=best_tof * 86400)
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

v_arrivalSecond = get_spice_planet_vectors(planet_id=planet2id, date=JulianArrivalBest) #Spice implemented
v_arrivalSecond = np.array(v_arrivalSecond[["vx", "vy", "vz"]].iloc[0]) * 1000


arrivalDate = julian_to_utc(JulianArrivalBest)

angleArr = np.degrees(np.arccos(np.dot(v2, v_arrivalSecond) / (np.linalg.norm(v2)*np.linalg.norm(v_arrivalSecond))))

daysToF, hoursToF, minutesToF, secsToF = get_Clear_ToF_Time(correctedToF=best_tof*86400)

r1r2angle = np.arccos(np.dot(r_first3, r_second) / (r1_norm * r2_norm))

print("Angle between", planet1name, "and", planet2name, "relative to the Sun at departure:", np.round(np.degrees(r1r2angle), 2), "°")
print("Optimal launch angle:", np.round((optimalAngle), 2))
print("UTC departure date:", utcBestLaunch, "Julian departure date:", jd)
print("Time of flight:", daysToF, "days", hoursToF, "hours", minutesToF, "minutes", round(secsToF), "seconds")
print("UTC arrival date:", arrivalDate, "Julian arrival date:", JulianArrivalBest)
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