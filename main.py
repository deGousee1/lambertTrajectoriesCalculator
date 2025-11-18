import numpy as np
import spiceypy as spice
from db import get_planet_orbPeriod
from lambert import get_ToF_estimate, get_Corrected_ToF_estimate, get_LambertV, get_Optimal_Launch_Angle, get_Delta_V
from plots import transfer_Angle_Scan, porkchop_plot, create_Result_PDF_File
from utils import get_Clear_ToF_Time, julian_to_utc, ask_for_Entry_Data, welcomeScreenprint
from ephemerides import get_spice_planet_vectors
import os
import sys

# Sprawdzenie czy program jest odpalony przez .exe czy przez IDE
if getattr(sys, 'frozen', False):
    base_dir = sys._MEIPASS  # base_dir kiedy program jest odpalony przez .exe
else:
    base_dir = os.path.dirname(os.path.abspath(__file__))  # base_dir kiedy program jest odpalony w IDE

# Ustalenie ścieżek do folderów wewnętrznych
working_dir = os.path.join(base_dir, 'workingFiles')
results_dir = os.path.join(base_dir, 'results')
kernels_dir = os.path.join(base_dir, "kernels")

# Utworzenie folderów jeśli ich nie ma
os.makedirs(working_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)
os.makedirs(kernels_dir, exist_ok=True)

# Załadowanie plików z bazy danych SPICE
spice.furnsh(os.path.join(kernels_dir, "de440.bsp"))
spice.furnsh(os.path.join(kernels_dir, "naif0012.tls"))
spice.furnsh(os.path.join(kernels_dir, "pck00010.tpc"))

# Wydrukowanie intra
welcomeScreenprint()

#Zdobycie danych wejściowych
date_julian, planet1name, planet1id, planet2name, planet2id, departOrbitHeight, arrivalOrbitHeight = ask_for_Entry_Data()

# Pobranie wektorów
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

secondsToF = get_ToF_estimate(planet2name, r1_norm) # Pierwsze obliczenie przybliżonego czasu lotu
ToFdays = secondsToF/86400
JulianArrivalETA = ToFdays + date_julian
correctedToF = get_Corrected_ToF_estimate(date_julian, JulianArrivalETA, planet1id, planet2id) # Drugie, w teorii dokładniejsze obliczenie czasu lotu
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

print("-------------------------------------------")
print("Calculating the given transfer maneuver... Please wait...")
print("-------------------------------------------")
#znalezienie przybliżonych okien startowych i idealnych katów do transferu
bestAngleDateJ, bestAngle, utcTransferWindow, worstAngle, utcWorstAngleDate = transfer_Angle_Scan(Tsyn, date_julian, planet1id, planet2id, planet2name, correctedToFdays, outward)

#określenie zakresu skanu i dokładności dla pierwszego porkchop plota
scanRange=Tsyn*0.1
if scanRange < 100:
    scanRange = 100
scanStep=round(scanRange*0.05)
if scanStep < 1:
    scanStep = 1
scanStepToF = 0.01*correctedToFdays
#pierwszy szeroki porkchop plot
porkchopNumber = 1
jd, best_tof, best_deltaV = porkchop_plot(scanRange, scanStep, scanStepToF,bestAngleDateJ, correctedToFdays, planet1id, planet2id, planet1name, planet2name, departOrbitHeight, arrivalOrbitHeight, porkchopNumber)

# Określenie zakresu skanu i dokładności dla drugiego wykresu
scanRange=20#round(scanRange*0.2)
scanStep=1
# Drugi dokładniejszy porkchop plot
porkchopNumber = 2
jd, best_tof, best_deltaV = porkchop_plot(scanRange, scanStep, scanStep, jd, best_tof, planet1id, planet2id, planet1name, planet2name, departOrbitHeight, arrivalOrbitHeight, porkchopNumber)
utcBestLaunch = julian_to_utc(jd)

first_v3 = get_spice_planet_vectors(planet1id, jd)
r_first3 = np.array(first_v3[["x", "y", "z"]].iloc[0]) * 1000

second_v = get_spice_planet_vectors(planet2id, jd)
r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000

best_tof = float(best_tof)
JulianArrivalBest = best_tof + jd
# Ostatni Lambert dla najlepszych parametrów daty startu i czasu lotu według porkchop plotów
v1, v2 = get_LambertV(JulianArrivalCorrected=JulianArrivalBest, date_julian=jd, planet1id=planet1id, planet2id=planet2id, correctedToF=best_tof * 86400)
v1_norm = np.linalg.norm(v1)
v2_norm = np.linalg.norm(v2)

# Obliczenie delta V potrzebnego do manewru
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
deltaV = departDeltaV + arrivalDeltaV

v_arrivalSecond = get_spice_planet_vectors(planet_id=planet2id, date=JulianArrivalBest)
v_arrivalSecond = np.array(v_arrivalSecond[["vx", "vy", "vz"]].iloc[0]) * 1000

arrivalDate = julian_to_utc(JulianArrivalBest)

angleArr = np.degrees(np.arccos(np.dot(v2, v_arrivalSecond) / (np.linalg.norm(v2)*np.linalg.norm(v_arrivalSecond))))

daysToF, hoursToF, minutesToF, secsToF = get_Clear_ToF_Time(correctedToF=best_tof*86400)

r1r2angle = np.arccos(np.dot(r_first3, r_second) / (r1_norm * r2_norm))

print("-------------------------------------------")
print("Transfer maneuver calculations completed!")
print("-------------------------------------------")
print("Maneuver data:")
#print("Angle between", planet1name, "and", planet2name, "relative to the Sun at departure:", np.round(np.degrees(r1r2angle), 2), "°")
#print("Optimal launch angle:", np.round((optimalAngle), 2))
print("UTC departure date:", utcBestLaunch, "Julian departure date:", jd)
print("Time of flight:", daysToF, "days", hoursToF, "hours", minutesToF, "minutes", round(secsToF), "seconds")
print("UTC arrival date:", arrivalDate, "Julian arrival date:", JulianArrivalBest)
print("Delta V needed for transfer from", planet1name, "orbit at height of", departOrbitHeight/1000, "km:", np.round(departDeltaV, 1), "m/s")
print("Delta V needed for capture at", planet2name, "orbit at", arrivalOrbitHeight/1000, "km:", np.round(arrivalDeltaV, 1) ,"m/s")
print(f"Ship's arrival velocity vector angle relative to the {planet2name} velocity vector at arrival: {np.round(angleArr, 1)}°\n")
#if angleArr>25:
#    print("WARNING! Catastrophically inefficient trajectory! The launch date is most likely set far from the Hohmann transfer window!")
#elif angleArr>15:
#    print("Warning! Inefficient trajectory! The launch date may not be set during the perfect transfer window.")
#elif angleArr>5:
#    print("Warning! Slightly inefficient trajectory! The launch date may not be set during the perfect transfer window.")
#else:
#    print("Good Hohmann transfer found!")
# Zakomentowany kod powyżej był moim pomysłem na debugowanie kodu. Przy idealnym transferze Hohmanna kąt między wektorem prędkosci statku a planety docelowej powinien być mały.
# Jeśli program znajdował złą datę transferu i "wymuszał" manewr to trajektoria nie przypominała idealnego Hohmanna i ten kod miał to wykrywać.
# Obecnie program nie spełnia swojego zadania bo przy transferach np. z Ziemi na Neptuna kąt i tak jest dość duży i program fałszywie zwraca błąd a poza tym mam inne narzędzia do debugowania jak porkchop ploty czy ten graf kątów między planetami

manData_path = os.path.join(working_dir, 'manData.txt')
# Utworzenie opisu manewru .txt
with open(manData_path, 'w') as f:
    f.write('Maneuver data:\n')
    f.write(f"UTC departure date: {utcBestLaunch} Julian departure date: {np.round(jd,2)}\n")
    f.write(f"Time of flight: {daysToF} days {hoursToF} hours {minutesToF} minutes\n")
    f.write(f"UTC arrival date: {arrivalDate} Julian arrival date: {np.round(JulianArrivalBest, 2)}\n")
    f.write(f"Delta V needed for transfer from {planet1name} orbit at height of {departOrbitHeight/1000 }km to {planet2name}: {np.round(departDeltaV, 1)} m/s\n")
    f.write(f"Delta V needed for capture at {planet2name} for an orbit height of {arrivalOrbitHeight/1000 }km: {np.round(arrivalDeltaV, 1)} m/s\n")
    f.write(f"Ship's arrival velocity vector angle relative to the {planet2name} velocity vector at arrival: {np.round(angleArr, 1)}°\n")
tWindowData_path = os.path.join(working_dir, 'tWindowData.txt')
# Utworzenie drugiego opisu w .txt
with open(tWindowData_path, 'w') as f:
    f.write('Transfer window data:\n')
    f.write(f"Transfer window date: {utcTransferWindow}\n")
    f.write(f"Optimal transfer angle: {np.round(optimalAngle, 1)}°\n")
    f.write(f"Worst transfer angle date: {utcWorstAngleDate}\n")
    f.write(f"Synodic time for a transfer from {planet1name} to {planet2name}: {np.round(Tsyn, 1)} days\n")
create_Result_PDF_File(planet1name, planet2name)

# Zapytanie czy otworzyć plik PDF z wynikami
openOrNot = input("Do you want to open a detailed PDF result file? (y/n)")
if openOrNot == 'y':
    results_path = os.path.join(results_dir, 'ResultFile.pdf')
    os.startfile(results_path)
    sys.exit()
elif openOrNot == 'n':
    sys.exit()
else:
    print(". . .")
    sys.exit()