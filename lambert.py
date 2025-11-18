import numpy as np
from db import get_planet_semimajor, get_planet_GM, get_planet_Radius, get_planet_orbPeriod
from ephemerides import get_spice_planet_vectors
from utils import stumpff_C, stumpff_S

def get_ToF_estimate(planet2name, r1_norm):
    # Funkcja oblicza wstępny czas lotu na podstawie faktycznej odległości planety startowej od słońca w chwili startu i półosi wielkiej orbity planety docelowej
    planetName = planet2name
    planet2semimajor = get_planet_semimajor(planetName)

    planetName = "Sun"
    sunGM = get_planet_GM(planetName)

    estToF = np.pi * np.sqrt( ( ( (r1_norm + planet2semimajor) / 2 ) ** 3 ) / sunGM)
    return estToF

def get_Corrected_ToF_estimate(date_julian, JulianArrivalETA, planet1id, planet2id):
    # Funkcja oblicza w teorii dokładniejsze przybliżenie czasu lotu na podstawie odległości planety startowej od słońca w chwili startu i odległości planety docelowej od słońca w chwili przylotu (którą mogę wyznaczyć na podstawie poprzedniego, mniej dokładnego przybliżenia)
    first_v = get_spice_planet_vectors(planet1id, date_julian)
    second_v = get_spice_planet_vectors(planet2id, JulianArrivalETA)
    r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000
    r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
    r1_norm = np.linalg.norm(r_first)
    r2_norm = np.linalg.norm(r_second)
    planetName = "Sun"
    sunGM = get_planet_GM(planetName)
    correctedToF = np.pi * np.sqrt( ( ( (r1_norm + r2_norm) / 2 ) ** 3 ) / sunGM)
    return correctedToF

def get_Optimal_Launch_Angle(planet2name, correctedToFdays, outward):
    # Funkcja oblicza optymalny kąt heliocentryczny do transferu
    planetName = planet2name
    planet2orbPeriod = float(get_planet_orbPeriod(planetName))
    omega=360/planet2orbPeriod
    deltaPhi=omega*correctedToFdays
    if outward:
        optimalAngle = 180 - deltaPhi
    else:
        optimalAngle = deltaPhi - 180
    return optimalAngle % 180

def get_LambertV(JulianArrivalCorrected, date_julian, planet1id, planet2id, correctedToF):
    #Zdobycie wektorów
    origin_vec = get_spice_planet_vectors(planet1id, date_julian)
    dest_vec = get_spice_planet_vectors(planet2id, JulianArrivalCorrected)
    pos_origin = np.array(origin_vec[["x", "y", "z"]].iloc[0]) * 1000
    pos_dest = np.array(dest_vec[["x", "y", "z"]].iloc[0]) * 1000
    #Wektor prędkosci planety docelowej dla sprawdzenia rozwiązania long way
    v_dest = np.array(dest_vec[["vx", "vy", "vz"]].iloc[0]) * 1000
    r1_norm = np.linalg.norm(pos_origin)
    r2_norm = np.linalg.norm(pos_dest)
    #Zdobycie parametru grawitacyjnego słońca
    planetName = "Sun"
    sunGM = get_planet_GM(planetName)
    #Obliczenie kąta między wektorami
    theta = np.arccos(np.dot(pos_origin, pos_dest) / (r1_norm * r2_norm))
    A = np.sin(theta) * np.sqrt( (r1_norm * r2_norm) / (1 - np.cos(theta)) )
    #Inicjalizacja zmiennych
    ToFLambertMid = 1
    yzMid = 2
    IterationCounter=0
    v1 = 0
    v2 = 0
    short_way_done = False
    #Problem Lamberta
    while not short_way_done:
        z1 = 0
        z2 = 20
        while abs(correctedToF - ToFLambertMid) > 0.001:
            Cz1 = stumpff_C(z1)
            Sz1 = stumpff_S(z1)
            Cz2 = stumpff_C(z2)
            Sz2 = stumpff_S(z2)

            yz1 = r1_norm + r2_norm + A * ((z1 * Sz1 - 1) / np.sqrt(Cz1))
            yz2 = r1_norm + r2_norm + A * ((z2 * Sz2 - 1) / np.sqrt(Cz2))

            ToFLambert1 = (((yz1 / Cz1) ** 1.5) * Sz1 + A * np.sqrt(yz1)) / np.sqrt(sunGM)
            ToFLambert2 = (((yz2 / Cz2) ** 1.5) * Sz2 + A * np.sqrt(yz2)) / np.sqrt(sunGM)
            iloczyn = (correctedToF - ToFLambert1) * (correctedToF - ToFLambert2)
            if iloczyn < 0:
                zMid = (z1+z2)/2
                CzMid = stumpff_C(zMid)
                SzMid = stumpff_S(zMid)
                yzMid = r1_norm + r2_norm + A * ((zMid * SzMid - 1) / np.sqrt(CzMid))
                ToFLambertMid = (((yzMid / CzMid) ** 1.5) * SzMid + A * np.sqrt(yzMid)) / np.sqrt(sunGM)
                iloczyn = (correctedToF - ToFLambert1) * (correctedToF - ToFLambertMid)
                if iloczyn > 0:
                    z1 = zMid
                else:
                    z2 = zMid
            #else:
                #print("1:", (correctedToF - ToFLambert1), "2:", (correctedToF - ToFLambert2))
                #print("z1:", z1, "z2:", z2) Zakomentowany kod przydawał się wcześniej do debugowania

            IterationCounter += 1
            # Wykrywanie zbyt dużej liczby iteracji (po zastosowaniu metody bisekcji nie powinno nigdy się aktywować, bo program znajduje parametr z w 30-35 iteracji
            if IterationCounter>1000:
                print("Error iteration count too high")


        #Wyznaczenie parametrów f, g i g'
        f=1 - (yzMid / r1_norm)
        g = A * np.sqrt(yzMid/sunGM)
        gp = 1-(yzMid/r2_norm)
        v1=(1/g)*(pos_dest-(f*pos_origin))
        v2=(1/g)*(gp*pos_dest-pos_origin)
        # Sprawdzenie czy podana trajektoria jest krótszą drogą
        angle = np.arccos(np.dot(v2, v_dest) / (np.linalg.norm(v2)*np.linalg.norm(v_dest)))
        if angle > np.pi/2:
            A = -A
            continue
        else:
            short_way_done = True
    return v1, v2
def get_vInfinity(planetVector, shipVector):
    vInfinity = shipVector - planetVector
    return vInfinity

def get_orbSpeed(orbitHeight, planetName):
    planetRadius = get_planet_Radius(planetName)
    planetGM = get_planet_GM(planetName)
    orbitSpeed = np.sqrt(planetGM/(orbitHeight+planetRadius))
    return orbitSpeed

def get_Peri_Speed(orbitHeight, planetName, vInfinity):
    planetRadius = get_planet_Radius(planetName)
    planetGM = get_planet_GM(planetName)
    periSpeed = np.sqrt(np.linalg.norm(vInfinity)**2 + 2*planetGM/(planetRadius+orbitHeight))
    return periSpeed

def get_Delta_V(planet2id, planet1id, JulianArrivalCorrected, v1, v2, planet1name, planet2name, departOrbitHeight, arrivalOrbitHeight, jd):
    v_arrivalSecond = get_spice_planet_vectors(planet2id, JulianArrivalCorrected)
    v_arrivalSecond = np.array(v_arrivalSecond[["vx", "vy", "vz"]].iloc[0]) * 1000


    first_v = get_spice_planet_vectors(planet1id, jd)
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

    departDeltaV = (departPeriSpeed - departOrbitSpeed)
    arrivalDeltaV = (arrivalPeriSpeed - arrivalOrbitSpeed)
    return departDeltaV, arrivalDeltaV