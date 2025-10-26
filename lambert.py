import numpy as np
from astropy import constants as const
from scipy.spatial.distance import kulczynski1
from tqdm import tqdm
import sys
from db import get_planet_semimajor, get_planet_GM, get_planet_Radius, get_planet_orbPeriod
from ephemerides import get_spice_planet_vectors
from utils import stumpff_C, stumpff_S

AU = const.au.value
DAY = 86400


def get_ToF_estimate(planet2name, r1_norm):
    planetName = planet2name
    planet2semimajor = get_planet_semimajor(planetName)

    planetName = "Sun"
    sunGM = get_planet_GM(planetName)

    estToF = np.pi * np.sqrt( ( ( (r1_norm + planet2semimajor) / 2 ) ** 3 ) / sunGM)
    return estToF

def get_Corrected_ToF_estimate(date_julian, JulianArrivalETA, planet1id, planet2id):
    first_v = get_spice_planet_vectors(planet1id, date_julian) #Spice implemented
    second_v = get_spice_planet_vectors(planet2id, JulianArrivalETA) #Spice implemented
    r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000
    r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
    r1_norm = np.linalg.norm(r_first)
    r2_norm = np.linalg.norm(r_second)
    planetName = "Sun"
    sunGM = get_planet_GM(planetName)
    correctedToF = np.pi * np.sqrt( ( ( (r1_norm + r2_norm) / 2 ) ** 3 ) / sunGM)
    return correctedToF

def get_Optimal_Launch_Angle(planet2name, correctedToFdays):
    planetName = planet2name
    planet2orbPeriod = float(get_planet_orbPeriod(planetName))
    omega=360/planet2orbPeriod
    deltaPhi=omega*correctedToFdays
    optimalAngle=180-deltaPhi
    return optimalAngle

def get_LambertV(JulianArrivalCorrected, date_julian, planet1id, planet2id, correctedToF):
    #Zdobycie wektorów
    origin_vec = get_spice_planet_vectors(planet1id, date_julian) #Spice implemented
    dest_vec = get_spice_planet_vectors(planet2id, JulianArrivalCorrected) #Spice implemented
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
    #z=9.1526593
    ToFLambert = 1
    yz = 2
    IterationCounter=0
    k1 = 10000
    k2 = 10000
    short_way_done = False
    while not short_way_done:
        z=0
        while abs(correctedToF - ToFLambert) > 0.001:
            #z = z - 0.0000000000001
            #yz = r1_norm+r2_norm(A*(z* stumpff_S(z) -1) / (np.sqrt( stumpff_C(z) )))
            #ToFLambert = ((yz ** 1.5)/stumpff_C(z)) * stumpff_S(z) + A * np.sqrt(yz/sunGM)
            #z -= 1e-12
            if correctedToF - ToFLambert < 0:
                z -=1e-12 * (correctedToF - ToFLambert)*k1
            else:
                z += 1e-12 * (correctedToF - ToFLambert)*k2
            if IterationCounter>1000000:
                #k1 = k1 / 10
                #k2 = k2 / 10
                z=0
            Cz = stumpff_C(z)
            Sz = stumpff_S(z)
            yz = r1_norm + r2_norm + A * ((z * Sz - 1) / np.sqrt(Cz))
            ToFLambert = (((yz / Cz) ** 1.5) * Sz + A * np.sqrt(yz)) / np.sqrt(sunGM) #ToFLambert = ((yz ** 1.5) * Sz / Cz + A * np.sqrt(yz)) / np.sqrt(sunGM)
            IterationCounter += 1
            #sys.stdout.write(
            #    f"\rLambert iteration: {IterationCounter} | ΔToF = {np.round(correctedToF - ToFLambert, 2)} | z = {z}"
            #)
            #sys.stdout.flush()

        f=1 - (yz / r1_norm)
        #g = np.sqrt(yz ** 3 / sunGM)  # zamiast A*sqrt(y/mu)
        g = A * np.sqrt(yz/sunGM)
        gp = 1-(yz/r2_norm)
        v1=(1/g)*(pos_dest-(f*pos_origin))
        v2=(1/g)*(gp*pos_dest-pos_origin)

        # Sprawdzamy faktyczny kierunek trajektorii
        angle = np.arccos(np.dot(v2, v_dest) / (np.linalg.norm(v2)*np.linalg.norm(v_dest)))
        if angle > np.pi/2:
            # long-way wykryty → odwróć A i powtórz iterację po Z
            A = -A
            continue
        else:
            short_way_done = True
    #print()
    #print("Lambert calculations done!")
    return v1, v2

#def get_vInfinity(planetVector,shipVector):
 #   vInfinity = shipVector - planetVector
  #  return vInfinity
def get_vInfinity(planetVector, shipVector):
    # wektor nadwyżki względem planety
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
    v_arrivalSecond = get_spice_planet_vectors(planet2id, JulianArrivalCorrected) #Spice implemented
    v_arrivalSecond = np.array(v_arrivalSecond[["vx", "vy", "vz"]].iloc[0]) * 1000


    first_v = get_spice_planet_vectors(planet1id, jd) #Spice implemented
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