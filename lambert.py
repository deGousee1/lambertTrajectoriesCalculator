import numpy as np
from astropy import constants as const
from db import get_planet_semimajor, get_planet_GM
from ephemerides import get_planet_vectors

AU = const.au.value


def get_ToF_estimate(planet2name, r1_norm):
    planetName = planet2name
    planet2semimajor = get_planet_semimajor(planetName)

    planetName = "Sun"
    sunGM = get_planet_GM(planetName)

    estToF = np.pi * np.sqrt( ( ( (r1_norm + planet2semimajor) / 2 ) ** 3 ) / sunGM)
    return estToF

def get_Corrected_ToF_estimate(date_julian, JulianArrivalETA, planet1id, planet2id):
    first_v = get_planet_vectors(planet1id, date_julian)
    second_v = get_planet_vectors(planet2id, JulianArrivalETA)
    r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * AU
    r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * AU
    r1_norm = np.linalg.norm(r_first)
    r2_norm = np.linalg.norm(r_second)
    planetName = "Sun"
    sunGM = get_planet_GM(planetName)
    correctedToF = np.pi * np.sqrt( ( ( (r1_norm + r2_norm) / 2 ) ** 3 ) / sunGM)
    print(first_v) #Tymczasowo do debugowania
    print(second_v) #Tymczasodo do debugowania
    return correctedToF