import numpy as np

from db import get_planet_semimajor, get_planet_GM


def get_ToF_estimate(planet1name, planet2name):
    planetName = planet1name
    planet1semimajor = get_planet_semimajor(planetName)

    planetName = planet2name
    planet2semimajor = get_planet_semimajor(planetName)

    planetName = "Sun"
    sunGM = get_planet_GM(planetName)

    estToF = np.pi * np.sqrt( ( ( (planet1semimajor + planet2semimajor) / 2 ) ** 3 ) / sunGM)
    return estToF
