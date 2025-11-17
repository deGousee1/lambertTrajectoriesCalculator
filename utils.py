import sys

import numpy as np
from astropy.time import Time
import warnings
from db import get_planet_Min_Orb_Height, get_planet_SOI
from erfa import ErfaWarning
warnings.filterwarnings('ignore', category=ErfaWarning)

def get_julian_date(date_str: str) -> float:
    t=Time(date_str, scale='utc')
    return t.jd

def julian_to_utc(jd: float) -> str:
    t = Time(jd, format='jd', scale='utc')
    return t.utc.iso.split(" ")[0]

def get_Clear_ToF_Time(correctedToF):
    daysToF = int(correctedToF // 86400)
    secondsRemaining = correctedToF % 86400
    hoursToF = int(secondsRemaining // 3600)
    secondsRemaining %= 3600
    minutesToF = int(secondsRemaining // 60)
    secsToF = secondsRemaining % 60
    return daysToF, hoursToF, minutesToF, secsToF

def stumpff_C(z):
    if z > 0:
        return (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        return (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        return 1/2

def stumpff_S(z):
    if z > 0: #Można zmusić program do użycia przyblizenia 1/6 jak z jest bardzo blisko 0. Zapobiegnie to błedom. To samo do C(z)
        return (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z)**3)
    elif z < 0:
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z)**3)
    else:
        return 1/6

def get_planet_id(planetname):
    if planetname == "Mercury":
        planetid = 1
    elif planetname == "Venus":
        planetid = 299
    elif planetname == "Earth":
        planetid = 399 #Można w przyszłości ustawić 3 zamiast 399 czyli barycentrum układu Ziemia-Księżyc a nie samą Ziemię
    elif planetname == "Mars":
        planetid = 4
    elif planetname == "Jupiter":
        planetid = 5
    elif planetname == "Saturn":
        planetid = 6
    elif planetname == "Uranus":
        planetid = 7
    elif planetname == "Neptune":
        planetid = 8
    #elif planetname == "Pluto":
    #    planetid = 9
    else:
        planetid = 10
    return planetid

def ask_for_Entry_Data():
    date = input("Date of departure (yyyy-mm-dd): ")
    date_julian = get_julian_date(date)

    planetName = input("First planet name: ")
    planet1name = planetName
    planet1id = get_planet_id(planetName)
    if planet1id == 10:
        print("Given planet name is incorrect. Please try again.")
        planetName = input("First planet name: ")
        planet1name = planetName
        planet1id = get_planet_id(planetName)
        if planet1id == 10:
            print("Given planet name is incorrect. The first letter must be a capital letter. Please try again.")
            planetName = input("First planet name: ")
            planet1name = planetName
            planet1id = get_planet_id(planetName)
            if planet1id == 10:
                print("Given planet name is incorrect. Use names from this list: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune. Please try again.")
                planetName = input("First planet name: ")
                planet1name = planetName
                planet1id = get_planet_id(planetName)
                if planet1id == 10:
                    sys.exit()
    departOrbMin = get_planet_Min_Orb_Height(planetName)
    departOrbMax = get_planet_SOI(planetName)

    planetName = input("Second planet name: ")
    planet2name = planetName
    planet2id = get_planet_id(planetName)
    if planet2id == 10:
        print("Given planet name is incorrect. Please try again.")
        planetName = input("Second planet name: ")
        planet2name = planetName
        planet2id = get_planet_id(planetName)
        if planet2id == 10:
            print("Given planet name is incorrect. The first letter must be a capital letter. Please try again.")
            planetName = input("Second planet name: ")
            planet2name = planetName
            planet2id = get_planet_id(planetName)
            if planet2id == 10:
                print("Given planet name is incorrect. Use names from this list: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune. Please try again.")
                planetName = input("Second planet name: ")
                planet2name = planetName
                planet2id = get_planet_id(planetName)
                if planet2id == 10:
                    sys.exit()
    if planet2id == planet1id:
        print("Both given planets are the same. Please try again.")
        sys.exit()
    arrivalOrbMin = get_planet_Min_Orb_Height(planetName)
    arrivalOrbMax = get_planet_SOI(planetName)

    departOrbitHeight = float(input("Departure Orbit Height (km): ")) * 1000
    arrivalOrbitHeight = float(input("Arrival Orbit Height (km): ")) * 1000
    print("Departure orbit min: ", departOrbMin)
    print("Arrival orbit min: ", arrivalOrbMin)
    print("Departure orbit max: ", departOrbMax)
    print("Arrival orbit max: ", arrivalOrbMax)
    if departOrbitHeight < departOrbMin*1000:
        departOrbitHeight = departOrbMin*1000
        print("Departure orbit height too low! New valid orbit height is set to", departOrbitHeight/1000, "km.")
    if departOrbitHeight > departOrbMax*1000:
        departOrbitHeight = departOrbMax*1000
        print("Departure orbit height too high! New valid orbit height is set to", departOrbitHeight/1000, "km.")
    if arrivalOrbitHeight < arrivalOrbMin*1000:
        arrivalOrbitHeight = arrivalOrbMin*1000
        print("Arrival orbit height too low! New valid orbit height is set to", arrivalOrbitHeight/1000, "km.")
    if arrivalOrbitHeight > arrivalOrbMax*1000:
        arrivalOrbitHeight = arrivalOrbMax*1000
        print("Arrival orbit height too high! New valid orbit height is set to", arrivalOrbitHeight/1000, "km.")
    print("Arrival orbit height: ", arrivalOrbitHeight)
    print("Departure orbit height: ", departOrbitHeight)
    return date_julian, planet1name, planet1id, planet2name, planet2id, departOrbitHeight, arrivalOrbitHeight