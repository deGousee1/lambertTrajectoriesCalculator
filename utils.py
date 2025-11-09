import numpy as np
from astropy.time import Time
import warnings
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
    if planetname == "Venus":
        planetid = 299
    if planetname == "Earth":
        planetid = 399 #Można w przyszłości ustawić 3 zamiast 399 czyli barycentrum układu Ziemia-Księżyc a nie samą Ziemię
    if planetname == "Mars":
        planetid = 4
    if planetname == "Jupiter":
        planetid = 5
    if planetname == "Saturn":
        planetid = 6
    if planetname == "Uranus":
        planetid = 7
    if planetname == "Neptune":
        planetid = 8
    if planetname == "Pluto":
        planetid = 9
    return planetid

def ask_for_Entry_Data():
    date = input("Date of departure (yyyy-mm-dd): ")
    date_julian = get_julian_date(date)

    planetName = input("First planet name: ")
    planet1name = planetName
    planet1id = get_planet_id(planetName)

    planetName = input("Second planet name: ")
    planet2name = planetName
    planet2id = get_planet_id(planetName)

    departOrbitHeight = float(input("Departure Orbit Height (km): ")) * 1000
    arrivalOrbitHeight = float(input("Arrival Orbit Height (km): ")) * 1000
    return date_julian, planet1name, planet1id, planet2name, planet2id, departOrbitHeight, arrivalOrbitHeight