import numpy as np
from astropy.time import Time

def get_julian_date(date_str: str) -> float:
    t=Time(date_str, scale='utc')
    return t.jd

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
    if z > 0:
        return (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z)**3)
    elif z < 0:
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z)**3)
    else:
        return 1/6

def get_planet_id(planetname):
    if planetname == "Mercury":
        planetid = 199
    if planetname == "Venus":
        planetid = 299
    if planetname == "Earth":
        planetid = 399 #Można w przyszłości ustawić 3 zamiast 399 czyli barycentrum układu Ziemia-Księżyc a nie samą Ziemię
    if planetname == "Mars":
        planetid = 499
    if planetname == "Jupiter":
        planetid = 599
    if planetname == "Saturn":
        planetid = 699
    if planetname == "Uranus":
        planetid = 799
    if planetname == "Neptune":
        planetid = 899
    if planetname == "Pluto":
        planetid = 999
    if planetname == "Voyager 1": #A taki easter egg ;)
        planetid = -31
    return planetid