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

def debug_vectors(v_planet, v_ship, r_planet, label=""):
    v_planet = np.array(v_planet)
    v_ship = np.array(v_ship)
    Vinf_vec = v_ship - v_planet
    norm_vp = np.linalg.norm(v_planet)
    norm_vs = np.linalg.norm(v_ship)
    norm_vinf = np.linalg.norm(Vinf_vec)
    cosang = np.dot(v_planet, v_ship) / (norm_vp * norm_vs)
    cosang = np.clip(cosang, -1.0, 1.0)
    ang_deg = np.degrees(np.arccos(cosang))
    print(f"--- DEBUG {label} ---")
    print("v_planet (m/s)    :", v_planet, "||", norm_vp, "m/s")
    print("v_ship  (m/s)     :", v_ship,  "||", norm_vs, "m/s")
    print("V_inf (ship-planet):", Vinf_vec, "||", norm_vinf, "m/s")
    print("angle between v_planet and v_ship:", np.round(ang_deg,2), "deg")
    # sanity:
    print("sanity: |v_ship| approx |v_planet| +/- |Vinf|? ->", np.round(norm_vs,3), "~", np.round(norm_vp + norm_vinf,3))
    print("-------------------------\n")