from astropy.time import Time

def get_julian_date(date_str: str) -> float:
    t=Time(date_str, scale='utc')
    return t.jd

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
    if planetname == "Neptun":
        planetid = 899
    if planetname == "Pluto":
        planetid = 999
    if planetname == "Voyager 1": #A taki easter egg ;)
        planetid = -31
    return planetid