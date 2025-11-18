import sqlite3
import os
import sys

if getattr(sys, 'frozen', False):
    base_dir = sys._MEIPASS  # folder, w którym jest .exe
else:
    base_dir = os.path.dirname(os.path.abspath(__file__))  # folder skryptu w IDE

#Ścieżka do bazy danych
db_path = os.path.join(base_dir, "db", "celestial_bodies.db")

def get_planet_semimajor(planetName):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT semimajor_axis FROM celestial_bodies WHERE name = ?", (planetName,))
    planetsemimajor = cursor.fetchone()
    conn.close()
    return planetsemimajor[0]

def get_planet_GM(planetName):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT mu FROM celestial_bodies WHERE name = ?", (planetName,))
    planetGM = cursor.fetchone()
    conn.close()
    return planetGM[0]

def get_planet_Radius(planetName):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT radius FROM celestial_bodies WHERE name = ?", (planetName,))
    planetRadius = cursor.fetchone()
    conn.close()
    return planetRadius[0]

def get_planet_orbPeriod(planetName):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT orbPeriod FROM celestial_bodies WHERE name = ?", (planetName,))
    planetOrbPeriod = cursor.fetchone()
    conn.close()
    return planetOrbPeriod[0]

def get_planet_SOI(planetName):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT soi_range FROM celestial_bodies WHERE name = ?", (planetName,))
    planetSoiRange = cursor.fetchone()
    conn.close()
    return planetSoiRange[0]

def get_planet_Min_Orb_Height(planetName):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT min_orbit_height FROM celestial_bodies WHERE name = ?", (planetName,))
    planetMinOrbHeight = cursor.fetchone()
    conn.close()
    return planetMinOrbHeight[0]