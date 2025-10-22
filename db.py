import sqlite3
import os

# Ścieżka do bazy danych
db_path = os.path.join("db", "celestial_bodies.db")

# Połączenie z bazą
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Wczytanie wszystkich danych z tabeli
cursor.execute("SELECT * FROM celestial_bodies")
rows = cursor.fetchall()

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