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

# Wyświetlenie nagłówków i danych w czytelnej formie
print(f"{'ID':<3} {'Name':<10} {'μ (m³/s²)':<15} {'Radius (m)':<12} {'Semimajor Axis (m)':<18}")
print("-" * 60)
for row in rows:
    id, name, mu, radius, semimajor_axis = row
    print(f"{id:<3} {name:<10} {mu:<15} {radius:<12} {semimajor_axis:<18}")

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