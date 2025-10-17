import sqlite3

# Tworzymy plik bazy (jeśli nie istnieje)
conn = sqlite3.connect("celestial_bodies.db")
cursor = conn.cursor()
cursor.execute("SELECT name, radius, mu FROM celestial_bodies WHERE mu IS NOT NULL")
for row in cursor.fetchall():
    print(row)

conn.close()