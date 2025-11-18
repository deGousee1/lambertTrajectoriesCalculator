import sys
import spiceypy as spice
import pandas as pd

def get_spice_planet_vectors(planet_id: int, date: float):
    # Konwersja czasu Juliańskiego na Ephemeris Time którego używa baza danych SPICE
    # Wcześniej korzystałem z innej bazy danych JPL Horizons ale wymagała ona pobierania danych z ich serwera co spowalniało program
    # Poprzednia baza wektorów używała daty juliańskiej więc zamiast zmieniać system dat w całym programie po prostu dodałem tą konwersję formatu daty
    et = (date - 2451545.0) * 86400.0

    target = str(planet_id)
    observer = "0"  # barycentrum układu słonecznego
    ref_frame = "ECLIPJ2000"
    abcorr = "NONE"

    try:
        state, lt = spice.spkezr(target, et, ref_frame, abcorr, observer)
    except Exception:
        print("No planet vectors available for this date")
        sys.exit()

    #Ustalenie formatu w jakim zwracane są dane (ustawione w taki format, jaki zwracała poprzednia baza danych bo pod to został napisany program)
    df = pd.DataFrame([{
        "x": state[0],
        "y": state[1],
        "z": state[2],
        "vx": state[3],
        "vy": state[4],
        "vz": state[5]
    }])
    return df