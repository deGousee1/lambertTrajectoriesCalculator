import sys
import numpy as np
from astropy.time import Time
from ephemerides import get_spice_planet_vectors
from lambert import get_LambertV, get_Optimal_Launch_Angle, get_Delta_V
import matplotlib.pyplot as plt
from utils import julian_to_utc

def transfer_Angle_Scan(Tsyn, date_julian, planet1id, planet2id, planet2name, correctedToFdays, outward):
    scanRange = Tsyn * 0.5
    scanStep = round(Tsyn * 0.01)
    start_jd_array = np.arange(date_julian - scanRange, date_julian + scanRange, scanStep)
    iterationGoal = len(start_jd_array)
    utc_dates = [Time(jd, format='jd').to_datetime() for jd in start_jd_array]
    angleLoopCounter = 0
    angleDifference = 0
    angle_matrix = np.zeros(len(start_jd_array))
    for i, jd_start_val in enumerate(start_jd_array):
        jd_start_val: float = jd_start_val
        first_v = get_spice_planet_vectors(planet1id, jd_start_val)  # Spice implemented
        second_v = get_spice_planet_vectors(planet2id, jd_start_val)  # Spice implemented

        r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000
        r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
        r1_norm = np.linalg.norm(r_first)
        r2_norm = np.linalg.norm(r_second)

        optimalAngle = abs(get_Optimal_Launch_Angle(planet2name, correctedToFdays, outward))
        realAngle = np.degrees(np.arccos(np.dot(r_first, r_second) / (r1_norm * r2_norm)))
        angleDifference = abs(optimalAngle - realAngle)
        angle_matrix[i] = angleDifference
        angleLoopCounter += 1
        sys.stdout.write(
            f"\rTransfer window search progress: {angleLoopCounter} of {round(iterationGoal)}"
        )
        sys.stdout.flush()
    print()
    plt.plot(utc_dates, angle_matrix, marker='o')
    plt.xlabel('Maneuver start date', fontsize=11)
    plt.ylabel('Angle difference')
    plt.title('Optimal transfer angle deviation plot')
    plt.show()
    min_idx = np.argmin(angle_matrix)
    min_angle = angle_matrix[min_idx]

    angleDifference = 0
    jd_start_val = date_julian - scanRange
    jd_start_val: float = jd_start_val
    currentBestAngle = 181
    worstAngle = 0
    bestAngleDateJ = 0
    currentBestAngleDateJ = 0
    worstAngleDateJ = 0
    windowFound = False
    optimalAngle = abs(get_Optimal_Launch_Angle(planet2name, correctedToFdays, outward))

    for jdDate in start_jd_array:  # Finding maximum angle difference
        jdDate: float = jdDate
        first_v = get_spice_planet_vectors(planet1id, jdDate)  # Spice implemented
        second_v = get_spice_planet_vectors(planet2id, jdDate)  # Spice implemented

        r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000
        r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
        v_firstA = np.array(first_v[["vx", "vy", "vz"]].iloc[0]) * 1000
        v_secondA = np.array(second_v[["vx", "vy", "vz"]].iloc[0]) * 1000
        r1_norm = np.linalg.norm(r_first)
        r2_norm = np.linalg.norm(r_second)

        realAngle = np.degrees(np.arccos(np.dot(r_first, r_second) / (r1_norm * r2_norm)))
        angleDifference = abs(optimalAngle - realAngle)
        if optimalAngle < 100 and optimalAngle > 80:
            if angleDifference > worstAngle:
                velocityVectorAngle = np.arccos(
                    np.dot(v_firstA, v_secondA) / (np.linalg.norm(v_firstA) * np.linalg.norm(v_secondA)))
                utcScanDateDebug = julian_to_utc(jdDate)
                utcAngleDateDebug = julian_to_utc(worstAngleDateJ)
                if np.degrees(velocityVectorAngle) > 100:
                    if angleDifference > worstAngle:
                        worstAngle = angleDifference
                        worstAngleDateJ = jdDate
        else:
            if angleDifference > worstAngle:
                worstAngle = angleDifference
                worstAngleDateJ = jdDate
        angleLoopCounter += 1
        sys.stdout.write(
            f"\rTransfer window search progress: {angleLoopCounter} of {round(iterationGoal)}"
        )
        sys.stdout.flush()
    utcWorstAngleDate = julian_to_utc(worstAngleDateJ)

    start_jd_array = np.arange(worstAngleDateJ, worstAngleDateJ + scanRange * 2, scanStep)

    for jdDate in start_jd_array:
        if windowFound == False:
            jdDate: float = jdDate
            first_v = get_spice_planet_vectors(planet1id, jdDate)  # Spice implemented
            second_v = get_spice_planet_vectors(planet2id, jdDate)  # Spice implemented

            r_first = np.array(first_v[["x", "y", "z"]].iloc[0]) * 1000
            r_second = np.array(second_v[["x", "y", "z"]].iloc[0]) * 1000
            r1_norm = np.linalg.norm(r_first)
            r2_norm = np.linalg.norm(r_second)
            optimalAngle = abs(get_Optimal_Launch_Angle(planet2name, correctedToFdays, outward))
            realAngle = np.degrees(np.arccos(np.dot(r_first, r_second) / (r1_norm * r2_norm)))
            angleDifference = abs(optimalAngle - realAngle)
            if angleDifference < currentBestAngle:
                currentBestAngle = angleDifference
                currentBestAngleDateJ = jdDate
                if outward == True:
                    if currentBestAngle < 5:
                        bestAngle = currentBestAngle
                        bestAngleDateJ = jdDate
                        windowFound = True
        angleLoopCounter += 1
        sys.stdout.write(
            f"\rTransfer window search progress: {angleLoopCounter} of {round(iterationGoal)}"
        )
        sys.stdout.flush()
    print()
    bestAngleDateJ = currentBestAngleDateJ
    bestAngle = currentBestAngle
    if outward == True:
        bestAngleDateJ = currentBestAngleDateJ
        bestAngle = currentBestAngle

    utcTransferWindow = julian_to_utc(bestAngleDateJ)
    return bestAngleDateJ, bestAngle, utcTransferWindow, worstAngle, utcWorstAngleDate

def porkchop_plot(scanRange, scanStep, scanStepToF,bestAngleDateJ, correctedToFdays, planet1id, planet2id, planet1name, planet2name, departOrbitHeight, arrivalOrbitHeight, porkchopNumber):
    if porkchopNumber == 1:
        start_jd_array = np.arange(bestAngleDateJ - scanRange, bestAngleDateJ + scanRange, scanStep)
        tof_days_array = np.arange(correctedToFdays * 0.5, correctedToFdays + scanRange, scanStepToF)
    else:
        start_jd_array = np.arange(bestAngleDateJ - scanRange, bestAngleDateJ + scanRange, scanStep)
        tof_days_array = np.arange(correctedToFdays - scanRange, correctedToFdays + scanRange, scanStep)
    print("Scan range: ", scanRange, "Scan step: ", scanStep, "Scan step ToF", scanStepToF)
    iterationGoal = len(start_jd_array) * len(tof_days_array)
    utc_dates = [Time(jd, format='jd').to_datetime() for jd in start_jd_array]
    firstLoopCounter = 0
    deltaV_matrix = np.zeros((len(tof_days_array), len(start_jd_array)))
    for i, tof in enumerate(tof_days_array):
        for j, jd_start_val in enumerate(start_jd_array):
            tof: float = tof
            jd_arrival = jd_start_val + tof
            v1, v2 = get_LambertV(
                JulianArrivalCorrected=jd_arrival,
                date_julian=jd_start_val,
                planet1id=planet1id,
                planet2id=planet2id,
                correctedToF=tof * 86400
            )
            # Obliczenia deltaV
            jd_arrival = jd_start_val + tof
            departDeltaV, arrivalDeltaV = get_Delta_V(planet2id=planet2id,
                                                      planet1id=planet1id
                                                      , v1=v1
                                                      , v2=v2
                                                      , planet1name=planet1name
                                                      , planet2name=planet2name
                                                      , departOrbitHeight=departOrbitHeight
                                                      , arrivalOrbitHeight=arrivalOrbitHeight
                                                      , JulianArrivalCorrected=jd_arrival
                                                      , jd=jd_start_val)
            deltaV = departDeltaV + arrivalDeltaV
            deltaV_matrix[i, j] = deltaV
            firstLoopCounter += 1
            sys.stdout.write(
                f"\rPorkchop iteration progress: {firstLoopCounter} of {round(iterationGoal)}"
            )
            sys.stdout.flush()

    print()
    print("First rough sieve porkchop graph done!")

    plt.figure(figsize=(10, 6))
    X, Y = np.meshgrid(utc_dates, tof_days_array)

    min_idx = np.unravel_index(np.argmin(deltaV_matrix), deltaV_matrix.shape)
    i_min, j_min = min_idx
    best_tof = tof_days_array[i_min]
    jd = float(start_jd_array[j_min])
    utcBestLaunch = julian_to_utc(jd)
    best_deltaV = deltaV_matrix[i_min, j_min]
    print("Best time of flight:", best_tof, "Best launch date:", utcBestLaunch, "Best deltaV possible:", best_deltaV)

    deltaV_matrix_masked = np.ma.masked_greater(deltaV_matrix, best_deltaV * 2)
    plt.contourf(X, Y, deltaV_matrix_masked, levels=50, cmap='viridis')
    plt.colorbar(label='Delta-V [m/s]')
    plt.xlabel('Data startu')
    plt.ylabel('Czas lotu [dni]')
    plt.title('Porkchop plot')
    plt.show()
    return jd, best_tof, best_deltaV