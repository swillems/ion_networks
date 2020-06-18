#!python

import numpy as np
import numba
import multiprocessing as mp










def get_mzs_and_frequencies(mzs, ppm):
    new_mzs = np.sort(mzs)
    lower = np.searchsorted(new_mzs, new_mzs * (1 - ppm / 10**6), "left")
    upper = np.searchsorted(new_mzs, new_mzs * (1 + ppm / 10**6), "right")
    return new_mzs, upper - lower


@numba.njit
def find_peaks(intensities, mzs, max_distance=0.9):
    peaks = np.zeros(int(mzs[-1]), np.int64)
    current_max_mz = 0
    current_max_int = 0
    current_max_index = 0
    for index, (intensity, mz) in enumerate(zip(intensities, mzs)):
        if mz > current_max_mz + max_distance:
            peaks[int(current_max_mz)] = current_max_index
            current_max_mz = mz
            current_max_int = intensity
            current_max_index = index
        elif intensity > current_max_int:
            current_max_mz = mz
            current_max_int = intensity
            current_max_index = index
    return peaks








































# for f in /media/proteomics2/RAW/Synapt/HDMSe/Covid19_eSwab_Dilution/*blanco*; do apex3d -pRawDirName "$f" -outputDirName /home/sander/data/covid_eswab_hdmse/ -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -bEnableCentroids 0; done
