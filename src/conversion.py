#!python

import numpy as np
import pandas as pd
from pyteomics import mgf


def read_mgf(in_file_name):
    # TODO: Docstring
    mz1s = []
    mz2s = []
    rts = []
    ints = []
    for spectrum in mgf.read(in_file_name):
        peak_count = len(spectrum["intensity array"])
        ints.append(spectrum["intensity array"])
        mz2s.append(spectrum["m/z array"])
        rts.append(np.repeat(spectrum["params"]["rtinseconds"] / 60, peak_count))
        mz1s.append(np.repeat(spectrum["params"]["pepmass"][0], peak_count))
    mz1s = np.concatenate(mz1s)
    mz2s = np.concatenate(mz2s)
    rts = np.concatenate(rts)
    ints = np.log(np.concatenate(ints))
    dimensions = ["MZ2", "RT", "LOGINT", "MZ1"]
    data = np.stack([mz2s, rts, ints, mz1s]).T
    df = pd.DataFrame(data, columns=dimensions)
    return df


def read_sonar(in_file_name):
    # TODO: Docstring
    data = pd.read_csv(
        in_file_name,
        engine="c",
        dtype=np.float,
        usecols=["Function", "m_z", "rt", "mobility", "area"]
    ).values
    data = data[np.searchsorted(data[:, 0], 2):, 1:]
    data[:, 2] = np.log(data[:, 2])
    data[:, 3] = 400 + data[:, 3] * (900 - 400) / 200
    dimensions = ["MZ2", "RT", "LOGINT", "MZ1"]
    return pd.DataFrame(data, columns=dimensions)


def read_hdmse(in_file_name):
    # TODO: Docstring
    data = pd.read_csv(
        in_file_name,
        engine="c",
        dtype=np.float,
        usecols=["Function", "m_z", "rt", "mobility", "area"]
    ).values
    data = data[np.searchsorted(data[:, 0], 2):, 1:]
    data[:, 2] = np.log(data[:, 2])
    dimensions = ["MZ2", "RT", "LOGINT", "DT"]
    return pd.DataFrame(data, columns=dimensions)


def read_swim_dia(in_file_name):
    # TODO: Docstring
    data = pd.read_csv(
        in_file_name,
        engine="c",
        dtype=np.float,
        usecols=["Function", "m_z", "rt", "mobility", "area"]
    ).values
    data[:, 2] = np.log(data[:, 2])
    dimensions = ["MZ2", "RT", "LOGINT", "DT"]
    return pd.DataFrame(data, columns=dimensions)


def write(data, out_file_name):
    # TODO: Docstring
    data.to_csv(out_file_name, index=False)
