#!python

import os
import sys
import logging
import json

import numpy as np
import pandas as pd
from pyteomics import mgf


BASE_PATH = os.path.dirname(os.path.dirname(__file__))
LIB_PATH = os.path.join(BASE_PATH, "lib")
DEFAULT_PARAMETER_PATH = os.path.join(LIB_PATH, "default_parameters")
DEFAULT_PARAMETER_FILES = {
    "create": "create_parameters.json",
    "evidence": "evidence_parameters.json"
}
DATA_TYPE_FILE_EXTENSIONS = {
    "dda": "mgf",
    "sonar": "csv",
    "hdmse": "csv",
    "swimdia": "csv",
}


def read_parameters_from_json_file(file_name=None, default=None):
    """
    Read a custom or default parameter file.

    Parameters
    ----------
    default : str
        The default parameters that should be loaded. Options are:
            "create"
            "evidence"
            None
    file_name : str
        The name of a .json file that contains parameters defined by the user.
        These will override the default parameters.

    Returns
    -------
    dict
        A dictionary with parameters.
    """
    if default is None:
        parameters = {}
    else:
        default_parameter_file_name = os.path.join(
            DEFAULT_PARAMETER_PATH,
            DEFAULT_PARAMETER_FILES[default]
        )
        with open(default_parameter_file_name, "r") as in_file:
            parameters = json.load(in_file)
    if file_name is not None:
        with open(file_name, "r") as in_file:
            user_defined_parameters = json.load(in_file)
        parameters.update(user_defined_parameters)
    return parameters


def get_file_names_with_extension(input_path, extension=""):
    """
    Get all file names with a specific extension from a list of files and
    folders.

    Parameters
    ----------
    input_path : iterable[str]
        An iterable with files or folders from which all files with a specific
        extension need to be selected.
    extension : str
        The extension of the files of interest.

    Returns
    -------
    list
        A sorted list with unique file names with the specific extension.
    """
    input_files = set()
    for current_path in input_path:
        if os.path.isfile(current_path):
            if current_path.endswith(extension):
                input_files.add(current_path)
        elif os.path.isdir(current_path):
            for current_file_name in os.listdir(current_path):
                if current_file_name.endswith(extension):
                    file_name = os.path.join(
                        current_path,
                        current_file_name
                    )
                    input_files.add(file_name)
    return sorted(input_files)


def read_data_from_file(
    data_type,
    file_name,
    logger=logging.getLogger('ion_network_log')
):
    """
    Convert an [input_file.*] file to a pd.DataFrame with as columns the
    dimensions associated with the data type.

    Parameters
    ----------
    data_type : str
        The data type of the [input_file.*] file. Options are:
            'dda'
            'sonar'
            'hdmse'
            'swimdia'
    file_name : str
        The file name of the DDA .mgf file (generated with ms-convert).
    logger : logging.logger
        The logger that indicates all progress.

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_MZ,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    if data_type == "dda":
        data = read_data_from_mgf_file(file_name, logger)
    elif data_type == "sonar":
        data = read_data_from_sonar_file(file_name, logger)
    elif data_type == "hdmse":
        data = read_data_from_hdmse_file(file_name, logger)
    elif data_type == "swimdia":
        data = read_data_from_swimdia_file(file_name, logger)
    return data


def read_data_from_mgf_file(
    file_name,
    logger=logging.getLogger('ion_network_log')
):
    """
    Convert an [mgf_input.mgf] file to a pd.DataFrame with as columns the
    PRECURSOR_RT, PRECURSOR_MZ, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.

    Parameters
    ----------
    file_name : str
        The file name of the DDA .mgf file (generated with ms-convert).
    logger : logging.logger
        The logger that indicates all progress.

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_MZ,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    logger.info(f"Reading {file_name}")
    mz1s = []
    mz2s = []
    rts = []
    ints = []
    for spectrum in mgf.read(file_name):
        peak_count = len(spectrum["intensity array"])
        ints.append(spectrum["intensity array"])
        mz2s.append(spectrum["m/z array"])
        rts.append(
            np.repeat(spectrum["params"]["rtinseconds"] / 60, peak_count)
        )
        mz1s.append(np.repeat(spectrum["params"]["pepmass"][0], peak_count))
    mz1s = np.concatenate(mz1s)
    mz2s = np.concatenate(mz2s)
    rts = np.concatenate(rts)
    ints = np.log(np.concatenate(ints))
    dimensions = [
        "FRAGMENT_MZ",
        "PRECURSOR_RT",
        "FRAGMENT_LOGINT",
        "PRECURSOR_MZ"
    ]
    data = np.stack([mz2s, rts, ints, mz1s]).T
    return pd.DataFrame(data, columns=dimensions)


def read_data_from_sonar_file(
    file_name,
    logger=logging.getLogger('ion_network_log')
):
    """
    Convert a [sonar_input.csv] file to a pd.DataFrame with as columns the
    PRECURSOR_RT, PRECURSOR_MZ, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.

    Parameters
    ----------
    file_name : str
        The file name of the SONAR .csv file (generated with Waters' Apex3d).
    logger : logging.logger
        The logger that indicates all progress.

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_MZ,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    logger.info(f"Reading {file_name}")
    data = pd.read_csv(
        file_name,
        engine="c",
        dtype=np.float,
        usecols=["Function", "m_z", "rt", "mobility", "area"]
    ).values
    data = data[np.searchsorted(data[:, 0], 2):, 1:]
    data[:, 2] = np.log(data[:, 2])
    data[:, 3] = 400 + data[:, 3] * (900 - 400) / 200
    dimensions = [
        "FRAGMENT_MZ",
        "PRECURSOR_RT",
        "FRAGMENT_LOGINT",
        "PRECURSOR_MZ"
    ]
    return pd.DataFrame(data, columns=dimensions)


def read_data_from_hdmse_file(
    file_name,
    logger=logging.getLogger('ion_network_log')
):
    """
    Convert a [hdmse_input.csv] file to a pd.DataFrame with as columns the
    PRECURSOR_RT, PRECURSOR_DT, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.

    Parameters
    ----------
    file_name : str
        The file name of the HDMSE .csv file (generated with Waters' Apex3d).
    logger : logging.logger
        The logger that indicates all progress.

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_DT,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    logger.info(f"Reading {file_name}")
    data = pd.read_csv(
        file_name,
        engine="c",
        dtype=np.float,
        usecols=["Function", "m_z", "rt", "mobility", "area"]
    ).values
    data = data[np.searchsorted(data[:, 0], 2):, 1:]
    data[:, 2] = np.log(data[:, 2])
    dimensions = [
        "FRAGMENT_MZ",
        "PRECURSOR_RT",
        "FRAGMENT_LOGINT",
        "PRECURSOR_DT"
    ]
    return pd.DataFrame(data, columns=dimensions)


def read_data_from_swimdia_file(
    file_name,
    logger=logging.getLogger('ion_network_log')
):
    """
    Convert a [swimdia_input.csv] file to a pd.DataFrame with as columns the
    PRECURSOR_RT, PRECURSOR_DT, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.

    Parameters
    ----------
    file_name : str
        The file name of the SWIM-DIA .csv file (generated with Waters' Apex3d).
    logger : logging.logger
        The logger that indicates all progress.

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_DT,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    logger.info(f"Reading {file_name}")
    data = pd.read_csv(
        file_name,
        engine="c",
        dtype=np.float,
        usecols=["Function", "m_z", "rt", "mobility", "area"]
    ).values
    data[:, 2] = np.log(data[:, 2])
    dimensions = [
        "FRAGMENT_MZ",
        "PRECURSOR_RT",
        "FRAGMENT_LOGINT",
        "PRECURSOR_DT"
    ]
    return pd.DataFrame(data, columns=dimensions)


def write_data_to_csv_file(
    data,
    out_file_name,
    logger=logging.getLogger('ion_network_log')
):
    """
    Save a pandas dataframe with ion coordinates to a file.

    Parameters
    ----------
    data : pd.DataFrame
        A pd.DataFrame with as columns the selection / separation dimensions.
    out_file_name : str
        The file name of the .csv file in which to save the data.
    logger : logging.logger
        The logger that indicates all progress.
    """
    logger.info(f"Writing {out_file_name}")
    data.to_csv(out_file_name, index=False)


class open_logger(object):
    """
    Create a logger to track all progress.
    """

    def __init__(
        self,
        log_file_name,
        log_level=logging.DEBUG,
        overwrite=False,
        parameters={}
    ):
        """
        Create a logger to track all progress.

        Parameters
        ----------
        log_file_name : str
            If a log_file_name is provided, the current log is appended to this
            file.
        log_level : int
            The level at which log messages are returned. by default this is
            logging.DEBUG.
        overwrite : bool
            If overwrite is True, the current log is not appended to the file
            name but overwrites it instead.
        """
        logger = logging.getLogger('ion_network_log')
        formatter = logging.Formatter('%(asctime)s > %(message)s')
        logger.setLevel(log_level)
        console_handler = logging.StreamHandler(stream=sys.stdout)
        console_handler.setLevel(log_level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        if log_file_name is None:
            log_file_name = parameters["log_file_name"]
        if log_file_name is not None:
            directory = os.path.dirname(log_file_name)
            if not os.path.exists(directory):
                os.makedirs(directory)
            if overwrite:
                mode = "w"
            else:
                mode = "a"
            file_handler = logging.FileHandler(
                log_file_name,
                mode=mode
            )
            file_handler.setLevel(log_level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            parameters["log_file_name"] = log_file_name
        else:
            parameters["log_file_name"] = ""
        logger.info("=" * 50)
        logger.info(
            "Executing command: ion_networks.py " + " ".join(sys.argv[1:])
        )
        logger.info("")
        self.logger = logger
        self.log_file_name = log_file_name

    def __enter__(self):
        return self.logger

    def __exit__(self, type, value, traceback):
        if type is not None:
            self.logger.exception("Errors occurred, execution incomplete!")
            handlers = self.logger.handlers[:]
            for handler in handlers:
                handler.close()
                self.logger.removeHandler(handler)
            sys.exit()
        else:
            self.logger.info("Successfully finished execution.")
