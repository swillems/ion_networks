#!python

# builtin
import os
import sys
import logging
import json
import time
# external
import numpy as np
import pandas as pd
import h5py
import pyteomics.mgf


BASE_PATH = os.path.dirname(os.path.dirname(__file__))
LIB_PATH = os.path.join(BASE_PATH, "lib")
DEFAULT_PARAMETER_PATH = os.path.join(LIB_PATH, "default_parameters")
DEFAULT_PARAMETER_FILES = {
    "convert": "convert_parameters.json",
    "create": "create_parameters.json",
    "evidence": "evidence_parameters.json",
    "interface": "interface_parameters.json",
    "database": "database_parameters.json",
}
DATA_TYPE_FILE_EXTENSIONS = {
    "DDA": ".mgf",
    "SONAR": "_Apex3DIons.csv",
    "HDMSE": "_Apex3DIons.csv",
    "SWIMDIA": "_Apex3DIons.csv",
    "DIAPASEF": "_centroids.hdf",
}


class open_logger(object):
    """
    Create a logger to track all progress.
    """

    def __init__(
        self,
        log_file_name,
        log_level=logging.INFO,
    ):
        """
        Create a logger to track all progress.

        Parameters
        ----------
        log_file_name : str
            If a log_file_name is provided, the current log is appended to this
            file.
        log_level : int
            The level at which log messages are returned. By default this is
            logging.INFO.
        """
        self.start_time = time.time()
        logger = logging.getLogger()
        formatter = logging.Formatter('%(asctime)s > %(message)s')
        logger.setLevel(log_level)
        if len(logger.handlers) == 0:
            console_handler = logging.StreamHandler(stream=sys.stdout)
            console_handler.setLevel(log_level)
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)
        if log_file_name != "":
            directory = os.path.dirname(log_file_name)
            if not os.path.exists(directory):
                os.makedirs(directory)
            file_handler = logging.FileHandler(log_file_name, mode="a")
            file_handler.setLevel(log_level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        self.logger = logger
        self.log_file_name = log_file_name

    def __enter__(self):
        self.logger.info("=" * 50)
        self.logger.info(
            "ion_networks.py " + " ".join(sys.argv[1:])
        )
        if self.log_file_name != "":
            self.logger.info(
                f"This log is being saved as: {self.log_file_name}."
            )
        self.logger.info("")
        return self.logger

    def __exit__(self, type, value, traceback):
        if type is not None:
            self.logger.exception("Errors occurred, execution incomplete!")
            sys.exit()
        else:
            self.logger.info("")
            self.logger.info("Successfully finished execution.")
            self.logger.info(f"Time taken: {time.time() - self.start_time}.")
        for handler in list(self.logger.handlers)[1:]:
            handler.close()
            self.logger.removeHandler(handler)


def read_parameters_from_json_file(file_name="", default=""):
    """
    Read a custom or default parameter file.

    Parameters
    ----------
    default : str
        The default parameters that should be loaded. Options are:
            "create"
            "evidence"
            "interface"
            ""
    file_name : str
        The name of a .json file that contains parameters defined by the user.
        These will override the default parameters.

    Returns
    -------
    dict
        A dictionary with parameters.
    """
    if default == "":
        parameters = {"log_file_name": ""}
    else:
        default_parameter_file_name = os.path.join(
            DEFAULT_PARAMETER_PATH,
            DEFAULT_PARAMETER_FILES[default]
        )
        with open(default_parameter_file_name, "r") as in_file:
            parameters = json.load(in_file)
    if file_name != "":
        with open(file_name, "r") as in_file:
            user_defined_parameters = json.load(in_file)
        parameters.update(user_defined_parameters)
    # TODO: Numba expects proper floats or integers, not a mixture
    # TODO: e.g. DT_error = 2.0, instead of DT_error = 2
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
    return sorted([os.path.abspath(file_name) for file_name in input_files])


def read_data_from_file(
    data_type,
    file_name,
    log_transform_intensity=True,
    logger=logging.getLogger()
):
    """
    Convert an [input_file.*] file to a pd.DataFrame with as columns the
    dimensions associated with the data type.

    Parameters
    ----------
    data_type : str
        The data type of the [input_file.*] file. Options are:
            'DDA'
            'SONAR'
            'HDMSE'
            'SWIMDIA'
            'DIAPASEF'
    file_name : str
        The file name containing centroided ions.
    log_transform_intensity : bool
        Transform the intensities to logarithmic values.
    logger : logging.logger
        The logger that indicates all progress.

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_MZ,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    if data_type == "DDA":
        read_function = read_data_from_mgf_file
    elif data_type == "SONAR":
        read_function = read_data_from_sonar_file
    elif data_type == "HDMSE":
        read_function = read_data_from_hdmse_file
    elif data_type == "SWIMDIA":
        read_function = read_data_from_swimdia_file
    elif data_type == "DIAPASEF":
        read_function = read_data_from_diapasef_file
    data = read_function(
        file_name,
        log_transform_intensity=log_transform_intensity,
        logger=logger
    )
    return data


def read_data_from_mgf_file(
    file_name,
    log_transform_intensity=True,
    logger=logging.getLogger()
):
    """
    Convert an [mgf_input.mgf] file to a pd.DataFrame with as columns the
    PRECURSOR_RT, PRECURSOR_MZ, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.

    Parameters
    ----------
    file_name : str
        The file name of the DDA .mgf file (generated with ms-convert).
    log_transform_intensity : bool
        Transform the intensities to logarithmic values.
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
    for spectrum in pyteomics.mgf.read(file_name):
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
    ints = np.concatenate(ints)
    if log_transform_intensity:
        ints = np.log2(ints)
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
    log_transform_intensity=True,
    logger=logging.getLogger()
):
    """
    Convert a [sonar_input.csv] file to a pd.DataFrame with as columns the
    PRECURSOR_RT, PRECURSOR_MZ, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.

    Parameters
    ----------
    file_name : str
        The file name of the SONAR .csv file (generated with Waters' Apex3d).
    log_transform_intensity : bool
        Transform the intensities to logarithmic values.
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
    if log_transform_intensity:
        data[:, 2] = np.log2(data[:, 2])
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
    log_transform_intensity=True,
    logger=logging.getLogger()
):
    """
    Convert a [hdmse_input.csv] file to a pd.DataFrame with as columns the
    PRECURSOR_RT, PRECURSOR_DT, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.

    Parameters
    ----------
    file_name : str
        The file name of the HDMSE .csv file (generated with Waters' Apex3d).
    log_transform_intensity : bool
        Transform the intensities to logarithmic values.
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
    if log_transform_intensity:
        data[:, 2] = np.log2(data[:, 2])
    dimensions = [
        "FRAGMENT_MZ",
        "PRECURSOR_RT",
        "FRAGMENT_LOGINT",
        "PRECURSOR_DT"
    ]
    return pd.DataFrame(data, columns=dimensions)


def read_data_from_swimdia_file(
    file_name,
    log_transform_intensity=True,
    logger=logging.getLogger()
):
    """
    Convert a [swimdia_input.csv] file to a pd.DataFrame with as columns the
    PRECURSOR_RT, PRECURSOR_DT, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.

    Parameters
    ----------
    file_name : str
        The file name of the SWIM-DIA .csv file (generated with Waters' Apex3d).
    log_transform_intensity : bool
        Transform the intensities to logarithmic values.
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
        usecols=["m_z", "rt", "mobility", "area"]
    ).values
    if log_transform_intensity:
        data[:, 2] = np.log2(data[:, 2])
    dimensions = [
        "FRAGMENT_MZ",
        "PRECURSOR_RT",
        "FRAGMENT_LOGINT",
        "PRECURSOR_DT"
    ]
    return pd.DataFrame(data, columns=dimensions)


def read_data_from_diapasef_file(
    file_name,
    min_intensity=1000,
    min_cluster_size=10,
    log_transform_intensity=True,
    logger=logging.getLogger()
):
    """
    Convert a [diapasef_input_centroids.hdf] file to a pd.DataFrame with as columns
    the PRECURSOR_RT, PRECURSOR_DT, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.

    Parameters
    ----------
    file_name : str
        The file name of the DIAPASEF _centroids.hdf file (generated with
        diapasef.py).
    min_intensity : float
        The minimimum intensity of an ion to retain it.
    min_cluster_size : int
        The minimimum cluster size of an ion to retain it.
    log_transform_intensity : bool
        Transform the intensities to logarithmic values.
    logger : logging.logger
        The logger that indicates all progress.

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_DT,
        PRECURSOR_MZ, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    logger.info(f"Reading {file_name}")
    with h5py.File(file_name, "r") as hdf_file:
        centroided_fragment_mzs = hdf_file["fragment_mz_values"][...]
        centroided_fragment_intensities = hdf_file["fragment_intensity_values"][...]
        centroided_precursor_mzs = hdf_file["precursor_mz_values"][...]
        centroided_precursor_dts = hdf_file["precursor_dt_values"][...]
        centroided_precursor_rts = hdf_file["precursor_rt_values"][...]
        cluster_sizes = hdf_file["cluster_sizes"][...]
    selection = (cluster_sizes > min_cluster_size)
    if min_intensity > 0:
        selection &= (centroided_fragment_intensities > min_intensity)
    selection = np.flatnonzero(selection)
    if log_transform_intensity:
        centroided_fragment_intensities = np.log2(
            centroided_fragment_intensities
        )
    return pd.DataFrame(
        np.stack(
            [
                centroided_fragment_mzs[selection],
                centroided_fragment_intensities[selection],
                centroided_precursor_mzs[selection],
                centroided_precursor_dts[selection],
                centroided_precursor_rts[selection] / 60,
            ]
        ).T,
        columns=[
            "FRAGMENT_MZ",
            "FRAGMENT_LOGINT",
            "PRECURSOR_MZ",
            "PRECURSOR_DT",
            "PRECURSOR_RT",
        ]
    )


def write_data_to_csv_file(
    data,
    out_file_name,
    logger=logging.getLogger()
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


def read_hdf_coordinates(
    hdf_wrapper_object,
    base_folder,
    dimensions,
    indices
):
    """
    Read coordinates from an HDF wrapper object.

    Parameters
    ----------
    hdf_wrapper_object : object
        An object that is an HDF wrapper, i.e. has a file_name that refers
        to a readable HDF file.
    base_folder : str
        The folder that contains the dimensions of interest.
    dimensions : str, iterable[str] or None
        The dimension(s) to retrieve from the HDF wrapper object. If this is
        None, a sorted list with all the available dimensions is returned.
    indices : ellipsis, slice, int, iterable[int] or iterable[bool]
        The indices that should be selected from the array. This is most
        performant when this is an ellipsis or slice, but fancy indexing
        with a mask, list or np.ndarray are also supported.

    Returns
    -------
    np.ndarray or list[np.ndarray]
        A (list of) numpy array(s) with the coordinates from the
        requested dimension(s). If dimensions is None, a sorted list with
        all the available dimensions is returned.
    """
    if dimensions is None:
        with h5py.File(hdf_wrapper_object.file_name, "r") as hdf_file:
            dimensions = sorted(hdf_file[base_folder])
    single_dimension = isinstance(dimensions, str)
    if single_dimension:
        dimensions = [dimensions]
    arrays = []
    try:
        iter(indices)
        fancy_indices = True
    except TypeError:
        fancy_indices = False
    with h5py.File(hdf_wrapper_object.file_name, "r") as hdf_file:
        for dimension in dimensions:
            array = hdf_file[base_folder][dimension]
            if fancy_indices:
                array = array[...]
            array = array[indices]
            arrays.append(array)
    if single_dimension:
        return arrays[0]
    else:
        return arrays
