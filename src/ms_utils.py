#!python

# builtin
import os
import sys
import logging
import json
import time
import contextlib
import multiprocessing
# external
import numpy as np
import pandas as pd
import h5py
import pyteomics.mgf
# local
from _version import __version__ as VERSION


BASE_PATH = os.path.dirname(os.path.dirname(__file__))
LIB_PATH = os.path.join(BASE_PATH, "lib")
DEFAULT_PARAMETER_PATH = os.path.join(LIB_PATH, "default_parameters")
DEFAULT_PARAMETER_FILES = {
    "convert": "convert_parameters.json",
    "create": "create_parameters.json",
    "evidence": "evidence_parameters.json",
    "interface": "interface_parameters.json",
    "database": "database_parameters.json",
    "annotation": "annotation_parameters.json",
}
DATA_TYPE_FILE_EXTENSIONS = {
    "DDA": ".mgf",
    "SONAR": "_Apex3DIons.csv",
    "HDMSE": "_Apex3DIons.csv",
    "SWIMDIA": "_Apex3DIons.csv",
    "DIAPASEF": "_centroids.hdf",
}
LOGGER = logging.getLogger("Ion-networks")
MAX_THREADS = 1


@contextlib.contextmanager
def open_logger(log_file_name, log_level=logging.INFO):
    # TODO: Docstring
    start_time = time.time()
    formatter = logging.Formatter('%(asctime)s > %(message)s')
    LOGGER.setLevel(log_level)
    if LOGGER.hasHandlers():
        LOGGER.handlers.clear()
    console_handler = logging.StreamHandler(stream=sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    LOGGER.addHandler(console_handler)
    if log_file_name is not None:
        if log_file_name == "":
            log_file_name = BASE_PATH
        else:
            log_file_name = os.path.abspath(log_file_name)
        if os.path.isdir(log_file_name):
            log_file_name = os.path.join(log_file_name, "log.txt")
        directory = os.path.dirname(log_file_name)
        if not os.path.exists(directory):
            os.makedirs(directory)
        file_handler = logging.FileHandler(log_file_name, mode="a")
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        LOGGER.addHandler(file_handler)
    LOGGER.info("=" * 50)
    LOGGER.info(f"COMMAND: ion_networks.py {' '.join(sys.argv[1:])}")
    LOGGER.info(f"VERSION: {VERSION}")
    LOGGER.info(f"LOGFILE: {log_file_name}")
    LOGGER.info("")
    try:
        yield LOGGER
        LOGGER.info("")
        LOGGER.info("Successfully finished execution")
    except:
        LOGGER.info("")
        LOGGER.exception("Something went wrong, execution incomplete!")
    finally:
        LOGGER.info(f"Time taken: {time.time() - start_time}")
        LOGGER.info("=" * 50)
        if LOGGER.hasHandlers():
            LOGGER.handlers.clear()


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
    if "threads" in parameters:
        global MAX_THREADS
        max_cpu_count = multiprocessing.cpu_count()
        threads = parameters["threads"]
        if threads > max_cpu_count:
            MAX_THREADS = max_cpu_count
        else:
            while threads <= 0:
                threads += max_cpu_count
            MAX_THREADS = threads
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
    if not isinstance(extension, str):
        for tmp_extension in extension:
            for file_name in get_file_names_with_extension(
                input_path,
                tmp_extension
            ):
                input_files.add(file_name)
    else:
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
    )
    return data


def read_data_from_mgf_file(
    file_name,
    log_transform_intensity=True,
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

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_MZ,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    LOGGER.info(f"Reading mgf file {file_name}")
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

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_MZ,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    LOGGER.info(f"Reading sonar file {file_name}")
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

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_DT,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    LOGGER.info(f"Reading hdmse file {file_name}")
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

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_DT,
        FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    LOGGER.info(f"Reading swimdia dile {file_name}")
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
):
    """
    Convert a [diapasef_input_centroids.hdf] file to a pd.DataFrame with as
    columns the PRECURSOR_RT, PRECURSOR_DT, FRAGMENT_MZ and FRAGMENT_LOGINT
    dimensions.

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

    Returns
    -------
    pd.DataFrame
        A pd.DataFrame with as columns the PRECURSOR_RT, PRECURSOR_DT,
        PRECURSOR_MZ, FRAGMENT_MZ and FRAGMENT_LOGINT dimensions.
    """
    LOGGER.info(f"Reading diapasef file {file_name}")
    with h5py.File(file_name, "r") as hdf_file:
        centroided_fragment_mzs = hdf_file["fragment_mz_values"][...]
        centroided_fragment_intensities = hdf_file[
            "fragment_intensity_values"
        ][...]
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


def read_centroided_csv_file(
    centroided_csv_file_name,
    parameters,
):
    """
    Read a centroided .csv file and return this as a pd.DataFrame.

    Parameters
    ----------
    centroided_csv_file_name : str
        The name of a .csv file with centroided ion peaks.
    parameters : dict
        A dictionary with optional parameters for the creation of an
        ion-network.

    Returns
    -------
    pd.Dataframe
        A pd.Dataframe with centroided ion peaks.

    Raises
    -------
    KeyError
        If the PRECURSOR_RT, FRAGMENT_MZ or FRAGMENT_LOGINT column is
        missing.
    """
    LOGGER.info(f"Reading centroided csv file {centroided_csv_file_name}")
    data = pd.read_csv(
        centroided_csv_file_name,
        engine="c",
    )
    if "PRECURSOR_RT" not in data:
        raise KeyError("No PRECURSOR_RT column present")
    if "FRAGMENT_MZ" not in data:
        raise KeyError("No FRAGMENT_MZ column present")
    if "FRAGMENT_LOGINT" not in data:
        raise KeyError("No FRAGMENT_LOGINT column present")
    data.sort_values(
        by=["PRECURSOR_RT", "FRAGMENT_MZ"],
        inplace=True
    )
    return data


def write_data_to_csv_file(
    data,
    out_file_name,
):
    """
    Save a pandas dataframe with ion coordinates to a file.

    Parameters
    ----------
    data : pd.DataFrame
        A pd.DataFrame with as columns the selection / separation dimensions.
    out_file_name : str
        The file name of the .csv file in which to save the data.
    """
    LOGGER.info(f"Writing to centroided csv file {out_file_name}")
    data.to_csv(out_file_name, index=False)


class HDF_File(object):
    # TODO: Docstring

    @property
    def directory(self):
        return os.path.dirname(self.file_name)

    @property
    def file_name(self):
        return self.__file_name

    @property
    def original_file_name(self):
        return self.read_attr("original_file_name")

    @property
    def creation_time(self):
        return self.read_attr("creation_time")

    @property
    def last_updated(self):
        return self.read_attr("last_updated")

    @property
    def version(self):
        try:
            return self.read_attr("version")
        except KeyError:
            return "0.0.0"

    @property
    def is_read_only(self):
        return self.__is_read_only

    def __init__(
        self,
        file_name,
        is_read_only=True,
        new_file=False,
    ):
        # TODO: Docstring
        self.__file_name = os.path.abspath(file_name)
        if not isinstance(new_file, bool):
            raise ValueError(
                f"HDF File {self.file_name} file is not defined as (un)existing"
            )
        if new_file:
            is_read_only = False
            if not os.path.exists(self.directory):
                os.makedirs(self.directory)
            with h5py.File(self.file_name, "w") as hdf_file:
                hdf_file.attrs["creation_time"] = time.asctime()
                hdf_file.attrs["version"] = VERSION
                hdf_file.attrs["original_file_name"] = self.__file_name
                self.__update_timestamp(hdf_file)
        else:
            with h5py.File(self.file_name, "r") as hdf_file:
                if self.version != VERSION:
                    LOGGER.warning(
                        f"WARNING: {self.file_name} was created with version"
                        f" {self.version} instead of {VERSION}"
                    )
        self.__is_read_only = is_read_only

    def __eq__(self, other):
        return self.file_name == other.file_name

    def __hash__(self):
        return hash(self.file_name)

    def __str__(self):
        return f"<{self.file_name}>"

    def __repr__(self):
        return str(self)

    def __get_parent_group(self, hdf_file, parent_group_name):
        if parent_group_name == "":
            parent_group = hdf_file
        else:
            parent_group = hdf_file[parent_group_name]
        return parent_group

    def __update_timestamp(self, hdf_file):
        hdf_file.attrs["last_updated"] = time.asctime()

    def read_group(self, parent_group_name=""):
        # TODO: Docstring
        with h5py.File(self.file_name, "r") as hdf_file:
            parent_group = self.__get_parent_group(hdf_file, parent_group_name)
            group = sorted(parent_group)
        return group

    def read_attr(self, attr_key=None, parent_group_name=""):
        # TODO: Docstring
        with h5py.File(self.file_name, "r") as hdf_file:
            parent_group = self.__get_parent_group(hdf_file, parent_group_name)
            if attr_key is not None:
                attr = parent_group.attrs[attr_key]
            else:
                attr = sorted(parent_group.attrs)
        return attr

    def read_dataset(
        self,
        dataset_name,
        parent_group_name="",
        indices=Ellipsis,
        return_length=False,
        return_dtype=False,
    ):
        # TODO: Docstring
        try:
            iter(indices)
        except TypeError:
            fancy_indices = False
        else:
            fancy_indices = True
        with h5py.File(self.file_name, "r") as hdf_file:
            parent_group = self.__get_parent_group(hdf_file, parent_group_name)
            array = parent_group[dataset_name]
            if return_length:
                return len(parent_group[dataset_name])
            if return_dtype:
                return len(parent_group[dataset_name].dtype)
            if fancy_indices:
                array = array[...]
            return array[indices]

    def write_group(self, group_name, parent_group_name="", overwrite=False):
        # TODO: Docstring
        if self.is_read_only:
            raise IOError(f"HDF {self.file_name} file is opened as read only")
        with h5py.File(self.file_name, "a") as hdf_file:
            parent_group = self.__get_parent_group(hdf_file, parent_group_name)
            if group_name not in parent_group:
                hdf_group = parent_group.create_group(group_name)
            elif overwrite:
                del parent_group[group_name]
                hdf_group = parent_group.create_group(group_name)
            else:
                return
            hdf_group.attrs["creation_time"] = time.asctime()
            self.__update_timestamp(hdf_file)

    def write_attr(self, attr_key, attr_value, parent_group_name=""):
        # TODO: Docstring
        if self.is_read_only:
            raise IOError(f"HDF {self.file_name} file is opened as read only")
        with h5py.File(self.file_name, "a") as hdf_file:
            parent_group = self.__get_parent_group(hdf_file, parent_group_name)
            if isinstance(attr_value, str):
                parent_group.attrs[attr_key] = attr_value
            else:
                try:
                    iter(attr_value)
                except TypeError:
                    parent_group.attrs[attr_key] = attr_value
                else:
                    parent_group.attrs[attr_key] = str(attr_value)
            self.__update_timestamp(hdf_file)

    def write_dataset(
        self,
        dataset_name,
        dataset,
        parent_group_name="",
        overwrite=True,
        compression="lzf"
    ):
        # TODO: Docstring
        if self.is_read_only:
            raise IOError(f"HDF {self.file_name} file is opened as read only")
        if isinstance(dataset, pd.core.frame.DataFrame):
            self.write_group(dataset_name, parent_group_name, overwrite)
            for column in dataset.columns:
                self.write_dataset(
                    column,
                    dataset[column].values,
                    dataset_name,
                    overwrite,
                    compression
                )
        else:
            with h5py.File(self.file_name, "a") as hdf_file:
                parent_group = self.__get_parent_group(
                    hdf_file,
                    parent_group_name
                )
                if overwrite and (dataset_name in parent_group):
                    del parent_group[dataset_name]
                if dataset_name not in parent_group:
                    if dataset.dtype.type == np.str_:
                        dataset = dataset.astype(np.dtype('O'))
                    if dataset.dtype == np.dtype('O'):
                        hdf_dataset = parent_group.create_dataset(
                            dataset_name,
                            data=dataset,
                            compression=compression,
                            dtype=h5py.string_dtype()
                        )
                    else:
                        hdf_dataset = parent_group.create_dataset(
                            dataset_name,
                            data=dataset,
                            compression=compression,
                        )
                    hdf_dataset.attrs["creation_time"] = time.asctime()
                    self.__update_timestamp(hdf_file)
