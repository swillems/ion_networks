#!python

import os
import logging
import functools
import time
import math

import h5py
import scipy
import numba
import pandas as pd
import numpy as np


class Network(object):
    """
    An ion-network with ions as nodes and possible precursor origin of fragments
    as edges.
    """

    def __init__(
        self,
        network_file_name,
        centroided_csv_file_name=None,
        parameters={},
        logger=logging.getLogger('ion_network_log')
    ):
        """
        Loads an ion-network. Alternatively, an ion-network is created if a
        .csv file with centroided data is provided.

        Parameters
        ----------
        network_file_name : str
            The file name of the ion-network.
        centroided_csv_file_name : str / None
            The name of a .csv file that contains the centroided data or None.
            WARNING: If not None, any existing ion-networks is overwritten with
            the data from this centroided .csv file.
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        logger : logging.logger
            The logger that indicates all progress.
        """
        self.file_name = network_file_name
        self.logger = logger
        if centroided_csv_file_name is not None:
            data = self.read_centroided_csv_file(
                centroided_csv_file_name,
                parameters
            )
            self.create_from_data(data, parameters)

    def read_centroided_csv_file(
        self,
        centroided_csv_file_name,
        parameters
    ):
        """
        Read a centroided .csv file and return as a pd.DataFrame

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
            If the RT, MZ2 or LOGINT column is missing.
        """
        self.logger.info(
            f"Reading centroided csv file {centroided_csv_file_name}."
        )
        data = pd.read_csv(
            centroided_csv_file_name,
            engine="c",
            dtype=np.float,
        )
        if "RT" not in data:
            raise KeyError("No RT column present")
        if "MZ2" not in data:
            raise KeyError("No MZ2 column present")
        if "LOGINT" not in data:
            raise KeyError("No LOGINT column present")
        data = data.sort_values(by=["RT", "MZ2"])
        return data

    def create_from_data(self, data, parameters):
        """
        Creates a .hdf file with nodes and edges from a pd.DataFrame.

        Parameters
        ----------
        data : pd.Dataframe
            A pd.Dataframe with centroided ion peaks.
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        """
        self.logger.info(f"Creating ion-network {self.file_name}.")
        directory = os.path.dirname(self.file_name)
        if not os.path.exists(directory):
            os.makedirs(directory)
        with h5py.File(self.file_name, "w") as network_file:
            self.write_nodes(network_file, data, parameters)
            self.write_edges(network_file, parameters)
            self.write_parameters(network_file, parameters)

    def write_nodes(self, network_file, data, parameters):
        """
        Saves all data as individual arrays to an .hdf file. Each array
        will be placed in a 'nodes' group and names according to its column
        name.

        Parameters
        ----------
        network_file : h5py.file
            An opened and writeable .hdf file representing the ion-network.
        data : pd.Dataframe
            A pd.Dataframe with centroided ion peaks.
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        """
        self.logger.info(f"Writing nodes of ion-network {self.file_name}.")
        node_group = network_file.create_group("nodes")
        for column in data.columns:
            node_group.create_dataset(
                column,
                data=data[column]
            )

    def write_edges(self, network_file, parameters):
        """
        Creates and saves all edges of an ion-network. They are saved in an .hdf
        file as indptr and indices of a scipy.sparse.csr matrix in en 'edge'
        group.

        Parameters
        ----------
        network_file : h5py.file
            An opened and writeable .hdf file representing the ion-network.
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        """
        self.logger.info(f"Writing edges of ion-network {self.file_name}.")
        edge_group = network_file.create_group("edges")
        dimensions = [
            dimension for dimension in self.dimensions if dimension not in [
                "RT",
                "MZ2",
                "LOGINT"
            ]
        ]
        max_deviations = [
            parameters[
                f"max_edge_deviation_{dimension}"
            ] for dimension in dimensions
        ]
        indptr, indices = self.create_sparse_edges(
            self.get_ion_coordinates("RT"),
            parameters[f"max_edge_deviation_RT"],
            tuple(self.get_ion_coordinates(dimensions)),
            tuple(max_deviations)
        )
        edge_group.create_dataset(
            "indptr",
            data=indptr
        )
        edge_group.create_dataset(
            "indices",
            data=indices,
            compression="lzf"
        )

    def write_parameters(self, network_file, parameters):
        """
        Saves all parameters of an ion-network. They are saved in an .hdf
        file as attributes of 'parameters' group.

        Parameters
        ----------
        network_file : h5py.file
            An opened and writeable .hdf file representing the ion-network.
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        """
        self.logger.info(f"Writing parameters of ion-network {self.file_name}.")
        parameter_group = network_file.create_group("parameters")
        parameter_group.attrs["creation_time"] = time.asctime()
        parameter_group.attrs["original_file_name"] = self.file_name
        parameter_group.attrs["node_count"] = len(
            network_file["edges"]["indptr"]
        ) - 1
        parameter_group.attrs["edge_count"] = len(
            network_file["edges"]["indices"]
        )
        for parameter_key, parameter_value in parameters.items():
            parameter_group.attrs[parameter_key] = parameter_value

    # @functools.lru_cache() # TODO LRU only works if dimensions always is a list
    def get_ion_coordinates(self, dimensions=None, indices=...):
        """
        Get an array with ion coordinates from the ion-network.

        Parameters
        ----------
        dimensions : str, iterable[str] or None
            The dimension(s) to retrieve from the ion-network. If this is None,
            a sorted list with all the available dimensions is returned.
        indices : ellipsis, slice, int, iterable[int] or iterable[bool]
            The indices that should be selected from the array. This is most
            performant when this is an ellipsis or slice, but fancy indexing
            with a mask, list or np.ndarray are also supported.

        Returns
        -------
        np.ndarray or list[np.ndarray]
            A (list of) numpy array(s) with the ion coordinates from the
            requested dimension(s). If dimensions is None, a sorted list with
            all the available dimensions is returned.
        """
        if dimensions is None:
            dimensions = self.dimensions
        arrays = []
        single_dimension = isinstance(dimensions, str)
        if single_dimension:
            dimensions = [dimensions]
        try:
            iter(indices)
            fancy_indices = True
        except TypeError:
            fancy_indices = False
        with h5py.File(self.file_name, "r") as network_file:
            for dimension in dimensions:
                array = network_file["nodes"][dimension]
                if fancy_indices:
                    array = array[...]
                array = array[indices]
                arrays.append(array)
        if single_dimension:
            return arrays[0]
        else:
            return arrays

    @functools.lru_cache()
    def get_edges(
        self,
        rows=...,
        columns=...,
        return_as_scipy_csr=True
    ):
        # TODO: Docstring
        # try:
        #     iter(rows)
        #     fancy_indptr = True
        # except TypeError:
        #     fancy_indptr = False
        with h5py.File(self.file_name, "r") as network_file:
            indptr = network_file["edges"]["indptr"]
            # if fancy_indptr:
            #     indptr = indptr[...]
            indptr = indptr[rows]
            indices = network_file["edges"]["indices"][...]
        if return_as_scipy_csr:
            matrix = scipy.sparse.csr_matrix(
                (
                    np.ones(self.edge_count, dtype=np.bool),
                    indices,
                    indptr
                ),
                shape=(self.node_count, self.node_count)
            )
            return matrix
        else:
            return indptr, indices

    def create_sparse_edges(
        self,
        rt_coordinate_array,
        max_rt_deviation,
        coordinate_arrays,
        max_deviations,
        symmetric=False
    ):
        # TODO: Docstring
        @numba.njit(fastmath=True)
        def numba_wrapper():
            if not symmetric:
                lower_limits = np.arange(len(rt_coordinate_array)) + 1
            else:
                lower_limits = np.searchsorted(
                    rt_coordinate_array,
                    rt_coordinate_array - max_rt_deviation,
                    "left"
                )
            upper_limits = np.searchsorted(
                rt_coordinate_array,
                rt_coordinate_array + max_rt_deviation,
                "right"
            )
            neighbors = []
            for index, (low_limit, high_limit) in enumerate(
                zip(
                    lower_limits,
                    upper_limits
                )
            ):
                small_deviations = np.repeat(True, high_limit - low_limit)
                for max_deviation, coordinate_array in zip(
                    max_deviations,
                    coordinate_arrays
                ):
                    neighbor_coordinates = coordinate_array[
                        low_limit: high_limit
                    ]
                    local_coordinate = coordinate_array[index]
                    small_coordinate_deviations = np.abs(
                        neighbor_coordinates - local_coordinate
                    ) <= max_deviation
                    small_deviations = np.logical_and(
                        small_deviations,
                        small_coordinate_deviations
                    )
                local_neighbors = low_limit + np.flatnonzero(small_deviations)
                neighbors.append(local_neighbors)
            return neighbors
        neighbors = numba_wrapper()
        indptr = np.empty(len(neighbors) + 1, dtype=np.int)
        indptr[0] = 0
        indptr[1:] = np.cumsum([len(n) for n in neighbors])
        indices = np.concatenate(neighbors)
        return indptr, indices

    def align_nodes(self, other, parameters):
        # TODO: Docstring
        self.logger.info(f"Aligning {self.file_name} with {other.file_name}.")
        dimensions = set(self.dimensions + other.dimensions)
        dimensions = [
            dimension for dimension in dimensions if dimension not in [
                "MZ2",
                "LOGINT"
            ]
        ]
        self_mz = self.get_ion_coordinates("MZ2")
        other_mz = other.get_ion_coordinates("MZ2")
        self_mz_order = np.argsort(self_mz)
        other_mz_order = np.argsort(other_mz)
        self_coordinates = tuple(
            [
                self_coordinate[
                    self_mz_order
                ] for self_coordinate in self.get_ion_coordinates(dimensions)
            ]
        )
        other_coordinates = tuple(
            [
                other_coordinate[
                    other_mz_order
                ] for other_coordinate in other.get_ion_coordinates(dimensions)
            ]
        )
        max_deviations = tuple(
            [
                parameters[
                    f"max_alignment_deviation_{dimension}"
                ] for dimension in dimensions
            ]
        )
        ppm = parameters["ppm"]
        max_mz_diff = 1 + ppm * 10**-6
        @numba.njit
        def numba_wrapper():
            low_limits = np.searchsorted(
                other_mz[other_mz_order],
                self_mz[self_mz_order] / max_mz_diff,
                "left"
            )
            high_limits = np.searchsorted(
                other_mz[other_mz_order],
                self_mz[self_mz_order] * max_mz_diff,
                "right"
            )
            results = []
            for index, (low_limit, high_limit) in enumerate(
                zip(
                    low_limits,
                    high_limits
                )
            ):
                small_deviations = np.repeat(True, high_limit - low_limit)
                for max_deviation, self_coordinate, other_coordinate in zip(
                    max_deviations,
                    self_coordinates,
                    other_coordinates
                ):
                    local_coordinate = self_coordinate[index]
                    alignment_coordinates = other_coordinate[
                        low_limit: high_limit
                    ]
                    small_coordinate_deviations = np.abs(
                        alignment_coordinates - local_coordinate
                    ) <= max_deviation
                    small_deviations = np.logical_and(
                        small_deviations,
                        small_coordinate_deviations
                    )
                matches = other_mz_order[
                    low_limit + np.flatnonzero(small_deviations)
                ]
                results.append(matches)
            return results
        results = numba_wrapper()
        first_indices = np.repeat(self_mz_order, [len(l) for l in results])
        second_indices = np.concatenate(results)
        return np.stack([first_indices, second_indices]).T

    def evidence(
        self,
        alignment_file_name,
        evidence_file_name=None,
        parameters=None
    ):
        # TODO: Docstring
        print("evidence", self, alignment_file_name, evidence_file_name, parameters)
        pass

    @property
    def node_count(self):
        """
        Get the number of nodes in the ion-network.

        Returns
        -------
        int
            The number of nodes in the ion-network.
        """
        with h5py.File(self.file_name, "r") as network_file:
            node_count = network_file["parameters"].attrs["node_count"]
        return node_count

    @property
    def dimensions(self):
        """
        Get a sorted list with all the dimensions of the nodes in the
        ion-network.

        Returns
        -------
        list[str]
            A sorted list with the names of all dimensions ion-network.
        """
        with h5py.File(self.file_name, "r") as network_file:
            dimensions = list(network_file["nodes"].keys())
        return dimensions

    @property
    def edge_count(self):
        """
        Get the number of edges in the ion-network.

        Returns
        -------
        int
            The number of edges in the ion-network.
        """
        with h5py.File(self.file_name, "r") as network_file:
            edge_count = network_file["parameters"].attrs["edge_count"]
        return edge_count

    @property
    def original_file_name(self):
        """
        Get the original file name of the ion-network.

        Returns
        -------
        str
            The original file name of the ion-network.
        """
        with h5py.File(self.file_name, "r") as network_file:
            original_file_name = network_file["parameters"].attrs[
                "original_file_name"
            ]
        return original_file_name

    @property
    def creation_time(self):
        """
        Get the creation time of the ion-network.

        Returns
        -------
        str
            The creation time of the ion-network as a string.
        """
        with h5py.File(self.file_name, "r") as network_file:
            creation_time = network_file["parameters"].attrs["creation_time"]
        return creation_time

    @property
    def key(self):
        """
        Get the unique key that identifies an ion-network.

        Returns
        -------
        str
            The unique key of an ion-network as a string.
        """
        return f"{self.edge_count}_{self.node_count}"

    def __str__(self):
        result = (
            f"Ion-network {self.file_name} with {self.node_count} nodes "
            f"and {self.edge_count} edges"
        )
        return result

    def __hash__(self):
        node_count = self.node_count
        edge_count = self.edge_count
        offset = math.ceil(math.log10(node_count + 1))
        return int(edge_count * 10**offset + node_count)

    def __cmp__(self, other):
        if self.edge_count != other.edge_count:
            if self.edge_count < other.edge_count:
                return -1
            return 1
        elif self.node_count != other.node_count:
            if self.node_count < other.node_count:
                return -1
            return 1
        return 0

    def __eq__(self, other):
        if self.__cmp__(other) == 0:
            return True
        return False

    def __ne__(self, other):
        if self.__cmp__(other) != 0:
            return True
        return False

    def __lt__(self, other):
        if self.__cmp__(other) < 0:
            return True
        return False

    def __gt__(self, other):
        if self.__cmp__(other) > 0:
            return True
        return False

    def __le__(self, other):
        if self.__cmp__(other) <= 0:
            return True
        return False

    def __ge__(self, other):
        if self.__cmp__(other) >= 0:
            return True
        return False
