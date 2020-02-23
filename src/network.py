#!python

# builtin
import logging
import os
import time
import math
# external
import h5py
import scipy
import scipy.sparse
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
        logger=logging.getLogger()
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
            If the PRECURSOR_RT, FRAGMENT_MZ or FRAGMENT_LOGINT column is
            missing.
        """
        self.logger.info(
            f"Reading centroided csv file {centroided_csv_file_name}."
        )
        data = pd.read_csv(
            centroided_csv_file_name,
            engine="c",
            dtype=np.float,
        )
        if "PRECURSOR_RT" not in data:
            raise KeyError("No PRECURSOR_RT column present")
        if "FRAGMENT_MZ" not in data:
            raise KeyError("No FRAGMENT_MZ column present")
        if "FRAGMENT_LOGINT" not in data:
            raise KeyError("No FRAGMENT_LOGINT column present")
        data = data.sort_values(by=["PRECURSOR_RT", "FRAGMENT_MZ"])
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
                data=data[column],
                compression="lzf"
            )

    def write_edges(self, network_file, parameters):
        """
        Creates and saves all edges of an ion-network. They are saved in an
        .hdf file as indptr and indices of a scipy.sparse.csr matrix in an
        'edge' group.

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
                "PRECURSOR_RT",
                "FRAGMENT_MZ",
                "FRAGMENT_LOGINT"
            ]
        ]
        max_absolute_errors = [
            parameters[
                f"max_edge_absolute_error_{dimension}"
            ] for dimension in dimensions
        ]
        indptr, indices = self.create_sparse_edges(
            self.get_ion_coordinates("PRECURSOR_RT"),
            parameters[f"max_edge_absolute_error_PRECURSOR_RT"],
            tuple(self.get_ion_coordinates(dimensions)),
            tuple(max_absolute_errors)
        )
        edge_group.create_dataset(
            "indptr",
            data=indptr,
            compression="lzf"
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
        network_file.attrs["creation_time"] = time.asctime()
        network_file.attrs["node_count"] = len(
            network_file["edges"]["indptr"]
        ) - 1
        network_file.attrs["edge_count"] = len(
            network_file["edges"]["indices"]
        )
        network_file.attrs["original_file_name"] = self.file_name
        for parameter_key, parameter_value in parameters.items():
            if parameter_key.startswith("max_edge_absolute_error"):
                if parameter_key[19:] not in self.dimensions:
                    continue
            network_file.attrs[parameter_key] = parameter_value

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
        single_dimension = isinstance(dimensions, str)
        if single_dimension:
            dimensions = [dimensions]
        arrays = []
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

    def get_edges(
        self,
        return_as_scipy_csr=True,
        symmetric=False,
        data_as_index=False
    ):
        with h5py.File(self.file_name, "r") as network_file:
            indptr = network_file["edges"]["indptr"][...]
            indices = network_file["edges"]["indices"][...]
        if return_as_scipy_csr or symmetric:
            matrix = scipy.sparse.csr_matrix(
                (
                    np.ones(self.edge_count, dtype=np.bool),
                    indices,
                    indptr
                ),
                shape=(self.node_count, self.node_count)
            )
            if symmetric and return_as_scipy_csr:
                return matrix + matrix.T
            elif return_as_scipy_csr:
                if data_as_index:
                    matrix.data = np.arange(matrix.nnz)
                return matrix
            else:
                second_indices = np.repeat(
                    np.arange(self.node_count),
                    np.diff(matrix.indptr)
                )
                return second_indices, matrix.indices
        else:
            second_indices = np.repeat(
                np.arange(self.node_count),
                np.diff(indptr)
            )
            return second_indices, indices

    def create_sparse_edges(
        self,
        rt_coordinate_array,
        max_rt_absolute_error,
        coordinate_arrays,
        max_absolute_errors,
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
                    rt_coordinate_array - max_rt_absolute_error,
                    "left"
                )
            upper_limits = np.searchsorted(
                rt_coordinate_array,
                rt_coordinate_array + max_rt_absolute_error,
                "right"
            )
            neighbors = []
            for index, (low_limit, high_limit) in enumerate(
                zip(
                    lower_limits,
                    upper_limits
                )
            ):
                small_absolute_errors = np.repeat(True, high_limit - low_limit)
                for max_absolute_error, coordinate_array in zip(
                    max_absolute_errors,
                    coordinate_arrays
                ):
                    neighbor_coordinates = coordinate_array[
                        low_limit: high_limit
                    ]
                    local_coordinate = coordinate_array[index]
                    small_coordinate_absolute_errors = np.abs(
                        neighbor_coordinates - local_coordinate
                    ) <= max_absolute_error
                    small_absolute_errors = np.logical_and(
                        small_absolute_errors,
                        small_coordinate_absolute_errors
                    )
                local_neighbors = low_limit + np.flatnonzero(small_absolute_errors)
                neighbors.append(local_neighbors)
            return neighbors
        neighbors = numba_wrapper()
        indptr = np.empty(len(neighbors) + 1, dtype=np.int)
        indptr[0] = 0
        indptr[1:] = np.cumsum([len(n) for n in neighbors])
        indices = np.concatenate(neighbors)
        return indptr, indices

    def align_nodes(self, other, parameters, return_as_scipy_csr=True):
        # TODO: Docstring
        self.logger.info(f"Aligning {self.file_name} with {other.file_name}.")
        dimensions = self.dimension_overlap(
            other,
            remove_fragment_mz=True,
            remove_fragment_logint=True,
            remove_precursor_rt=parameters["calibrate_PRECURSOR_RT"]
        )
        self_mz = self.get_ion_coordinates("FRAGMENT_MZ")
        other_mz = other.get_ion_coordinates("FRAGMENT_MZ")
        self_mz_order = np.argsort(self_mz)
        other_mz_order = np.argsort(other_mz)
        self_coordinates = [
            self_coordinate[
                self_mz_order
            ] for self_coordinate in self.get_ion_coordinates(dimensions)
        ]
        other_coordinates = [
            other_coordinate[
                other_mz_order
            ] for other_coordinate in other.get_ion_coordinates(dimensions)
        ]
        max_absolute_errors = [
            parameters[
                f"max_alignment_absolute_error_{dimension}"
            ] for dimension in dimensions
        ]
        ppm = parameters["max_alignment_ppm_FRAGMENT_MZ"]
        max_mz_diff = 1 + ppm * 10**-6
        if parameters["calibrate_PRECURSOR_RT"]:
            self_coordinates += [
                self.get_ion_coordinates("PRECURSOR_RT")[
                    self_mz_order
                ]
            ]
            other_coordinates += [
                other.calibrate_precursor_rt(self, parameters)[
                    other_mz_order
                ]
            ]
            max_absolute_errors += [
                parameters["max_alignment_absolute_error_PRECURSOR_RT"]
            ]
        self_coordinates = tuple(self_coordinates)
        other_coordinates = tuple(other_coordinates)
        max_absolute_errors = tuple(max_absolute_errors)
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
                small_absolute_errors = np.repeat(True, high_limit - low_limit)
                for max_absolute_error, self_coordinate, other_coordinate in zip(
                    max_absolute_errors,
                    self_coordinates,
                    other_coordinates
                ):
                    local_coordinate = self_coordinate[index]
                    alignment_coordinates = other_coordinate[
                        low_limit: high_limit
                    ]
                    small_coordinate_absolute_errors = np.abs(
                        alignment_coordinates - local_coordinate
                    ) <= max_absolute_error
                    small_absolute_errors = np.logical_and(
                        small_absolute_errors,
                        small_coordinate_absolute_errors
                    )
                matches = other_mz_order[
                    low_limit + np.flatnonzero(small_absolute_errors)
                ]
                results.append(matches)
            return results
        results = numba_wrapper()
        self_indices = np.repeat(self_mz_order, [len(l) for l in results])
        other_indices = np.concatenate(results)
        if return_as_scipy_csr:
            return scipy.sparse.csr_matrix(
                (
                    np.ones(self_indices.shape[0], dtype=np.bool),
                    (self_indices, other_indices)
                ),
                shape=(self.node_count, other.node_count)
            )
        return self_indices, other_indices

    def quick_align(self, other, ppm):
        # TODO:
        self_mzs = self.get_ion_coordinates("FRAGMENT_MZ")
        other_mzs = other.get_ion_coordinates("FRAGMENT_MZ")
        self_mz_order = np.argsort(self_mzs)
        other_mz_order = np.argsort(other_mzs)
        max_mz_diff = 1 + ppm * 10**-6
        low_limits = np.searchsorted(
            self_mzs[self_mz_order],
            other_mzs[other_mz_order] / max_mz_diff,
            "left"
        )
        high_limits = np.searchsorted(
            self_mzs[self_mz_order],
            other_mzs[other_mz_order] * max_mz_diff,
            "right"
        )
        other_rt_order = np.argsort(other_mz_order)
        self_indices = np.concatenate(
            [
                self_mz_order[l:h] for l, h in zip(
                    low_limits[other_rt_order],
                    high_limits[other_rt_order]
                )
            ]
        )
        other_indices = np.repeat(
            np.arange(len(other_rt_order)),
            high_limits[other_rt_order] - low_limits[other_rt_order]
        )
        selection = longest_increasing_subsequence(self_indices)
        self_indices_mask = np.empty(len(selection) + 2, dtype=int)
        self_indices_mask[0] = 0
        self_indices_mask[1: -1] = self_indices[selection]
        self_indices_mask[-1] = len(self_mzs) - 1
        other_indices_mask = np.empty(len(selection) + 2, dtype=int)
        other_indices_mask[0] = 0
        other_indices_mask[1: -1] = other_indices[selection]
        other_indices_mask[-1] = len(other_mzs) - 1
        return self_indices_mask, other_indices_mask

    def calibrate_precursor_rt(self, other, parameters):
        # TODO:
        self.logger.info(
            f"Calibrating PRECURSOR_RT of {self.file_name} with "
            f"{other.file_name}."
        )
        self_indices, other_indices = self.quick_align(
            other,
            parameters["calibration_ppm_FRAGMENT_MZ"]
        )
        self_rts = self.get_ion_coordinates("PRECURSOR_RT")
        other_rts = other.get_ion_coordinates(
            "PRECURSOR_RT",
            indices=other_indices
        )
        calibrated_self_rts = []
        for self_start_index, self_end_index, other_rt_start, other_rt_end in zip(
            self_indices[:-1],
            self_indices[1:],
            other_rts[:-1],
            other_rts[1:]
        ):
            self_rt_start = self_rts[self_start_index]
            self_rt_end = self_rts[self_end_index]
            if self_rt_start == self_rt_end:
                new_rts = np.repeat(other_rt_start, self_end_index - self_start_index)
            else:
                slope = (other_rt_end - other_rt_start) / (self_rt_end - self_rt_start)
                new_rts = other_rt_start + slope * (
                    self_rts[self_start_index: self_end_index] - self_rt_start
                )
            calibrated_self_rts.append(new_rts)
        calibrated_self_rts.append([other_rts[-1]])
        calibrated_self_rts = np.concatenate(calibrated_self_rts)
        return calibrated_self_rts

    def dimension_overlap(
        self,
        other,
        remove_fragment_mz=False,
        remove_fragment_logint=False,
        remove_precursor_rt=False
    ):
        """
        Get the overlapping dimensions of this ion-network with another
        ion-network.

        Parameters
        ----------
        other : ion_network
            The second ion-network.
        remove_fragment_mz : bool
            If True, remove 'FRAGMENT_MZ' from the overlapping dimensions.
        remove_fragment_logint : bool
            If True, remove 'FRAGMENT_LOGINT' from the overlapping dimensions.
        remove_precursor_rt : bool
            If True, remove 'PRECURSOR_RT' from the overlapping dimensions.

        Returns
        -------
        list[str]
            A sorted list with all the overlapping dimensions.
        """
        dimensions = set(self.dimensions + other.dimensions)
        if remove_fragment_mz:
            dimensions.remove("FRAGMENT_MZ")
        if remove_fragment_logint:
            dimensions.remove("FRAGMENT_LOGINT")
        if remove_precursor_rt:
            dimensions.remove("PRECURSOR_RT")
        return sorted(dimensions)

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
            node_count = network_file.attrs["node_count"]
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
            dimensions = list(network_file["nodes"])
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
            edge_count = network_file.attrs["edge_count"]
        return edge_count

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
            creation_time = network_file.attrs["creation_time"]
        return creation_time

    @property
    def original_file_name(self):
        """
        Get the original file name of the ion-network.

        Returns
        -------
        str
            The original file name of the ion-network as a string.
        """
        with h5py.File(self.file_name, "r") as network_file:
            original_file_name = network_file.attrs["original_file_name"]
        return original_file_name

    @property
    def file_name_base(self):
        """
        Get the file name base of the ion-network.

        Returns
        -------
        str
            The file name base of the ion-network as a string.
        """
        return os.path.basename(self.original_file_name)[:-9]

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
        return (
            f"Ion-network {self.file_name} with {self.node_count} nodes "
            f"and {self.edge_count} edges"
        )

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


@numba.njit(fastmath=True)
def longest_increasing_subsequence(sequence):
    # TODO:
    M = np.repeat(0, len(sequence) + 1)
    P = np.repeat(0, len(sequence))
    max_subsequence_length = 0
    for current_index, current_element in enumerate(sequence):
        low_bound = 1
        high_bound = max_subsequence_length
        while low_bound <= high_bound:
            mid = (low_bound + high_bound) // 2
            if sequence[M[mid]] <= current_element:
                low_bound = mid + 1
            else:
                high_bound = mid - 1
        subsequence_length = low_bound
        P[current_index] = M[subsequence_length - 1]
        M[subsequence_length] = current_index
        if subsequence_length > max_subsequence_length:
            max_subsequence_length = subsequence_length
    longest_increasing_subsequence = np.repeat(0, max_subsequence_length)
    index = M[max_subsequence_length]
    for current_index in range(max_subsequence_length - 1, -1, -1):
        longest_increasing_subsequence[current_index] = index
        index = P[index]
    # good = np.repeat(True, max_subsequence_length)
    # good[1:] = sequence[longest_increasing_subsequence[:-1]] != sequence[longest_increasing_subsequence[1:]]
    # return longest_increasing_subsequence[good]
    return longest_increasing_subsequence
