#!python

# builtin
import os
# external
import numpy as np
import scipy
import numba
import h5py
import scipy.sparse
# local
import ms_database
import ms_utils


class HDF_MS_Run_File(ms_utils.HDF_File):
    # TODO: Docstring

    @staticmethod
    def convert_reference_to_trimmed_file_name(reference):
        # TODO: Docstring
        if not isinstance(reference, str):
            try:
                reference = reference.file_name
            except AttributeError:
                raise ValueError("Invalid reference for HDF file")
        reference_components = reference.split(".")
        if not (
            (
                len(reference_components) >= 3
            ) and (
                reference_components[-1] in ["hdf", "csv"]
            ) and (
                reference_components[-2] in ["inet", "evidence", "annotation"]
            )
        ):
            raise ValueError("Invalid reference for HDF file")
        return ".".join(reference_components[:-2])

    @property
    def run_name(self):
        # TODO: Docstring
        return ".".join(os.path.basename(self.file_name).split(".")[:-2])


class Network(HDF_MS_Run_File):
    """
    An ion-network with ions as nodes and possible precursor origin of
    fragments as edges.
    """

    def __init__(
        self,
        reference,
        is_read_only=True,
        new_file=False,
    ):
        # TODO: Docstring
        file_name = self.convert_reference_to_trimmed_file_name(reference)
        super().__init__(
            f"{file_name}.inet.hdf",
            is_read_only=is_read_only,
            new_file=new_file,
        )

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
        ms_utils.LOGGER.info(f"Creating ion-network {self.file_name}")
        self.write_nodes(data)
        self.write_edges(parameters)
        self.write_parameters(parameters)

    def write_nodes(self, data):
        """
        Saves all data as individual arrays to an .hdf file. Each array
        will be placed in a 'nodes' group and names according to its column
        name.

        Parameters
        ----------
        data : pd.Dataframe
            A pd.Dataframe with centroided ion peaks.
        """
        ms_utils.LOGGER.info(f"Writing nodes of ion-network {self.file_name}")
        regular_columns = [
            column for column in data.columns if not column.startswith("#")
        ]
        self.create_dataset("nodes", data[regular_columns])
        if len(regular_columns) != data.shape[1]:
            comment_columns = [
                column for column in data.columns if column.startswith("#")
            ]
            self.create_dataset("comments", data[comment_columns])

    def write_edges(self, parameters):
        """
        Creates and saves all edges of an ion-network. They are saved in an
        .hdf file as indptr and indices of a scipy.sparse.csr matrix in an
        'edge' group. Edges are always saved such that rt1 <= rt2.

        Parameters
        ----------
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        """
        ms_utils.LOGGER.info(f"Writing edges of ion-network {self.file_name}")
        self.create_group("edges")
        dimensions = self.get_precursor_dimensions(
            remove_precursor_rt=True
        )
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
        self.create_dataset(
            "indptr",
            indptr,
            parent_group_name="edges"
        )
        self.create_dataset(
            "indices",
            indices,
            parent_group_name="edges"
        )

    def write_parameters(self, parameters):
        """
        Saves all parameters of an ion-network. They are saved in an .hdf
        file as attributes of 'parameters' group.

        Parameters
        ----------
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        """
        ms_utils.LOGGER.info(
            f"Writing parameters of ion-network {self.file_name}"
        )
        self.create_attr(
            "node_count",
            self.get_dataset("FRAGMENT_MZ", "nodes", return_length=True)
        )
        self.create_attr(
            "edge_count",
            self.get_dataset("indices", "edges", return_length=True)
        )
        for parameter_key, parameter_value in parameters.items():
            if parameter_key.startswith("max_edge_absolute_error"):
                if parameter_key[19:] not in self.dimensions:
                    continue
            self.create_attr(parameter_key, parameter_value)

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

    def get_ion_comments(self, dimensions=None, indices=...):
        """
        Get an array with ion comments from the ion-network.

        Parameters
        ----------
        dimensions : str, iterable[str] or None
            The comment dimension(s) to retrieve from the ion-network. If this
            is None, a sorted list with all the available comment dimensions is
            returned.
        indices : ellipsis, slice, int, iterable[int] or iterable[bool]
            The indices that should be selected from the array. This is most
            performant when this is an ellipsis or slice, but fancy indexing
            with a mask, list or np.ndarray are also supported.

        Returns
        -------
        np.ndarray or list[np.ndarray]
            A (list of) numpy array(s) with the ion comments from the
            requested comment dimension(s). If dimensions is None, a sorted
            list with all the available comment dimensions is returned.
        """
        if dimensions is None:
            dimensions = self.comment_dimensions
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
                array = network_file["comments"][dimension]
                if fancy_indices:
                    array = array[...]
                array = array[indices]
                arrays.append(array)
        if single_dimension:
            return arrays[0]
        else:
            return arrays
        # TODO make generic factory function with wrappers?

    def get_edges(
        self,
        return_as_scipy_csr=True,
        symmetric=False,
        data_as_index=False,
        indptr_and_indices=False,
    ):
        # TODO: Docstring
        with h5py.File(self.file_name, "r") as network_file:
            indptr = network_file["edges"]["indptr"][...]
            indices = network_file["edges"]["indices"][...]
        if indptr_and_indices:
            return indptr, indices
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

    @staticmethod
    @numba.njit(cache=True)
    def create_sparse_edges(
        rt_coordinate_array,
        max_rt_absolute_error,
        coordinate_arrays,
        max_absolute_errors,
        symmetric=False
    ):
        """
        Create all edges of an ion-network and return as indpt and indices of
        a sparse matrix.

        Parameters
        ----------
        rt_coordinate_array : np.ndarray
            An array with the individual rt coordinates of each ion.
        max_rt_absolute_error : float
            The maximum allowed rt error between two ions.
        coordinate_arrays : tuple[np.ndarray]
            A tuple with the individual coordinates of each ion per dimension.
        max_absolute_errors : typle[float]
            A tuple with the maximum allowed errors for each dimension.
        symmetric : bool
            If False, only create edges in one direction so that they satisfy
            rt1 <= rt2.

        Returns
        -------
        tuple[np.ndarray]
            A tuple with as first element the indptr and as second
            element the indices, equal to those of a scipy.csr.matrix.
        """
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
        indptr = np.empty(len(neighbors) + 1, np.int64)
        lens = np.empty(len(neighbors), np.int64)
        for i, n in enumerate(neighbors):
            lens[i] = len(n)
        indptr[0] = 0
        indptr[1:] = np.cumsum(lens)
        indices = np.empty(indptr[-1], np.int64)
        for s, e, n in zip(indptr[:-1], indptr[1:], neighbors):
            indices[s: e] = n
        return indptr, indices

    def align_nodes(self, other, parameters, return_as_scipy_csr=True):
        # TODO: Docstring
        ms_utils.LOGGER.info(f"Aligning {self.file_name} with {other.file_name}")
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
        # TODO: Move to seperate function and define type
        ms_utils.LOGGER.info(
            f"Matching ion coordinates from "
            f"{self.file_name} and {other.file_name}"
        )
        @numba.njit(cache=True)
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
                for (
                    max_absolute_error,
                    self_coordinate,
                    other_coordinate
                ) in zip(
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

    @staticmethod
    @numba.njit(cache=True)
    def quick_align(
        self_mzs,
        other_mzs,
        self_mz_order,
        other_mz_order,
        other_rt_order,
        ppm
    ):
        # TODO: Docstring
        max_mz_diff = 1 + ppm * 10**-6
        low_limits = np.searchsorted(
            self_mzs[self_mz_order],
            other_mzs[other_mz_order] / max_mz_diff,
            "left"
        )[other_rt_order]
        high_limits = np.searchsorted(
            self_mzs[self_mz_order],
            other_mzs[other_mz_order] * max_mz_diff,
            "right"
        )[other_rt_order]
        diffs = high_limits - low_limits
        ends = np.cumsum(diffs)
        self_indices = np.empty(ends[-1], np.int64)
        for l, h, e, d in zip(low_limits, high_limits, ends, diffs):
            self_indices[e - d: e] = self_mz_order[l: h]
        other_indices = np.repeat(
            np.arange(len(other_rt_order)),
            high_limits - low_limits
        )
        selection = longest_increasing_subsequence(self_indices)
        self_indices_mask = np.empty(len(selection) + 2, np.int64)
        self_indices_mask[0] = 0
        self_indices_mask[1: -1] = self_indices[selection]
        self_indices_mask[-1] = len(self_mzs) - 1
        other_indices_mask = np.empty(len(selection) + 2, np.int64)
        other_indices_mask[0] = 0
        other_indices_mask[1: -1] = other_indices[selection]
        other_indices_mask[-1] = len(other_mzs) - 1
        return self_indices_mask, other_indices_mask

    def calibrate_precursor_rt(self, other, parameters):
        # TODO: Docstring
        ms_utils.LOGGER.info(
            f"Calibrating PRECURSOR_RT of {self.file_name} with "
            f"{other.file_name}"
        )
        self_mzs = self.get_ion_coordinates("FRAGMENT_MZ")
        other_mzs = other.get_ion_coordinates("FRAGMENT_MZ")
        self_mz_order = np.argsort(self_mzs)
        other_mz_order = np.argsort(other_mzs)
        other_rt_order = np.argsort(other_mz_order)
        # TODO index rts (ordered 0-1, allow e.g. 0.1 distance)?
        self_indices, other_indices = self.quick_align(
            self_mzs,
            other_mzs,
            self_mz_order,
            other_mz_order,
            other_rt_order,
            parameters["calibration_ppm_FRAGMENT_MZ"]
        )
        self_rts = self.get_ion_coordinates("PRECURSOR_RT")
        other_rts = other.get_ion_coordinates(
            "PRECURSOR_RT",
            indices=other_indices
        )
        calibrated_self_rts = []
        for (
            self_start_index,
            self_end_index,
            other_rt_start,
            other_rt_end
        ) in zip(
            self_indices[:-1],
            self_indices[1:],
            other_rts[:-1],
            other_rts[1:]
        ):
            self_rt_start = self_rts[self_start_index]
            self_rt_end = self_rts[self_end_index]
            if self_rt_start == self_rt_end:
                new_rts = np.repeat(
                    other_rt_start,
                    self_end_index - self_start_index
                )
            else:
                slope = (
                    other_rt_end - other_rt_start
                ) / (
                    self_rt_end - self_rt_start
                )
                new_rts = other_rt_start + slope * (
                    self_rts[self_start_index: self_end_index] - self_rt_start
                )
            calibrated_self_rts.append(new_rts)
        calibrated_self_rts.append([other_rts[-1]])
        calibrated_self_rts = np.concatenate(calibrated_self_rts)
        return calibrated_self_rts

    def get_dimensions(
        self,
        remove_fragment_mz=False,
        remove_fragment_logint=False,
        remove_precursor_rt=False
    ):
        """
        Get the dimensions of this ion-network.

        Parameters
        ----------
        remove_fragment_mz : bool
            If True, remove 'FRAGMENT_MZ' from the dimensions.
        remove_fragment_logint : bool
            If True, remove 'FRAGMENT_LOGINT' from the dimensions.
        remove_precursor_rt : bool
            If True, remove 'PRECURSOR_RT' from the dimensions.

        Returns
        -------
        list[str]
            A sorted list with all the dimensions.
        """
        dimensions = set(self.dimensions)
        if remove_fragment_mz:
            dimensions.remove("FRAGMENT_MZ")
        if remove_fragment_logint:
            dimensions.remove("FRAGMENT_LOGINT")
        if remove_precursor_rt:
            dimensions.remove("PRECURSOR_RT")
        return sorted(dimensions)

    def get_precursor_dimensions(
        self,
        remove_precursor_rt=False
    ):
        """
        Get the precursor dimensions of this ion-network.

        Parameters
        ----------
        remove_precursor_rt : bool
            If True, remove 'PRECURSOR_RT' from the dimensions.

        Returns
        -------
        list[str]
            A sorted list with all the precursor dimensions.
        """
        dimensions = set(
            [
                dim for dim in self.dimensions if dim.startswith("PRECURSOR")
            ]
        )
        if remove_precursor_rt:
            dimensions.remove("PRECURSOR_RT")
        return sorted(dimensions)

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
    def dimensions(self):
        return self.get_group_list("nodes")

    @property
    def comment_dimensions(self):
        if "comments" in self.get_group_list():
            dimensions = self.get_group_list("comments")
        else:
            dimensions = []
        return dimensions

    @property
    def node_count(self):
        return self.get_attr("node_count")

    @property
    def edge_count(self):
        return self.get_attr("edge_count")


class Evidence(HDF_MS_Run_File):
    """
    An evidence set containing positive and negative edge evidence, as well as
    node evidence.
    """

    def __init__(
        self,
        reference,
        is_read_only=True,
        new_file=False,
    ):
        # TODO: Docstring
        file_name = self.convert_reference_to_trimmed_file_name(reference)
        super().__init__(
            f"{file_name}.evidence.hdf",
            is_read_only=is_read_only,
            new_file=new_file,
        )
        self.ion_network = Network(reference)
        self.ion_network.evidence = self

    def mutual_collect_evidence_from(
        self,
        other,
        parameters={},
        **kwargs,
    ):
        # TODO: Docstring
        if not self.need_to_set_evidence(other, parameters):
            return
        pairwise_alignment = self.ion_network.align_nodes(
            other.ion_network,
            parameters
        )
        pairwise_alignment_T = pairwise_alignment.T.tocsr()
        if "edges" in kwargs:
            self_edges = kwargs["edges"]
        else:
            self_edges = self.ion_network.get_edges()
        other_edges = other.ion_network.get_edges()
        self_ion_indices, other_ion_indices = pairwise_alignment.nonzero()
        self.align_edges(
            other,
            pairwise_alignment,
            pairwise_alignment_T,
            self_ion_indices,
            self_edges,
            other_edges,
            parameters
        )
        other.align_edges(
            self,
            pairwise_alignment_T,
            pairwise_alignment,
            other_ion_indices,
            other_edges,
            self_edges,
            parameters
        )

    def need_to_set_evidence(self, other, parameters):
        # TODO: Docstring
        if self.is_evidenced_with(other) or other.is_evidenced_with(self):
            if parameters["force_overwrite"]:
                pass
                # self.remove_evidence_from(other)
                # other.remove_evidence_from(self)
            else:
                ms_utils.LOGGER.info(
                    f"Evidence for {self.file_name} and {other.file_name} "
                    f"has already been collected"
                )
                return False
        return True

    def is_evidenced_with(self, other):
        # TODO: Docstring
        return other.run_name in self.get_group_list()

    def align_edges(
        self,
        other,
        pairwise_alignment,
        pairwise_alignment_T,
        aligned_ion_indices,
        self_edges,
        other_edges,
        parameters
    ):
        # TODO: Docstring
        ms_utils.LOGGER.info(
            f"Collecting evidence for {self.ion_network.file_name} from "
            f"{other.ion_network.file_name}"
        )
        parent_group_name = other.ion_network.run_name
        self.create_group(parent_group_name)
        positive = pairwise_alignment * other_edges * pairwise_alignment_T
        positive = (positive + positive.T).multiply(self_edges)
        positive_mask = (self_edges.astype(np.int8) + positive).data == 2
        alignment_mask = np.diff(pairwise_alignment.indptr) != 0
        left_node_indices, right_node_indices = self_edges.nonzero()
        negative_mask = alignment_mask[
            left_node_indices
        ] & alignment_mask[
            right_node_indices
        ] & ~positive_mask
        self.create_dataset(
            "positive_edges",
            positive_mask,
            parent_group_name=parent_group_name,
        )
        self.create_dataset(
            "negative_edges",
            negative_mask,
            parent_group_name=parent_group_name,
        )
        self.create_dataset(
            "aligned_nodes",
            aligned_ion_indices,
            parent_group_name=parent_group_name,
        )
        dimension_overlap = self.ion_network.dimension_overlap(
            other.ion_network
        )
        for parameter_key, parameter_value in parameters.items():
            if parameter_key.startswith("max_alignment_absolute_error_"):
                if parameter_key[24:] not in dimension_overlap:
                    continue
            self.create_attr(
                parameter_key,
                parameter_value,
                parent_group_name=parent_group_name
            )

    def get_alignment(
        self,
        other,
        return_as_scipy_csr=False
    ):
        """
        Get the pairwise alignment between two ion-networks.

        Parameters
        ----------
        other : evidence
            The other evidence set against which to align.
        return_as_scipy_csr : bool
            If True, return the alignment as a scipy.sparse.csr_matrix,
            otherwise, return the alignment as a 2-dimensional np.ndarray.

        Returns
        ----------
        np.ndarray or scipy.sparse.csr_matrix(bool)
            A 2-dimensional array of shape (2, n) with n the number of nodes
            aligned. The first column are the indices of the first ion-network
            and the second column contains the indices of the second
            ion-network. If return_as_scipy_csr is True, a sparse csr matrix
            is created from this array before it is returned.
        """
        self_ions = self.get_aligned_nodes_from_group(
            other.ion_network.run_name,
            return_as_mask=False
        )
        other_ions = other.get_aligned_nodes_from_group(
            self.ion_network.run_name,
            return_as_mask=False
        )
        if return_as_scipy_csr:
            return scipy.sparse.csr_matrix(
                (
                    np.ones(self_ions.shape[0], dtype=np.bool),
                    (self_ions, other_ions)
                ),
                shape=(
                    self.ion_network.node_count,
                    other.ion_network.node_count
                )
            )
        return np.stack([self_ions, other_ions]).T

    def get_aligned_nodes_from_group(self, group_name="", return_as_mask=True):
        # TODO: Docstring
        if group_name == "":
            ions = np.zeros(self.ion_network.node_count, dtype=np.int)
            for group_name in self.get_group_list():
                ion_mask = self.get_dataset(
                    "aligned_nodes",
                    parent_group_name=group_name,
                )
                ions[ion_mask] += 1
        else:
            ions = self.get_dataset(
                "aligned_nodes",
                parent_group_name=group_name,
            )
            if return_as_mask:
                mask = np.repeat(False, self.ion_network.node_count)
                mask[ions] = True
                ions = mask
        return ions

    def get_edge_mask_from_group(self, group_name="", positive=True):
        # TODO: Docstring
        if positive:
            mask_name = "positive_edges"
        else:
            mask_name = "negative_edges"
        if group_name == "":
            edges = np.zeros(self.ion_network.edge_count, dtype=np.int)
            for group_name in self.get_group_list():
                edges += self.get_dataset(
                    mask_name,
                    parent_group_name=group_name,
                )
        else:
            edges = self.get_dataset(
                mask_name,
                parent_group_name=group_name,
            )
        return edges

    @property
    def network_keys(self):
        """
        Get a sorted list with all the ion-network keys providing evidence.

        Returns
        -------
        list[str]
            A sorted list with the keys of all ion-network.
        """
        return self.get_group_list()

    @property
    def evidence_count(self):
        """
        Get the number of ion-network providing evidence.

        Returns
        -------
        int
            The number of ion-network providing evidence.
        """
        return len(self.network_keys)


class Annotation(HDF_MS_Run_File):
    # TODO: Docstring

    def __init__(
        self,
        reference,
        is_read_only=True,
        new_file=False,
    ):
        # TODO: Docstring
        file_name = self.convert_reference_to_trimmed_file_name(reference)
        super().__init__(
            f"{file_name}.annotation.hdf",
            is_read_only=is_read_only,
            new_file=new_file,
        )
        self.evidence = Evidence(reference)
        self.evidence.annotation = self
        self.ion_network = Network(reference)
        self.ion_network.annotation = self

    def create_annotations(self, database, parameters):
        if database is not None:
            if isinstance(database, str):
                database = ms_database.Database(database)
        self.write_parameters(database, parameters)
        self.write_candidates(database, parameters)

    def write_parameters(self, database, parameters):
        # TODO: Docstring
        ms_utils.LOGGER.info(f"Writing parameters to {self.file_name}")
        self.create_attr("database_file_name", database.file_name)
        for key, value in parameters.items():
            self.create_attr(key, value)

    def write_candidates(self, database, parameters):
        # TODO: Docstring
        ms_utils.LOGGER.info(f"Writing candidates to {self.file_name}")
        (
            low_peptide_indices,
            high_peptide_indices
        ) = self.get_candidate_peptide_indices_for_nodes(
            database,
            parameters
        )
        self.create_group("node_candidates")
        self.create_dataset(
            "low_peptide_indices",
            low_peptide_indices,
            parent_group_name="node_candidates"
        )
        self.create_dataset(
            "high_peptide_indices",
            high_peptide_indices,
            parent_group_name="node_candidates"
        )
        (
            edge_indptr,
            edge_indices
        ) = self.get_candidate_peptide_indices_for_edges(
            database,
            low_peptide_indices,
            high_peptide_indices
        )
        self.create_group("edge_candidates")
        self.create_dataset(
            "indptr",
            edge_indptr,
            parent_group_name="edge_candidates"
        )
        self.create_dataset(
            "indices",
            edge_indices,
            parent_group_name="edge_candidates"
        )

    def get_candidate_peptide_indices_for_edges(
        self,
        database,
        low_peptide_indices,
        high_peptide_indices
    ):
        ms_utils.LOGGER.info(f"Writing edge candidates to {self.file_name}")
        database_peptides = database.get_fragment_coordinates("peptide_index")
        indptr, indices = self.ion_network.get_edges(
            indptr_and_indices=True
        )
        max_batch = np.max(
            np.diff(indptr) * (high_peptide_indices - low_peptide_indices)
        )
        # TODO: remove init and just compile upfront?
        self.__get_candidate_peptide_indices_for_edges(
            indptr[:10],
            indices,
            low_peptide_indices,
            high_peptide_indices,
            database_peptides,
            max_batch
        )
        (
            edge_indptr,
            edge_indices
        ) = self.__get_candidate_peptide_indices_for_edges(
            indptr,
            indices,
            low_peptide_indices,
            high_peptide_indices,
            database_peptides,
            max_batch
        )
        return edge_indptr, edge_indices

    @staticmethod
    @numba.njit(cache=True)
    def __get_candidate_peptide_indices_for_edges(
        indptr,
        indices,
        low_peptide_indices,
        high_peptide_indices,
        database_peptides,
        max_batch
    ):
        # TODO: Docstring
        result_indptr = np.empty(indptr[-1], np.int64)
        result_indices = np.empty(max_batch, np.int64)
        current_index = 0
        for start, end, low, high in zip(
            indptr[:-1],
            indptr[1:],
            low_peptide_indices,
            high_peptide_indices,
        ):
            if (low == high) or (start == end):
                result_indptr[start:end] = current_index
                continue
            if (
                (end - start) * (high - low) + current_index
            ) >= result_indices.shape[0]:
                new_result_indices = np.empty(
                    max_batch + result_indices.shape[0],
                    np.int64
                )
                new_result_indices[:result_indices.shape[0]] = result_indices
                result_indices = new_result_indices
            peptide_candidates = database_peptides[low: high]
            peptide_candidates_set = set(peptide_candidates)
            neighbors = indices[start: end]
            for i, neighbor in enumerate(neighbors):
                neighbor_low = low_peptide_indices[neighbor]
                neighbor_high = high_peptide_indices[neighbor]
                if neighbor_low == neighbor_high:
                    result_indptr[start + i] = current_index
                    continue
                neighbor_peptide_candidates = database_peptides[
                    neighbor_low: neighbor_high
                ]
                for neighbor_peptide_candidate in neighbor_peptide_candidates:
                    if neighbor_peptide_candidate in peptide_candidates_set:
                        result_indices[
                            current_index
                        ] = neighbor_peptide_candidate
                        current_index += 1
                result_indptr[start + i] = current_index
        result_indptr[1:] = result_indptr[:-1]
        result_indptr[0] = 0
        return result_indptr, result_indices[:current_index]

    def get_candidate_peptide_indices_for_nodes(
        self,
        database,
        parameters
    ):
        # TODO: Docstring
        ms_utils.LOGGER.info(f"Writing node candidates to {self.file_name}")
        max_ppm = parameters["annotation_ppm"]
        self_mzs = self.ion_network.get_ion_coordinates("FRAGMENT_MZ")
        mz_order = np.argsort(self_mzs)
        database_mzs = database.get_fragment_coordinates("mz")
        low_limits = np.searchsorted(
            np.log(database_mzs) * 10**6,
            np.log(self_mzs[mz_order]) * 10**6 - max_ppm,
            "left"
        )
        high_limits = np.searchsorted(
            np.log(database_mzs) * 10**6,
            np.log(self_mzs[mz_order]) * 10**6 + max_ppm,
            "right"
        )
        inv_order = np.argsort(mz_order)
        return low_limits[inv_order], high_limits[inv_order]


@numba.njit(fastmath=True, cache=True)
def longest_increasing_subsequence(sequence):
    # TODO:Docstring
    M = np.zeros(len(sequence) + 1, np.int64)
    P = np.zeros(len(sequence), np.int64)
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
    return longest_increasing_subsequence
