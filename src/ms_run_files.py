#!python

# builtin
import os
import ast
import multiprocessing.pool
# external
import numpy as np
import scipy
import numba
import scipy.sparse
import sklearn.linear_model
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
        if len(reference_components) < 3:
            raise ValueError("Invalid reference for HDF file")
        if reference_components[-1] not in ["hdf", "csv"]:
            raise ValueError("Invalid reference for HDF file")
        if reference_components[-2] not in ["inet", "evidence", "annotation"]:
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

    @property
    def dimensions(self):
        return ast.literal_eval(self.read_attr("dimensions"))

    @property
    def precursor_dimensions(self):
        return [d for d in self.dimensions if d.startswith("PRECURSOR_")]

    @property
    def ion_comments(self):
        return ast.literal_eval(self.read_attr("ion_comments"))

    @property
    def node_count(self):
        return self.read_attr("node_count")

    @property
    def edge_count(self):
        return self.read_attr("edge_count")

    @property
    def centroids_file_name(self):
        return self.read_attr("centroids_file_name")

    def __init__(
        self,
        reference,
        is_read_only=True,
        new_file=False,
    ):
        # TODO: Docstring
        file_name_prefix = self.convert_reference_to_trimmed_file_name(
            reference
        )
        super().__init__(
            f"{file_name_prefix}.inet.hdf",
            is_read_only=is_read_only,
            new_file=new_file,
        )

    def __write_nodes(self, data):
        """
        Saves all data as individual arrays to an .hdf file. Each array
        will be placed in a 'nodes' group and names according to its column
        name.

        Parameters
        ----------
        data : pd.Dataframe
            A pd.Dataframe with centroided ion peaks.
        """
        ms_utils.LOGGER.info(
            f"Writing {data.shape[0]} nodes of ion-network {self.file_name}..."
        )
        precursor_rts = data["PRECURSOR_RT"].values
        if np.any(precursor_rts[:-1] > precursor_rts[1:]):
            data = data[np.argsort(precursor_rts)]
        regular_columns = [
            column for column in data.columns if not column.startswith("#")
        ]
        self.write_dataset("nodes", data[regular_columns])
        self.write_attr("node_count", data.shape[0])
        self.write_attr("dimensions", sorted(regular_columns))
        comment_columns = [
            column for column in data.columns if column.startswith("#")
        ]
        if len(comment_columns) > 0:
            self.write_dataset("ion_comments", data[comment_columns])
            self.write_attr("ion_comments", sorted(comment_columns))
        else:
            self.write_group("ion_comments")
            self.write_attr("ion_comments", "[]")

    def __create_edges(self, parameters):
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
        thread_count = ms_utils.MAX_THREADS
        precursor_coordinates = np.stack(
            self.get_ion_coordinates(self.precursor_dimensions)
        ).T
        precursor_errors = {
            dimension: parameters["precursor_errors"][
                dimension
            ] for dimension in self.precursor_dimensions
        }
        tuned = parameters["precursor_errors_auto_tune"]
        if tuned:
            precursor_errors.update(
                self.__tune_precursor_errors(
                    parameters["precursor_errors_auto_tune_ion_count"],
                    parameters["precursor_errors_auto_tune_ppm"],
                    parameters["precursor_errors_auto_tune_target_mz"],
                    precursor_errors,
                    parameters["precursor_errors_auto_tune_bandwidth"],
                    parameters["precursor_errors_auto_tune_noise_range"],
                )
            )
        max_errors = np.array(
            [
                precursor_errors[
                    dimension
                ] for dimension in self.precursor_dimensions
            ]
        )
        rt_coordinates = self.get_ion_coordinates("PRECURSOR_RT")
        lower_limits = np.arange(len(rt_coordinates)) + 1
        upper_limits = np.searchsorted(
            rt_coordinates,
            rt_coordinates + parameters["precursor_errors"]["PRECURSOR_RT"],
            "right"
        )
        ms_utils.LOGGER.info(f"Creating edges of ion-network {self.file_name}")
        with multiprocessing.pool.ThreadPool(thread_count) as p:
            results = p.starmap(
                align_coordinates,
                [
                    (
                        np.arange(
                            (i * self.node_count) // thread_count,
                            ((i + 1) * self.node_count) // thread_count
                        ),
                        lower_limits,
                        upper_limits,
                        precursor_coordinates,
                        precursor_coordinates,
                        max_errors,
                    ) for i in range(thread_count)
                ]
            )
        indptr = np.empty(self.node_count + 1, np.int64)
        indptr[0] = 0
        indptr[1:] = np.cumsum(np.concatenate([r[0] for r in results]))
        indices = np.concatenate([r[1] for r in results])
        del results
        self.__write_edges(indptr, indices, precursor_errors, tuned)

    def __write_edges(self, indptr, indices, precursor_errors, tuned):
        ms_utils.LOGGER.info(
            f"Writing {indices.shape[0]} edges of ion-network {self.file_name}"
        )
        self.write_group("edges")
        self.write_dataset(
            "indptr",
            indptr,
            parent_group_name="edges"
        )
        self.write_dataset(
            "indices",
            indices,
            parent_group_name="edges"
        )
        self.write_attr("edge_count", len(indices))
        self.write_attr(
            "precursor_errors",
            precursor_errors
        )
        self.write_attr(
            "precursor_errors_tuned",
            tuned
        )

    def __tune_precursor_errors(
        self,
        to_select_per_sample,
        ppm,
        mz_distance,
        prior_errors,
        precursor_bandwidth_tune,
        usable,
    ):
        ms_utils.LOGGER.info(
            f"Tuning precursor errors of ion-network {self.file_name}"
        )
        estimations = {}
        selected_indices = np.argpartition(
            self.get_ion_coordinates("FRAGMENT_LOGINT"),
            - to_select_per_sample
        )[-to_select_per_sample:]
        fragment_mzs = self.get_ion_coordinates("FRAGMENT_MZ")
        selected_indices = selected_indices[
            np.argsort(fragment_mzs[selected_indices])
        ]
        fragment_mzs = fragment_mzs[selected_indices]
        lower_limits = np.searchsorted(
            fragment_mzs,
            (fragment_mzs + mz_distance) / (1 + ppm * 10**-6),
            "left"
        )
        upper_limits = np.searchsorted(
            fragment_mzs,
            (fragment_mzs + mz_distance) * (1 + ppm * 10**-6),
            "right"
        )
        indptr = np.zeros(self.node_count + 1, np.int64)
        indptr[selected_indices + 1] = upper_limits - lower_limits
        indptr = np.cumsum(indptr)
        order = np.argsort(selected_indices)
        indices = np.concatenate(
            [
                selected_indices[low: high] for low, high in zip(
                    lower_limits[order],
                    upper_limits[order]
                )
            ]
        )
        use_slice = slice(
            int(indices.shape[0] * usable[0]),
            int(indices.shape[0] * usable[1])
        )
        for dimension in self.precursor_dimensions:
            precursor_array = self.get_ion_coordinates(dimension)
            precursor_differences = np.sort(
                np.abs(
                    np.repeat(
                        precursor_array,
                        np.diff(indptr)
                    ) - precursor_array[indices]
                )
            )
            bandwidth = prior_errors[dimension] / precursor_bandwidth_tune
            frequency = np.searchsorted(
                precursor_differences,
                precursor_differences + bandwidth
            ) - np.searchsorted(
                precursor_differences,
                precursor_differences - bandwidth
            )
            max_frequency_index = np.argmax(frequency)
            if max_frequency_index > indices.shape[0] * usable[0]:
                ms_utils.LOGGER.info(
                    f"{dimension} error tuning of ion-network "
                    f"{self.file_name} failed"
                )
                continue
            try:
                ransac_regressor = sklearn.linear_model.RANSACRegressor(
                    random_state=1
                )
                ransac_regressor.fit(
                    precursor_differences[use_slice].reshape(-1, 1),
                    frequency[use_slice].reshape(-1, 1),
                )
            except ValueError:
                ms_utils.LOGGER.info(
                    f"{dimension} error tuning of ion-network "
                    f"{self.file_name} failed"
                )
                continue
            min_index = max_frequency_index + np.argmin(
                ransac_regressor.predict(
                    precursor_differences[max_frequency_index:].reshape(-1, 1)
                ).flatten() < frequency[max_frequency_index:]
            )
            if max_frequency_index > min_index:
                ms_utils.LOGGER.info(
                    f"{dimension} error tuning of ion-network "
                    f"{self.file_name} failed"
                )
                continue
            tuning = precursor_differences[min_index]
            if (tuning > prior_errors[dimension]) or tuning < bandwidth:
                ms_utils.LOGGER.info(
                    f"{dimension} error tuning of ion-network "
                    f"{self.file_name} was out of bounds"
                )
                continue
            estimations[dimension] = tuning
            ms_utils.LOGGER.info(
                f"{dimension} error tuning of ion-network "
                f"{self.file_name}: {tuning}"
            )
        return estimations

    def create(
        self,
        centroids_file_name,
        parameters
    ):
        """
        Creates a .hdf file with nodes and edges from a pd.DataFrame.

        Parameters
        ----------
        dataframe : pd.Dataframe
            A pd.Dataframe with centroided fragment ion coordinates.
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        """
        ms_utils.LOGGER.info(
            f"Creating ion-network {self.file_name}..."
        )
        dataframe = ms_utils.read_centroided_csv_file(
            centroids_file_name,
            parameters,
        )
        self.write_attr("centroids_file_name", centroids_file_name)
        self.__write_nodes(dataframe)
        self.__create_edges(parameters)

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
        arrays = [
            self.read_dataset(
                dimension,
                parent_group_name="nodes",
                indices=indices
            ) for dimension in dimensions
        ]
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
        arrays = [
            self.read_dataset(
                dimension,
                parent_group_name="comments",
                indices=indices
            ) for dimension in dimensions
        ]
        if single_dimension:
            return arrays[0]
        else:
            return arrays

    def get_edges(
        self,
        symmetric=False,
        return_pointers=False,
        return_as_scipy_csr=False,
        return_as_pairs=False,
    ):
        # TODO: Docstring
        indptr = self.read_dataset("indptr", "edges")
        indices = self.read_dataset("indices", "edges")
        if symmetric:
            indptr, indices, pointers = make_symmetric(indptr, indices)
        elif return_pointers:
            pointers = np.arange(self.edge_count)
        if return_as_scipy_csr:
            if not return_pointers:
                pointers = np.ones(self.edge_count, dtype=np.bool)
            return scipy.sparse.csr_matrix(
                (
                    pointers,
                    indices,
                    indptr
                ),
                shape=(self.node_count, self.node_count)
            )
        if return_as_pairs:
            indptr = np.repeat(
                np.arange(self.node_count),
                np.diff(indptr)
            )
        if return_pointers:
            return indptr, indices, pointers
        else:
            return indptr, indices


class Evidence(HDF_MS_Run_File):
    """
    An evidence set containing positive and negative edge evidence, as well as
    node evidence.
    """

    @property
    def runs(self):
        return self.read_group("runs")

    @property
    def run_count(self):
        return self.read_attr("run_count")

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
        if new_file:
            self.write_group("runs")
            self.write_group("total")
            self.write_dataset(
                "positive_edges",
                np.zeros(self.ion_network.edge_count, np.int16),
                parent_group_name="total"
            )
            self.write_dataset(
                "negative_edges",
                np.zeros(self.ion_network.edge_count, np.int16),
                parent_group_name="total"
            )
            self.write_dataset(
                "nodes",
                np.zeros(self.ion_network.node_count, np.int16),
                parent_group_name="total"
            )
            self.write_attr("run_count", 0)

    def __contains__(self, other):
        if isinstance(other, str):
            return other in self.runs
        else:
            return other.run_name in self.runs

    def __getitem__(self, key):
        if key not in self:
            raise KeyError("Evidence {key} is missing from {self}")
        full_file_name = os.path.join(
            self.directory,
            key + ".evidence.hdf"
        )
        return Evidence(full_file_name)

    def __iter__(self):
        self.__iterator_pointer = 0
        return self

    def __next__(self):
        if self.__iterator_pointer < self.run_count:
            run_name = self.runs[self.__iterator_pointer]
            self.__iterator_pointer += 1
            return self[run_name]
        else:
            raise StopIteration

    def __align_nodes(self, other, parameters):
        # TODO: Docstring
        thread_count = ms_utils.MAX_THREADS
        errors = parameters["fragment_errors"]
        dimensions = sorted(
            set(self.ion_network.dimensions) & set(other.ion_network.dimensions)
        )
        dimensions.remove("FRAGMENT_LOGINT")
        self_coordinates = np.stack(
            self.ion_network.get_ion_coordinates(dimensions)
        ).T
        other_coordinates = np.stack(
            other.ion_network.get_ion_coordinates(dimensions)
        ).T
        fragment_index = dimensions.index("FRAGMENT_MZ")
        self_mz_order = np.argsort(self_coordinates[:, fragment_index])
        other_mz_order = np.argsort(other_coordinates[:, fragment_index])
        self_coordinates = self_coordinates[self_mz_order]
        other_coordinates = other_coordinates[other_mz_order]
        if parameters["calibrate_FRAGMENT_MZ"]:
            other_mzs, fragment_ppm_correction = other.__calibrate_fragment_mz(
                self,
                other_coordinates[:, fragment_index],
                self_coordinates[:, fragment_index],
                parameters
            )
            other_coordinates[:, fragment_index] = other_mzs
        else:
            fragment_ppm_correction = 0
        if parameters["calibrate_PRECURSOR_RT"]:
            other_rts = other.__calibrate_precursor_rt(
                self,
                fragment_ppm_correction,
                parameters
            )[other_mz_order]
            other_coordinates[:, dimensions.index("PRECURSOR_RT")] = other_rts
        self_coordinates[:, fragment_index] = np.log(
            self_coordinates[:, fragment_index]
        ) * 10**6
        other_coordinates[:, fragment_index] = np.log(
            other_coordinates[:, fragment_index]
        ) * 10**6
        lower_limits = np.searchsorted(
            other_coordinates[:, fragment_index],
            self_coordinates[:, fragment_index] - errors["FRAGMENT_MZ"],
            "left"
        )
        upper_limits = np.searchsorted(
            other_coordinates[:, fragment_index],
            self_coordinates[:, fragment_index] + errors["FRAGMENT_MZ"],
            "right"
        )
        max_errors = np.array([errors[dimension] for dimension in dimensions])
        ms_utils.LOGGER.info(
            f"Aligning nodes of {self.file_name} and {other.file_name}"
        )
        with multiprocessing.pool.ThreadPool(thread_count) as p:
            results = p.starmap(
                align_coordinates,
                [
                    (
                        np.arange(
                            (i * self.ion_network.node_count) // thread_count,
                            ((i + 1) * self.ion_network.node_count) // thread_count
                        ),
                        lower_limits,
                        upper_limits,
                        self_coordinates,
                        other_coordinates,
                        max_errors,
                    ) for i in range(thread_count)
                ]
            )
        indptr = np.empty(self.ion_network.node_count + 1, np.int64)
        indptr[0] = 0
        indptr[1:] = np.cumsum(np.concatenate([r[0] for r in results]))
        self_indices = np.repeat(self_mz_order, np.diff(indptr))
        other_indices = other_mz_order[np.concatenate([r[1] for r in results])]
        del results
        fragment_errors = {
            dimension: errors[dimension] for dimension in dimensions
        }
        if parameters["fragment_errors_auto_tune_minimal_improvement"] > -1:
            selection, fragment_errors = self.__auto_tune_fragment_errors(
                other,
                self_indices,
                other_indices,
                dimensions,
                self_coordinates[np.argsort(self_mz_order)],
                other_coordinates[np.argsort(other_mz_order)],
                parameters["fragment_errors_auto_tune_minimal_improvement"],
                fragment_errors
            )
            self_indices = self_indices[selection]
            other_indices = other_indices[selection]
        self_unique = np.bincount(self_indices) == 1
        other_unique = np.bincount(other_indices) == 1
        unique = self_unique[self_indices] & other_unique[other_indices]
        self_intensities = self.ion_network.get_ion_coordinates(
            "FRAGMENT_LOGINT",
            self_indices[unique]
        )
        other_median_intensity = other.ion_network.get_ion_coordinates(
            "FRAGMENT_LOGINT",
            other_indices[unique]
        )
        intensity_correction = np.median(
            self_intensities - other_median_intensity
        )
        self.__write_nodes(
            other,
            self_indices,
            unique,
            fragment_errors,
            intensity_correction,
            fragment_ppm_correction
        )
        other.__write_nodes(
            self,
            other_indices,
            unique,
            fragment_errors,
            -intensity_correction,
            -fragment_ppm_correction
        )

    def __write_nodes(
        self,
        other,
        self_indices,
        unique,
        fragment_errors,
        intensity_correction,
        fragment_ppm_correction
    ):
        ms_utils.LOGGER.info(
            f"Writing {np.sum(unique)} unique nodes and {np.sum(~unique)} "
            f"ambiguous nodes from {self.file_name} with {other.file_name}"
        )
        self.write_group("nodes", parent_group_name=f"runs/{other.run_name}")
        self.write_dataset(
            "unique",
            self_indices[unique],
            parent_group_name=f"runs/{other.run_name}/nodes"
        )
        self.write_dataset(
            "ambiguous",
            self_indices[~unique],
            parent_group_name=f"runs/{other.run_name}/nodes"
        )
        self.write_attr(
            "fragment_errors",
            fragment_errors,
            parent_group_name=f"runs/{other.run_name}/nodes"
        )
        self.write_attr(
            "unique_count",
            np.sum(unique),
            parent_group_name=f"runs/{other.run_name}/nodes"
        )
        self.write_attr(
            "ambiguous_count",
            np.sum(~unique),
            parent_group_name=f"runs/{other.run_name}/nodes"
        )
        total_mask = np.ones(self.ion_network.node_count, np.bool_)
        total_mask[self_indices] = False
        self.write_attr(
            "unaligned_count",
            np.sum(total_mask),
            parent_group_name=f"runs/{other.run_name}/nodes"
        )
        self.write_attr(
            "median_intensity_correction",
            intensity_correction,
            parent_group_name=f"runs/{other.run_name}/nodes"
        )
        self.write_attr(
            "fragment_ppm_correction",
            intensity_correction,
            parent_group_name=f"runs/{other.run_name}/nodes"
        )
        total = self.get_nodes()
        total[self_indices[unique]] += 1
        self.write_dataset(
            "nodes",
            total,
            parent_group_name=f"total"
        )
        self.write_attr(
            "count",
            np.sum(total),
            parent_group_name=f"total/nodes"
        )

    def __align_edges(
        self,
        other,
        parameters,
        self_indptr_=None,
        self_indices_=None,
        self_pointers_=None,
    ):
        ms_utils.LOGGER.info(
            f"Aligning edges of {self.file_name} and {other.file_name}"
        )
        if (self_indptr_ is None) or (self_indices_ is None) or (self_pointers_ is None):
            self_indptr_, self_indices_, self_pointers_ = self.ion_network.get_edges(
                symmetric=True,
                return_pointers=True
            )
        thread_count = ms_utils.MAX_THREADS
        other_indptr, other_indices = other.ion_network.get_edges()
        self_alignment = self.get_nodes(other)
        other_alignment = other.get_nodes(self)
        self_alignment_mask = np.zeros(self.ion_network.node_count, np.bool_)
        self_alignment_mask[self_alignment] = True
        alignment = np.empty_like(self_indptr_)
        alignment[self_alignment] = other_alignment
        with multiprocessing.pool.ThreadPool(thread_count) as p:
            results = p.starmap(
                align_edges,
                [
                    (
                        self_alignment[i::thread_count],
                        self_indptr_,
                        self_indices_,
                        self_pointers_,
                        other_indptr,
                        other_indices,
                        alignment,
                        self_alignment_mask,
                    ) for i in range(thread_count)
                ]
            )
        self_positive_pointers = np.concatenate([r[0] for r in results])
        other_positive_pointers = np.concatenate([r[1] for r in results])
        del results
        self_indptr, self_indices = self.ion_network.get_edges()
        self_negative_mask = np.repeat(
            self_alignment_mask,
            np.diff(self_indptr)
        )
        self_negative_mask &= self_alignment_mask[self_indices]
        self_negative_mask[self_positive_pointers] = False
        other_alignment_mask = np.zeros(other.ion_network.node_count, np.bool_)
        other_alignment_mask[other_alignment] = True
        other_negative_mask = np.repeat(
            other_alignment_mask,
            np.diff(other_indptr)
        )
        other_negative_mask &= other_alignment_mask[other_indices]
        other_negative_mask[other_positive_pointers] = False
        self.__write_edges(
            other,
            self_positive_pointers,
            np.flatnonzero(self_negative_mask),
        )
        other.__write_edges(
            self,
            other_positive_pointers,
            np.flatnonzero(other_negative_mask),
        )

    def __calibrate_fragment_mz(self, other, self_mzs, other_mzs, parameters):
        ms_utils.LOGGER.info(
            f"Calibrating PRECURSOR_MZ of {self.file_name} with "
            f"{other.file_name}"
        )
        ppm = parameters["calibration_ppm_FRAGMENT_MZ"]
        mass_defect = parameters["calibration_mass_defect"]
        self_lower = np.searchsorted(
            self_mzs,
            self_mzs * (1 - ppm / 10**6),
            "left"
        )
        self_upper = np.searchsorted(
            self_mzs,
            self_mzs * (1 + ppm / 10**6),
            "right"
        )
        self_peaks = find_peak_indices(
            self_mzs,
            self_upper - self_lower,
            mass_defect
        )
        other_lower = np.searchsorted(
            other_mzs,
            other_mzs * (1 - ppm / 10**6),
            "left"
        )
        other_upper = np.searchsorted(
            other_mzs,
            other_mzs * (1 + ppm / 10**6),
            "right"
        )
        other_peaks = find_peak_indices(
            other_mzs,
            other_upper - other_lower,
            mass_defect
        )
        if self_peaks.shape[0] > other_peaks.shape[0]:
            self_peaks = self_peaks[:other_peaks.shape[0]]
        elif other_peaks.shape[0] > self_peaks.shape[0]:
            other_peaks = other_peaks[:self_peaks.shape[0]]
        fragment_ppm_correction = np.median(
            (self_mzs[self_peaks] - other_mzs[other_peaks]) * 10**6 / self_mzs[self_peaks]
        )
        return (
            self_mzs * (1 - fragment_ppm_correction / 10**6),
            fragment_ppm_correction
        )

    def __auto_tune_fragment_errors(
        self,
        other,
        self_indices,
        other_indices,
        dimensions,
        self_coordinates,
        other_coordinates,
        minimal_improvement,
        fragment_errors,
    ):
        ms_utils.LOGGER.info(
            f"Auto-tuning fragment errors of {self.file_name} and {other.file_name}"
        )
        selection = np.ones(self_indices.shape[0], np.bool_)
        best_count = 0
        previous_best_count = 0
        best_dimension = ""
        arrays = []
        orders = []
        for i, dimension in enumerate(dimensions):
            array = np.abs(
                self_coordinates[
                    self_indices,
                    i
                ] - other_coordinates[
                    other_indices,
                    i
                ]
            )
            arrays.append(array)
            orders.append(np.argsort(array))
        while True:
            bad_selection = np.array([], np.int64)
            for i, dimension in enumerate(dimensions):
                if dimension == best_dimension:
                    continue
                array = arrays[i][selection]
                order = orders[i][selection]
                apex, counts = get_unique_apex_and_count(
                    self_indices[order],
                    other_indices[order],
                    return_all_counts=False
                )
                if counts[0] > best_count:
                    best_count = counts[0]
                    best_dimension = dimension
                    bad_selection = order[apex:]
                    best_value = array[order[apex]]
            if best_count > previous_best_count * (1 + minimal_improvement):
                selection[bad_selection] = False
                previous_best_count = best_count
                fragment_errors[best_dimension] = best_value
                ms_utils.LOGGER.info(
                    f"Updating {dimension} fragment errors between "
                    f"{self.file_name} and {other.file_name}: {best_value}"
                )
            else:
                break
        return selection, fragment_errors

    def __calibrate_precursor_rt(
        self,
        other,
        fragment_ppm_correction,
        parameters
    ):
        # TODO: Docstring
        ms_utils.LOGGER.info(
            f"Calibrating PRECURSOR_RT of {self.file_name} with "
            f"{other.file_name}"
        )
        self_mzs = self.ion_network.get_ion_coordinates("FRAGMENT_MZ")
        other_mzs = other.ion_network.get_ion_coordinates("FRAGMENT_MZ") * (
            1 + fragment_ppm_correction / 10**6
        )
        self_mz_order = np.argsort(self_mzs)
        other_mz_order = np.argsort(other_mzs)
        other_rt_order = np.argsort(other_mz_order)
        # TODO index rts (ordered 0-1, allow e.g. 0.1 distance)?
        self_indices, other_indices = quick_align(
            self_mzs,
            other_mzs,
            self_mz_order,
            other_mz_order,
            other_rt_order,
            parameters["calibration_ppm_FRAGMENT_MZ"]
        )
        self_rts = self.ion_network.get_ion_coordinates("PRECURSOR_RT")
        other_rts = other.ion_network.get_ion_coordinates(
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

    def __write_edges(
        self,
        other,
        positive_pointers,
        negative_pointers,
    ):
        ms_utils.LOGGER.info(
            f"Writing {positive_pointers.shape[0]} positive edges "
            f"and {negative_pointers.shape[0]} negative edges from "
            f"{self.file_name} with {other.file_name}"
        )
        self.write_group("edges", parent_group_name=f"runs/{other.run_name}")
        self.write_dataset(
            "positive",
            positive_pointers,
            parent_group_name=f"runs/{other.run_name}/edges"
        )
        self.write_dataset(
            "negative",
            negative_pointers,
            parent_group_name=f"runs/{other.run_name}/edges"
        )
        positive_count = positive_pointers.shape[0]
        self.write_attr(
            "positive_count",
            positive_count,
            parent_group_name=f"runs/{other.run_name}/edges"
        )
        negative_count = negative_pointers.shape[0]
        self.write_attr(
            "negative_count",
            negative_count,
            parent_group_name=f"runs/{other.run_name}/edges"
        )
        self.write_attr(
            "unaligned_count",
            self.ion_network.edge_count - negative_count - positive_count,
            parent_group_name=f"runs/{other.run_name}/edges"
        )
        total_positive = self.get_edges()
        total_positive[positive_pointers] += 1
        self.write_dataset(
            "positive_edges",
            total_positive,
            parent_group_name=f"total"
        )
        self.write_attr(
            "count",
            np.sum(total_positive),
            parent_group_name=f"total/positive_edges"
        )
        total_negative = self.get_edges(positive=False)
        total_negative[negative_pointers] += 1
        self.write_dataset(
            "negative_edges",
            total_negative,
            parent_group_name=f"total"
        )
        self.write_attr(
            "count",
            np.sum(total_negative),
            parent_group_name=f"total/negative_edges"
        )

    def align(
        self,
        other,
        parameters,
        indptr=None,
        indices=None,
        pointers=None,
        **kwargs,
    ):
        # TODO: Docstring
        if not parameters["force_overwrite"]:
            if (self in other) and (other in self):
                ms_utils.LOGGER.info(
                    f"{self.ion_network.file_name} and "
                    f"{other.ion_network.file_name} have already been aligned"
                )
                return
        ms_utils.LOGGER.info(
            f"Aligning {self.ion_network.file_name} and "
            f"{other.ion_network.file_name}"
        )
        self.write_group(other.run_name, parent_group_name="runs")
        other.write_group(self.run_name, parent_group_name="runs")
        self.write_attr("run_count", self.run_count + 1)
        other.write_attr("run_count", other.run_count + 1)
        self.__align_nodes(other, parameters)
        self.__align_edges(
            other,
            parameters,
            self_indptr_=indptr,
            self_indices_=indices,
            self_pointers_=pointers,
        )

    def get_nodes(
        self,
        other=None,
        return_as_mask=False,
        kind="unique"
    ):
        # TODO: Docstring
        if other is None:
            return self.read_dataset(
                "nodes",
                parent_group_name="total",
            )
        else:
            ions = self.read_dataset(
                kind,
                parent_group_name=f"runs/{other.run_name}/nodes",
            )
            if return_as_mask:
                mask = np.empty(self.ion_network.node_count, np.bool_)
                mask[ions] = True
                ions = mask
            return ions

    def get_edges(
        self,
        other=None,
        return_as_mask=False,
        positive=True
    ):
        # TODO: Docstring
        if other is None:
            return self.read_dataset(
                "positive_edges" if positive else "negative_edges",
                parent_group_name="total",
            )
        else:
            edges = self.read_dataset(
                "positive" if positive else "negative",
                parent_group_name="runs/{other.run_name}/edges",
            )
            if return_as_mask:
                mask = np.empty(self.ion_network.edge_count, np.bool_)
                mask[edges] = True
                edges = mask
            return edges


class Annotation(HDF_MS_Run_File):
    # TODO: Docstring

    def __init__(
        self,
        reference,
        is_read_only=True,
        new_file=False,
    ):
        # TODO: Docstring
        file_name = self.convert_reference_to_trimmed_file_name(
            reference
        )
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
        self.write_attr("database_file_name", database.file_name)
        for key, value in parameters.items():
            self.write_attr(key, value)

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
        self.write_group("node_candidates")
        self.write_dataset(
            "low_peptide_indices",
            low_peptide_indices,
            parent_group_name="node_candidates"
        )
        self.write_dataset(
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
        self.write_group("edge_candidates")
        self.write_dataset(
            "indptr",
            edge_indptr,
            parent_group_name="edge_candidates"
        )
        self.write_dataset(
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
    @numba.njit()
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
        ms_utils.LOGGER.info(
            f"Writing node candidates to {self.file_name}"
        )
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


@numba.njit(nogil=True, cache=True)
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


@numba.njit(nogil=True, cache=True)
def increase_buffer(buffer, max_batch=10**7):
    new_buffer = np.empty(buffer.shape[0] + max_batch, np.int64)
    new_buffer[:len(buffer)] = buffer
    return new_buffer


@numba.njit(nogil=True, cache=True)
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
    selection = longest_increasing_subsequence(self_indices)
    self_indices_mask = np.empty(len(selection) + 2, np.int64)
    self_indices_mask[0] = 0
    self_indices_mask[1: -1] = self_indices[selection]
    del self_indices
    self_indices_mask[-1] = len(self_mzs) - 1
    other_indices_mask = np.empty(len(selection) + 2, np.int64)
    other_indices_mask[0] = 0
    other_indices = np.repeat(
        np.arange(len(other_rt_order)),
        high_limits - low_limits
    )
    other_indices_mask[1: -1] = other_indices[selection]
    other_indices_mask[-1] = len(other_mzs) - 1
    return self_indices_mask, other_indices_mask


@numba.njit(nogil=True, cache=True)
def align_coordinates(
    queries,
    lower_limits,
    upper_limits,
    self_coordinates,
    other_coordinates,
    max_errors,
    # kind="euclidean"
):
    indptr = np.zeros(len(queries), np.int64)
    indices = np.empty(10**7, np.int64)
    total = 0
    for index, query in enumerate(queries):
        low_limit = lower_limits[query]
        high_limit = upper_limits[query]
        candidate_count = high_limit - low_limit
        if candidate_count == 0:
            continue
        elif (candidate_count + total) >= len(indices):
            indices = increase_buffer(indices)
        dists = other_coordinates[low_limit: high_limit] - self_coordinates[query]
        # TODO: what if error==0?
        # if kind == "euclidean":
        dists /= max_errors
        dists = dists**2
        projected_dists = np.sum(dists, axis=1)
        projected_dists = np.sqrt(projected_dists)
        candidates = low_limit + np.flatnonzero(projected_dists <= 1)
        # elif kind == "manhattan":
        #     projected_dists = np.all(dists < max_errors, axis=1)
        #     candidates = low_limit + np.flatnonzero(projected_dists)
        candidate_count = len(candidates)
        indices[total: total + candidate_count] = candidates
        indptr[index] = candidate_count
        total += candidate_count
    return (indptr, indices[:total])


@numba.njit(nogil=True, cache=True)
def make_symmetric(indptr, indices):
    # TODO: multithread?
    offsets = np.cumsum(np.bincount(indices))
    indptr_ = indptr.copy()
    indptr_[1:1 + offsets.shape[0]] += offsets
    indptr_[1 + offsets.shape[0]:] += offsets[-1]
    indices_ = np.empty(indptr_[-1], np.int64)
    pointers_ = np.empty_like(indices_)
    offsets = indptr_[:-1] + np.diff(indptr)
    for index in range(indptr.shape[0] - 1):
        start = indptr[index]
        end = indptr[index + 1]
        current_indices = indices[start: end]
        pointers = np.arange(start, end)
        start_ = indptr_[index]
        end_ = start_ + current_indices.shape[0]
        indices_[start_: end_] = current_indices
        pointers_[start_: end_] = pointers
        current_offsets = offsets[current_indices]
        indices_[current_offsets] = index
        pointers_[current_offsets] = pointers
        offsets[current_indices] += 1
    return indptr_, indices_, pointers_


@numba.njit(nogil=True, cache=True)
def align_edges(
    queries,
    self_indptr,
    self_indices,
    self_pointers,
    other_indptr,
    other_indices,
    alignment,
    alignment_mask,
):
    self_pointers_ = np.empty(10**7, np.int64)
    other_pointers_ = np.empty(10**7, np.int64)
    pointer_offset = 0
    for index in queries:
        possible_start = self_indptr[index]
        possible_end = self_indptr[index + 1]
        if possible_start == possible_end:
            continue
        current_index = alignment[index]
        current_start = other_indptr[current_index]
        current_end = other_indptr[current_index + 1]
        if current_start == current_end:
            continue
        possible_indices = self_indices[possible_start: possible_end]
        possible_mask = alignment_mask[possible_indices]
        if not np.any(possible_mask):
            continue
        possible_indices = alignment[possible_indices[possible_mask]]
        possible_pointers = self_pointers[possible_start: possible_end][
            possible_mask
        ]
        current_indices = other_indices[current_start: current_end]
        candidates1 = np.searchsorted(
            current_indices,
            possible_indices,
            "left"
        )
        candidates2 = np.searchsorted(
            current_indices,
            possible_indices,
            "right"
        )
        overlap = np.flatnonzero(candidates2 != candidates1)
        overlap_count = len(overlap)
        if len(overlap) == 0:
            continue
        elif (overlap_count + pointer_offset) >= len(self_pointers_):
            self_pointers_ = increase_buffer(self_pointers_)
            other_pointers_ = increase_buffer(other_pointers_)
        self_pointers_[
            pointer_offset: pointer_offset + overlap_count
        ] = possible_pointers[overlap]
        other_pointers_[
            pointer_offset: pointer_offset + overlap_count
        ] = current_start + candidates1[overlap]
        pointer_offset += overlap_count
    return self_pointers_[:pointer_offset], other_pointers_[:pointer_offset]


@numba.njit(cache=True)
def find_peak_indices(input_array, output_array, max_distance):
    peaks = np.zeros(int(input_array[-1]), np.int64)
    current_max_mz = 0
    current_max_int = 0
    current_max_index = 0
    for index, (intensity, mz) in enumerate(zip(output_array, input_array)):
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


@numba.njit(nogil=True, cache=True)
def get_unique_apex_and_count(
    ordered_self_indices,
    ordered_other_indices,
    return_all_counts=True
):
    counts = np.zeros_like(ordered_self_indices)
    self_max = np.max(ordered_self_indices)
    other_max = np.max(ordered_other_indices)
    unique_pair = np.zeros(counts.shape[0], np.bool_)
    self_frequencies = np.zeros(self_max + 1, np.int64)
    other_frequencies = np.zeros(other_max + 1, np.int64)
    self_indptr = np.empty(self_max + 2, np.int64)
    self_indptr[0] = 0
    self_indptr[1:] = np.cumsum(np.bincount(ordered_self_indices))
    self_order = np.argsort(ordered_self_indices)
    other_indptr = np.empty(other_max + 2, np.int64)
    other_indptr[0] = 0
    other_indptr[1:] = np.cumsum(np.bincount(ordered_other_indices))
    other_order = np.argsort(ordered_other_indices)
    unique_count = 0
    max_count = 0
    apex = 0
    for i in range(counts.shape[0]):
        self_index = ordered_self_indices[i]
        other_index = ordered_other_indices[i]
        if (
            self_frequencies[self_index] == 0
        ) & (
            other_frequencies[other_index] == 0
        ):
            unique_count += 1
            unique_pair[i] = True
            if unique_count > max_count:
                apex = i
                max_count = unique_count
        else:
            self_locs = self_order[
                self_indptr[self_index]: self_indptr[self_index + 1]
            ]
            if np.any(unique_pair[self_locs]):
                unique_count -= 1
            other_locs = other_order[
                other_indptr[other_index]: other_indptr[other_index + 1]
            ]
            if np.any(unique_pair[other_locs]):
                unique_count -= 1
            unique_pair[self_locs] = False
            unique_pair[other_locs] = False
        self_frequencies[self_index] += 1
        other_frequencies[other_index] += 1
        counts[i] = unique_count
    if not return_all_counts:
        counts = counts[apex: apex + 1]
    return apex, counts
