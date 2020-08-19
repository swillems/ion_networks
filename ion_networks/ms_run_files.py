#!python

# builtin
import os
import ast
import multiprocessing.pool
import csv
# external
import numpy as np
import scipy
import scipy.sparse
import sklearn.linear_model
import numexpr as ne
# local
try:
    from . import ms_database
    from . import ms_utils
    from . import numba_functions
except (ImportError, ModuleNotFoundError):
    import ms_database
    import ms_utils
    import numba_functions


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


class HDF_Network_File(HDF_MS_Run_File):
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
                numba_functions.align_coordinates,
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
            indptr, indices, pointers = numba_functions.make_symmetric(
                indptr,
                indices,
            )
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


class HDF_Evidence_File(HDF_MS_Run_File):
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
        self.ion_network = HDF_Network_File(reference)
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
            f"{key}.evidence.hdf"
        )
        return HDF_Evidence_File(full_file_name)

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
                numba_functions.align_coordinates,
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
            parent_group_name="total"
        )
        self.write_attr(
            "count",
            np.sum(total),
            parent_group_name="total/nodes"
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
                numba_functions.align_edges,
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
        self_peaks = numba_functions.find_peak_indices(
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
        other_peaks = numba_functions.find_peak_indices(
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
                apex, counts = numba_functions.get_unique_apex_and_count(
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
        self_indices, other_indices = numba_functions.quick_align(
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
            parent_group_name="total"
        )
        self.write_attr(
            "count",
            np.sum(total_positive),
            parent_group_name="total/positive_edges"
        )
        total_negative = self.get_edges(positive=False)
        total_negative[negative_pointers] += 1
        self.write_dataset(
            "negative_edges",
            total_negative,
            parent_group_name="total"
        )
        self.write_attr(
            "count",
            np.sum(total_negative),
            parent_group_name="total/negative_edges"
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

    def create_mgf(self, parameters):
        ms_utils.LOGGER.info(f"Creating mgf of {self.file_name}")
        selected_edges = ne.evaluate(
            parameters["edge_threshold"],
            local_dict={
                "positive_edges": self.get_edges(),
                "negative_edges": self.get_edges(positive=False),
                "evidence_run_count": self.run_count
            },
            global_dict={},
        )
        indptr, indices, edge_pointers = self.ion_network.get_edges(
            symmetric=True,
            return_pointers=True
        )
        clusters = numba_functions.cluster_network(
            indptr,
            indices,
            edge_pointers,
            selected_edges,
        )
        cluster_indices = np.argsort(clusters)
        cluster_indptr = np.empty(np.max(clusters + 2), np.int64)
        cluster_indptr[0] = 0
        cluster_indptr[1:] = np.cumsum(np.bincount(clusters))
        mzs, ints = self.ion_network.get_ion_coordinates(
            ["FRAGMENT_MZ", "FRAGMENT_LOGINT"]
        )
        precursor_dimensions = self.ion_network.get_ion_coordinates(
            self.ion_network.precursor_dimensions
        )
        rts = precursor_dimensions[
            self.ion_network.precursor_dimensions.index("PRECURSOR_RT")
        ]
        with open(
            os.path.join(self.directory, f"{self.run_name}.mgf"),
            "w"
        ) as infile:
            for cluster_index in np.flatnonzero(
                np.diff(cluster_indptr) > parameters["minimum_peaks"]
            ):
                infile.write("BEGIN IONS\n")
                cluster = cluster_indices[
                    cluster_indptr[cluster_index]: cluster_indptr[cluster_index + 1]
                ]
                # if expand:
                #     new_indices = [cluster]
                #     for index in cluster:
                #         neighbors = indices[
                #             indptr[index]: indptr[index + 1]
                #         ]
                #         pointers = edge_pointers[
                #             indptr[index]: indptr[index + 1]
                #         ]
                #         selected = selected_edges[pointers]
                #         new_indices.append(neighbors[selected])
                #     cluster = np.unique(np.concatenate(new_indices))
                local_mzs = np.round(mzs[cluster], 4)
                local_ints = np.round(2**ints[cluster], 2)
                local_rts = rts[cluster]
                precursor_strings = []
                for precursor_dimension, array in zip(
                    self.ion_network.precursor_dimensions,
                    precursor_dimensions
                ):
                    precursor_string = np.round(
                        np.average(array[cluster]),
                        2
                    )
                    precursor_strings.append(
                        f"{precursor_dimension}={precursor_string}"
                    )
                infile.write(
                    f"TITLE=cluster_index.{cluster_index}.{cluster_index}. "
                    f"File=\"{self.file_name}\" "
                    f"Precursor_dimensions:\"{' '.join(precursor_strings)}\" "
                    f"NativeID:\"sample=1 period=1 cycle={cluster_index-1} "
                    f"experiment=1\"\n"
                )
                infile.write(
                    f"RTINSECONDS={np.round(np.average(local_rts) * 60, 2)}\n"
                )
                infile.write("PEPMASS=1000\n")
                infile.write("CHARGE=2+\n")
                mz_order = np.argsort(local_mzs)
                for i in mz_order:
                    infile.write(f"{local_mzs[i]} {local_ints[i]}\n")
                infile.write("END IONS\n")

    def annotate(
        self,
        database,
        out_file_name,
        parameters
    ):
        parameters = ms_utils.read_parameters_from_json_file(default="annotation")
        threads = ms_utils.MAX_THREADS
        ms_utils.LOGGER.info(f"Reading {self.ion_network.file_name}")
        inet_mzs = self.ion_network.get_ion_coordinates("FRAGMENT_MZ")
        mz_order = np.argsort(inet_mzs)
        spectra_log_mzs = np.log(inet_mzs[mz_order]) * 10**6
        indptr, indices, edge_pointers = self.ion_network.get_edges(
            symmetric=True,
            return_pointers=True
        )
        ms_utils.LOGGER.info(f"Reading database {database.file_name}")
        peptide_pointers = database.get_fragment_coordinates("peptide_index")
        database_log_mzs = np.log(
            database.get_fragment_coordinates("mz")
        ) * 10**6
        ms_utils.LOGGER.info(
            f"Matching fragments of {self.file_name} with {database.file_name}"
        )
        # TODO mass calibrate
        low_limits = np.searchsorted(
            database_log_mzs,
            spectra_log_mzs - parameters["annotation_ppm"],
            "left"
        )
        high_limits = np.searchsorted(
            database_log_mzs,
            spectra_log_mzs + parameters["annotation_ppm"],
            "right"
        )
        inv_order = np.argsort(mz_order)
        ms_utils.LOGGER.info(
            f"Annotating fragments of {self.file_name} with {database.file_name}"
        )
        selected_edges = ne.evaluate(
            parameters["edge_threshold"],
            local_dict={
                "positive_edges": self.get_edges(),
                "negative_edges": self.get_edges(positive=False),
                "evidence_run_count": self.run_count
            },
            global_dict={},
        )
        with multiprocessing.pool.ThreadPool(threads) as p:
            results = p.starmap(
                numba_functions.annotate_network,
                [
                    (
                        np.arange(i, self.ion_network.node_count, threads),
                        indptr,
                        indices,
                        edge_pointers,
                        selected_edges,
                        low_limits[inv_order],
                        high_limits[inv_order],
                        peptide_pointers,
                    ) for i in range(threads)
                ]
            )
        scores = np.concatenate([r[0] for r in results])
        fragments = np.concatenate([r[1] for r in results])
        ion_indices = np.concatenate([r[2] for r in results])
        count_results = np.concatenate([r[3] for r in results])
        candidate_counts = np.concatenate([r[4] for r in results])
        neighbor_counts = np.concatenate([r[5] for r in results])
        del results
        modified_scores, fdr_values = ms_utils.calculate_modified_score(
            scores,
            count_results,
            neighbor_counts,
            database,
            peptide_pointers[fragments]
        )
        self.export_annotated_csv(
            scores,
            fragments,
            ion_indices,
            count_results,
            candidate_counts,
            neighbor_counts,
            database,
            peptide_pointers,
            out_file_name,
            export_decoys=parameters['export_decoys'],
            fdr_filter=parameters['fdr_filter'],
            fdr_values=fdr_values,
            modified_scores=modified_scores
        )

    def export_annotated_csv(
        self,
        scores,
        fragments,
        ion_indices,
        count_results,
        candidate_counts,
        neighbor_counts,
        database,
        peptide_pointers,
        out_file_name,
        export_decoys,
        fdr_filter,
        fdr_values,
        modified_scores,
    ):
        ms_utils.LOGGER.info(f"Exporting {out_file_name}")
        peptides = peptide_pointers[fragments]
        decoys = database.read_dataset("decoy", "peptides")
        peptide_modifications = database.read_dataset(
            "modifications",
            "peptides"
        )
        peptide_sequences = database.read_dataset("sequence", "peptides")
        fragment_ion_numbers = database.get_fragment_coordinates("ionnumber")
        fragment_is_y_ion = database.get_fragment_coordinates("y_ion")
        self_coordinates = self.ion_network.get_ion_coordinates(
            self.ion_network.dimensions
        )
        with open(out_file_name, "w") as raw_outfile:
            outfile = csv.writer(raw_outfile)
            header = ["Fragment_index"]
            header += self.ion_network.dimensions
            header += [
                "Fragment_ion_type",
                "Fragment_ion_number",
                "Peptide_sequence",
                "Peptide_mods",
                "Peptide_length",
                "Likelihood",
                "Count",
                "Candidates",
                "Neighbors",
                "Modified_score",
                "Decoy"
            ]
            outfile.writerow(header)
            for i, ion_index in enumerate(ion_indices):
                fdr = fdr_values[i]
                if fdr > fdr_filter:
                    continue
                fragment_index = fragments[i]
                peptide_index = peptides[i]
                if (not export_decoys):
                    if decoys[peptide_index]:
                        continue
                row = [ion_index]
                row += [
                    self_coordinate[
                        ion_index
                    ] for self_coordinate in self_coordinates
                ]
                peptide_sequence = peptide_sequences[peptide_index]
                row += [
                    "Y" if fragment_is_y_ion[fragment_index] else "B",
                    fragment_ion_numbers[fragment_index],
                    peptide_sequence,
                    peptide_modifications[peptide_index],
                    len(peptide_sequence),
                    scores[i],
                    count_results[i],
                    candidate_counts[i],
                    neighbor_counts[i],
                    modified_scores[i],
                    decoys[peptide_index],
                ]
                outfile.writerow(row)


class HDF_Annotation_File(HDF_MS_Run_File):
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
        self.evidence = HDF_Evidence_File(reference)
        self.evidence.annotation = self
        self.ion_network = HDF_Network_File(reference)
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
