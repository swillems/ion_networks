#!python

import os

import h5py
import pandas as pd
import numpy as np
import logging
import numba
import time


class Network(object):

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
            The logger that indicates all progress
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
        """
        self.logger.info(
            f"Reading centroided csv file {centroided_csv_file_name}."
        )
        data = pd.read_csv(
            centroided_csv_file_name,
            engine="c",
            dtype=np.float,
        )
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
        rts, dts = self.get_ion_coordinates(["RT", "DT"])
        indptr, indices = self.create_sparse_edges(
            rts,
            dts,
            parameters["max_dt_diff"],
            parameters["max_rt_diff"],
            symmetric=False
        )
        edge_group.create_dataset(
            "indptr",
            data=indptr
        )
        edge_group.create_dataset(
            "indices",
            data=indices
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
        self.creation_time = time.time()
        parameter_group.attrs["creation_time"] = self.creation_time
        for parameter_key, parameter_value in parameters.items():
            parameter_group.attrs[parameter_key] = parameter_value

    def get_ion_coordinates(self, dimensions, indices=...):
        """
        Get an array with ion coordinates from the ion-network.

        Parameters
        ----------
        dimensions : str or list[str]
            The dimension(s) to retrieve from the ion-network.
        indices : ellipsis, slice, int, iterable[int] or iterable[bool]
            The indices that should be selected from the array. This is most
            performant when this is an ellipsis or slice, but fancy indexing
            with a mask, list or np.ndarray are also supported.

        Returns
        -------
        np.ndarray or list[np.ndarray]
            A (list of) numpy array(s) with the ion coordinates from the
            requested dimension(s).
        """
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

    def get_edge_indptr_and_indices(
        self,
        rows=...,
        columns=...,
        return_as_scipy_csr=False
    ):
        # TODO
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
        return indptr, indices

    def create_sparse_edges(
        self,
        rts,
        dts,
        max_dt_diff,
        max_rt_diff,
        symmetric
    ):
        # TODO
        @numba.njit(fastmath=True)
        def numba_wrapper():
            if not symmetric:
                lower_limits = np.arange(len(rts)) + 1
            else:
                lower_limits = np.searchsorted(
                    rts,
                    rts - max_rt_diff,
                    "left"
                )
            upper_limits = np.searchsorted(
                rts,
                rts + max_rt_diff,
                "right"
            )
            neighbors = []
            for dt, low_limit, high_limit in zip(
                dts,
                lower_limits,
                upper_limits
            ):
                small_dt_diffs = np.abs(
                    dts[low_limit: high_limit] - dt
                ) <= max_dt_diff
                local_neighbors = low_limit + np.flatnonzero(small_dt_diffs)
                neighbors.append(local_neighbors)
            return neighbors
        neighbors = numba_wrapper()
        indptr = np.empty(len(neighbors) + 1, dtype=np.int)
        indptr[0] = 0
        indptr[1:] = np.cumsum([len(n) for n in neighbors])
        indices = np.concatenate(neighbors)
        return indptr, indices

    def align(
        self,
        other,
        alignment_file,
        parameters=None
    ):
        # TODO
        # TODO update self.parameters and other.parameters?
        self.logger.info(f"Aligning {self.file_name} with {other.file_name}.")
        first_indices, second_indices = self.pairwise_node_align(
            other,
            parameters["ppm"],
            parameters["max_rt_diff"],
            parameters["max_dt_diff"]
        )
        alignment_file.create_dataset(
            "first_indices",
            data=first_indices
        )
        alignment_file.create_dataset(
            "second_indices",
            data=second_indices
        )

    def pairwise_node_align(self, other, ppm, max_rt_diff, max_dt_diff):
        # TODO
        mz1, rt1, dt1 = self.get_ion_coordinates(["MZ2", "RT", "DT"])
        mz2, rt2, dt2 = other.get_ion_coordinates(["MZ2", "RT", "DT"])
        @numba.njit
        def numba_wrapper():
            max_mz_diff = (1 + ppm * 10**-6)
            mz1_order = np.argsort(mz1)
            mz2_order = np.argsort(mz2)
            low_limits = np.searchsorted(
                mz2[mz2_order],
                mz1[mz1_order] / max_mz_diff,
                "left"
            )
            high_limits = np.searchsorted(
                mz2[mz2_order],
                mz1[mz1_order] * max_mz_diff,
                "right"
            )
            rt2o = rt2[mz2_order]
            dt2o = dt2[mz2_order]
            results = []
            for rt, dt, low, high in zip(
                rt1[mz1_order],
                dt1[mz1_order],
                low_limits,
                high_limits
            ):
                rtd = rt2o[low: high] - rt
                dtd = dt2o[low: high] - dt
                good = (np.abs(rtd) < max_rt_diff) & (np.abs(dtd) < max_dt_diff)
                matches = mz2_order[low + np.flatnonzero(good)]
                results.append(matches)
            return results, mz1_order
        results, mz1_order = numba_wrapper()
        first_indices = np.repeat(mz1_order, [len(l) for l in results])
        second_indices = np.concatenate(results)
        return first_indices, second_indices

    def evidence(
        self,
        alignment_file_name,
        evidence_file_name=None,
        parameters=None
    ):
        # TODO
        print("evidence", self, alignment_file_name, evidence_file_name, parameters)
        pass


if __name__ == "__main__":
    pass
