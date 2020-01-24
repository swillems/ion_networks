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
        logger=logging.getLogger('Log')
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
        if logger is None:
            raise ValueError("No logger was provided")
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
        indptr, indices = create_sparse_edges(
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

    def get_edge_indptr_and_indices(self):
        # TODO
        with h5py.File(self.file_name, "r") as network_file:
            indices = network_file["edges"]["indices"][...]
            indptr = network_file["edges"]["indptr"][...]
        return indptr, indices

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
        self.logger.info(f"Writing edges of ion-network {self.file_name}.")
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
        indices : ellipsis, slice, iterable[int] or iterable[bool]
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
            sliced = False
        except TypeError:
            sliced = True
        with h5py.File(self.file_name, "r") as network_file:
            for dimension in dimensions:
                array = network_file["nodes"][dimension]
                if not sliced:
                    array = array[...]
                array = array[indices]
                arrays.append(array)
        if single_dimension:
            return arrays[0]
        else:
            return arrays

    def align(
        self,
        other,
        alignment_file_name=None,
        parameters=None
    ):
        # TODO
        print("align", self, other, alignment_file_name, parameters)
        pass

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










def create_sparse_edges(
    rts,
    dts,
    max_dt_diff,
    max_rt_diff,
    symmetric
):
    @numba.njit(fastmath=True)
    def numba_wrapper():
        if not symmetric:
            low_limits = np.arange(len(rts)) + 1
        else:
            low_limits = np.searchsorted(
                rts,
                rts - max_rt_diff,
                "left"
            )
        high_limits = np.searchsorted(
            rts,
            rts + max_rt_diff,
            "right"
        )
        neighbors = []
        for dt, l, h in zip(dts, low_limits, high_limits):
            small_dt_diffs = np.abs(dts[l: h] - dt) <= max_dt_diff
            local_neighbors = l + np.flatnonzero(small_dt_diffs)
            neighbors.append(local_neighbors)
        return neighbors
    neighbors = numba_wrapper()
    indptr = np.empty(len(neighbors) + 1, dtype=np.int)
    indptr[0] = 0
    indptr[1:] = np.cumsum([len(n) for n in neighbors])
    indices = np.concatenate(neighbors)
    return indptr, indices
