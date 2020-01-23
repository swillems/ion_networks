#!python

import os

import h5py
import scipy
import pandas as pd
import numpy as np
import multiprocessing as mp
import logging

import numba

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
            node_group = network_file.create_group("nodes")
            edge_group = network_file.create_group("edges")
            self.write_nodes(node_group, data, parameters)
            self.write_edges(edge_group, parameters)
            self.set_hash_id(parameters)

    def write_nodes(self, node_group, data, parameters):
        """
        Saves all data as individual datasets in the node_group of an .hdf file.

        Parameters
        ----------
        node_group : h5py.group
            The location within an opened .hdf file where to save the data.
        data : pd.Dataframe
            A pd.Dataframe with centroided ion peaks.
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        """
        self.logger.info(f"Writing nodes of ion-network {self.file_name}.")
        for column in data.columns:
            node_group.create_dataset(
                column,
                data=data[column]
            )

    def write_edges(self, edge_group, parameters):
        # TODO
        self.logger.info(f"Writing edges of ion-network {self.file_name}.")
        rts, dts = self.get_ion_coordinates(["RT", "DT"])
        indptr, indices = get_sparse_edges(
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

    def set_hash_id(self, parameters):
        # TODO
        # print(f"{self.node_count}_{self.edge_count}")
        pass

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










def get_sparse_edges(
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
