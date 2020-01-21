#!python

import os
import sys
import logging

import h5py
import scipy
import pandas as pd
import numpy as np
import multiprocessing as mp


class Network(object):

    def __init__(
        self,
        network_file_name,
        centroided_csv_file_name=None,
        parameters=None
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
            WARNING: If provided, any existing ion-networks is overwritten with
            the data from this centroided .csv file.
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        """
        self.parameters = parameters
        self.file_name = network_file_name
        if centroided_csv_file_name is not None:
            data = self.read_centroided_csv_file(centroided_csv_file_name)
            self.create_from_data(data)

    def read_centroided_csv_file(self, centroided_csv_file_name):
        # TODO
        data = pd.read_csv(
            centroided_csv_file_name,
            engine="c",
            dtype=np.float,
        )
        return data

    def create_from_data(self, data):
        # TODO
        directory = os.path.dirname(self.file_name)
        if not os.path.exists(directory):
            os.makedirs(directory)
        with h5py.File(self.file_name, "w") as network_file:
            node_group = network_file.create_group("nodes")
            edge_group = network_file.create_group("edges")
            self.write_nodes(node_group, data)
            self.write_edges(edge_group)

    def write_nodes(self, node_group, data):
        # TODO
        for column in data.columns:
            values = data[column]
            attribute_group = node_group.create_group(column)
            attribute_group.create_dataset(
                "raw",
                data=values
            )

    def write_edges(self, edge_group):
        # TODO
        pass

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
