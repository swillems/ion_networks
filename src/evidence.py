#!python

import logging
import time
import functools

import h5py
import scipy
import numpy as np


class Evidence(object):
    # TODO: Docstring

    def __init__(
        self,
        evidence_file_name,
        ion_network=None,
        alignment=None,
        evidence_ion_networks=[],
        parameters={},
        logger=logging.getLogger('ion_network_log')
    ):
        # TODO: Docstring
        # """
        # Loads an alignment. Alternatively, an alignment is created if list of
        # ion-networks is provided.
        #
        # Parameters
        # ----------
        # alignment_file_name : str
        #     The file name of the alignment.
        # ion_networks : None or iterable[ion-network]
        #     The ion-networks to align. If this is None, the alignment file
        #     is opened as read only.
        # parameters : dict
        #     A dictionary with optional parameters for the alignment of
        #     ion-networks.
        # logger : logging.logger
        #     The logger that indicates all progress
        # """
        self.file_name = evidence_file_name
        self.logger = logger
        if len(evidence_ion_networks) is not None:
            self.evidence_ion_networks(
                ion_network,
                alignment,
                evidence_ion_networks,
                parameters
            )

    def evidence_ion_networks(
        self,
        ion_network,
        alignment,
        evidence_ion_networks,
        parameters
    ):
        # TODO: Docstring
        # """
        # Pairwise align multiple ion-networks against each other.
        #
        # Parameters
        # ----------
        # ion_networks : iterable[ion_network]
        #     The ion-networks that will all be pairwise aligned ageainst each
        #     other.
        # parameters : dict
        #     A dictionary with optional parameters for the alignment of
        #     ion-networks.
        # """
        with h5py.File(
            self.file_name,
            parameters["file_mode"]
        ) as evidence_file:
            ion_network_edges = ion_network.get_edges()
            for second_ion_network in evidence_ion_networks:
                pairwise_alignment = alignment.get_node_alignment(
                    ion_network,
                    second_ion_network,
                    return_as_scipy_csr=True
                )
                second_ion_network_edges = second_ion_network.get_edges()
                indirect = pairwise_alignment.T * second_ion_network_edges * pairwise_alignment
                positive = indirect.multiply(ion_network_edges)
                available = ion_network_edges * pairwise_alignment.T * pairwise_alignment
                negative = available - positive


#
# import numpy as np
# import pandas as pd
# import os
# import network
# import importlib
# import matplotlib
# from matplotlib import pyplot as plt
# from timeit import timeit
# import logging
# import sys
# import seaborn as sns
# import importlib
# import alignment
# import h5py
# formatter = logging.Formatter(
#     '%(asctime)s > %(message)s'
# )
# logger = logging.getLogger('network_log')
# logger.setLevel(logging.DEBUG)
# console_handler = logging.StreamHandler(stream=sys.stdout)
# console_handler.setLevel(logging.DEBUG)
# console_handler.setFormatter(formatter)
# logger.addHandler(console_handler)
#
#
#
#
# importlib.reload(network)
# importlib.reload(alignment)
# inets = []
# in_folder = "/home/sander/Documents/Proteomics/data/ion_networks/ecoli_sonar/ion_networks"
# for file_name in sorted(os.listdir(in_folder)):
#     in_file_name = os.path.join(in_folder, file_name)
#     inet = network.Network(
#         in_file_name
#     )
#     inets.append(inet)
# al = alignment.Alignment(
#     "/home/sander/Documents/Proteomics/data/ion_networks/ecoli_sonar/alignment.hdf"
# #     "/home/sander/Documents/Proteomics/data/ion_networks/dda/dda_sonar_test_align.hdf"
# )
#
#
#
#
# %matplotlib notebook
#
# first_mz2, first_rt1, first_mz1, first_logint = inets[0].get_ion_coordinates(["MZ2", "RT", "MZ1", "LOGINT"])
# second_mz2, second_rt2, second_mz1, second_logint = inets[1].get_ion_coordinates(["MZ2", "RT", "MZ1", "LOGINT"])
# a = al.get_alignment(inets[0], inets[1], return_as_scipy_csr=False)
#
# sns.jointplot(first_mz1, first_mz2, kind="hex", gridsize=500)
# sns.jointplot(first_logint[a[:,0]], second_logint[a[:,1]], kind="hex", gridsize=500)
