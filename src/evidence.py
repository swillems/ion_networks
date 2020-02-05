#!python

import logging
import os
import time

import h5py
import numpy as np

from network import Network


class Evidence(object):
    """
    An evidence set containing positive and negative edge evidence, as well as
    node evidence.
    """

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
        #     The logger that indicates all progress.
        # """
        self.file_name = evidence_file_name
        self.logger = logger
        if len(evidence_ion_networks) > 0:
            self.set_ion_networks_evidence(
                ion_network,
                alignment,
                evidence_ion_networks,
                parameters
            )

    def set_ion_networks_evidence(
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
        directory = os.path.dirname(self.file_name)
        if not os.path.exists(directory):
            os.makedirs(directory)
        with h5py.File(
            self.file_name,
            parameters["file_mode"]
        ) as evidence_file:
            ion_network_edges = ion_network.get_edges()
            ion_network_edges_as_int = ion_network_edges.astype(np.int)
            left_node_indices, right_node_indices = ion_network_edges.nonzero()
            if "ion_network_key" not in evidence_file.attrs:
                evidence_file.attrs["ion_network_key"] = ion_network.key
            for second_ion_network in evidence_ion_networks:
                if ion_network == second_ion_network:
                    continue
                if second_ion_network.key in evidence_file:
                    if parameters["force_overwrite"]:
                        del evidence_file[second_ion_network.key]
                    else:
                        continue
                self.logger.info(
                    f"Evidencing {ion_network.file_name} with "
                    f"{second_ion_network.file_name}"
                )
                evidence_group = evidence_file.create_group(
                    second_ion_network.key
                )
                pairwise_alignment = alignment.get_alignment(
                    ion_network,
                    second_ion_network,
                    return_as_scipy_csr=True
                )
                positive = pairwise_alignment * second_ion_network.get_edges() * pairwise_alignment.T
                positive = (positive + positive.T).multiply(ion_network_edges)
                positive_mask = (ion_network_edges_as_int + positive).data == 2
                alignment_mask = np.diff(pairwise_alignment.indptr) > 0
                negative_mask = (
                    alignment_mask[left_node_indices] & alignment_mask[right_node_indices]
                ) & (~positive_mask)
                evidence_group.create_dataset(
                    "positive_edge_mask",
                    data=positive_mask,
                    compression="lzf",
                )
                evidence_group.create_dataset(
                    "negative_edge_mask",
                    data=negative_mask,
                    compression="lzf",
                )
                evidence_group.create_dataset(
                    "node_mask",
                    data=alignment_mask,
                    compression="lzf",
                )
                evidence_group.attrs["creation_time"] = time.asctime()

    def get_evidence(
        self,
        network_keys=None,
        kind=["positive_edges", "negative_edges", "nodes"],
        return_total=False
    ):
        # TODO: Docstring
        if network_keys is None:
            network_keys = self.network_keys
        if isinstance(network_keys, Network):
            # FIXME: Does not work?
            network_keys = [network_keys.key]
        elif isinstance(network_keys, str):
            network_keys = [network_keys]
        else:
            network_keys = [
                network_key if isinstance(
                    network_key, str
                ) else network_key.key for network_key in network_keys
            ]
        arrays = {
            "positive_edges": [],
            "negative_edges": [],
            "nodes": [],
        }
        single_kind = False
        if isinstance(kind, str):
            kind = [kind]
            single_kind = True
        with h5py.File(self.file_name, "r") as evidence_file:
            for network_key in network_keys:
                if "positive_edges" in kind:
                    arrays["positive_edges"].append(
                        evidence_file[network_key]["positive_edge_mask"][...]
                    )
                if "negative_edges" in kind:
                    arrays["negative_edges"].append(
                        evidence_file[network_key]["negative_edge_mask"][...]
                    )
                if "nodes" in kind:
                    arrays["nodes"].append(
                        evidence_file[network_key]["node_mask"][...]
                    )
        if return_total:
            arrays = {
                key: np.sum(value, axis=0) for key, value in arrays.items()
            }
        if (len(network_keys) == 1) or return_total:
            if return_total:
                if single_kind:
                    return arrays[kind[0]]
                return tuple(arrays[k] for k in kind)
            if single_kind:
                return arrays[kind[0]][0]
            return tuple(arrays[k][0] for k in kind)
        else:
            if single_kind:
                return arrays[kind[0]]
            return tuple(arrays[k] for k in kind)

    @property
    def network_keys(self):
        """
        Get a sorted list with all the ion-network keys providing evidence.

        Returns
        -------
        list[str]
            A sorted list with the keys of all ion-network.
        """
        with h5py.File(self.file_name, "r") as network_file:
            ion_networks = list(network_file)
            # FIXME: Can be deleted
            ion_networks.remove("parameters")
        return ion_networks

    @property
    def network_count(self):
        """
        Get the number of ion-network providing evidence.

        Returns
        -------
        int
            The number of ion-network providing evidence.
        """
        return len(self.network_keys)
