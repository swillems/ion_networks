#!python

# builtin
import logging
import os
import time
# external
import h5py
import numpy as np
import scipy
# local
from network import Network


class Evidence(object):
    """
    An evidence set containing positive and negative edge evidence, as well as
    node evidence.
    """

    def __init__(
        self,
        evidence_file_name,
        ion_network,
        parameters={},
        logger=logging.getLogger('ion_network_log')
    ):
        # TODO: Docstring
        self.file_name = evidence_file_name
        self.logger = logger
        self.ion_network = ion_network

    def mutual_collect_evidence_from(
        self,
        other,
        parameters={},
        logger=logging.getLogger('ion_network_log')
    ):
        # TODO: Docstring
        directory = os.path.dirname(self.file_name)
        if not os.path.exists(directory):
            os.makedirs(directory)
        if self.is_evidenced_with(other) or other.is_evidenced_with(self):
            # TODO: Should always be mutual evidence
            if parameters["force_overwrite"]:
                self.remove_evidence_from(other)
                other.remove_evidence_from(self)
            else:
                logger.info(
                    f"Evidence for {self.file_name} and {other.file_name} "
                    f"has already been collected"
                )
                return
        pairwise_alignment = self.ion_network.align_nodes(
            other.ion_network,
            parameters
        )
        pairwise_alignment_T = pairwise_alignment.T.tocsr()
        self_edges = self.ion_network.get_edges()
        other_edges = other.ion_network.get_edges()
        self.align_edges(
            other,
            pairwise_alignment,
            pairwise_alignment_T,
            self_edges,
            other_edges,
        )
        other.align_edges(
            self,
            pairwise_alignment_T,
            pairwise_alignment,
            other_edges,
            self_edges,
        )

    def is_evidenced_with(self, other):
        # TODO: Docstring
        result = False
        with h5py.File(self.file_name, "a") as evidence_file:
            if other.ion_network.file_name_base in evidence_file:
                result = True
        return result

    def remove_evidence_from(self, other):
        # TODO: Docstring
        if self.is_evidenced_with(other):
            with h5py.File(self.file_name, "a") as evidence_file:
                del evidence_file[other.ion_network.file_name_base]

    def align_edges(
        self,
        other,
        pairwise_alignment,
        pairwise_alignment_T,
        self_edges,
        other_edges,
    ):
        # TODO: Docstring
        self.logger.info(
            f"Collecting evidence for {self.ion_network.file_name} from "
            f"{other.ion_network.file_name}."
        )
        with h5py.File(self.file_name, "a") as evidence_file:
            evidence_group = evidence_file.create_group(
                other.ion_network.file_name_base
            )
            self_edges_as_int = self_edges.astype(np.int)
            left_node_indices, right_node_indices = self_edges.nonzero()
            positive = pairwise_alignment * other_edges * pairwise_alignment_T
            positive = (positive + positive.T).multiply(self_edges)
            positive_mask = (self_edges_as_int + positive).data == 2
            alignment_mask = np.diff(pairwise_alignment.indptr) > 0
            negative_mask = (
                alignment_mask[left_node_indices] & alignment_mask[right_node_indices]
            ) & (~positive_mask)
            evidence_group.create_dataset(
                "positive_edges",
                data=positive_mask,
                compression="lzf",
            )
            evidence_group.create_dataset(
                "negative_edges",
                data=negative_mask,
                compression="lzf",
            )
            evidence_group.create_dataset(
                "aligned_nodes",
                data=pairwise_alignment_T.indices,
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
            # TODO: earlier summation saves memory!
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

    # def get_alignment(
    #     self,
    #     first_ion_network,
    #     second_ion_network,
    #     return_as_scipy_csr=False
    # ):
    #     """
    #     Get the pairwise alignment between two ion-networks.
    #
    #     Parameters
    #     ----------
    #     first_ion_network : ion_network
    #         The first ion-networks of the alignment.
    #     second_ion_network : ion_network
    #         The second ion-networks of the alignment.
    #     return_as_scipy_csr : bool
    #         If True, return the alignment as a scipy.sparse.csr_matrix,
    #         otherwise, return the alignment as a 2-dimensional np.ndarray.
    #
    #     Returns
    #     ----------
    #     np.ndarray or scipy.sparse.csr_matrix(bool)
    #         A 2-dimensional array of shape (2, n) with n the number of nodes
    #         aligned. The first column are the indices of the first ion-network
    #         and the second column contains the indices of the second
    #         ion-network. If return_as_scipy_csr is True, a sparse csr matrix
    #         is created from this array before it is returned.
    #     """
    #     swap = False
    #     if first_ion_network > second_ion_network:
    #         first_ion_network, second_ion_network = second_ion_network, first_ion_network
    #         swap = True
    #     with h5py.File(self.file_name, "r") as alignment_file:
    #         alignment = alignment_file[first_ion_network.key][
    #             second_ion_network.key
    #         ][...]
    #     if return_as_scipy_csr:
    #         alignment = scipy.sparse.csr_matrix(
    #             (
    #                 np.ones(alignment.shape[0], dtype=np.bool),
    #                 (alignment[:, 0], alignment[:, 1])
    #             ),
    #             shape=(
    #                 first_ion_network.node_count,
    #                 second_ion_network.node_count
    #             )
    #         )
    #         if swap:
    #             alignment = alignment.T.tocsr()
    #     elif swap:
    #         alignment = alignment[:, ::-1]
    #     return alignment


# def align_ion_networks(self, ion_networks, parameters):
#     """
#     Pairwise align multiple ion-networks against each other.
#
#     Parameters
#     ----------
#     ion_networks : iterable[ion_network]
#         The ion-networks that will all be pairwise aligned ageainst each
#         other.
#     parameters : dict
#         A dictionary with optional parameters for the alignment of
#         ion-networks.
#     """
#     with h5py.File(
#         self.file_name,
#         parameters["file_mode"]
#     ) as alignment_file:
#         ion_networks = sorted(ion_networks)
#         for index, first_ion_network in enumerate(
#             ion_networks[:-1]
#         ):
#             if first_ion_network.key in alignment_file:
#                 first_group = alignment_file[first_ion_network.key]
#             else:
#                 first_group = alignment_file.create_group(
#                     first_ion_network.key
#                 )
#             for second_ion_network in ion_networks[index + 1:]:
#                 if second_ion_network.key in first_group:
#                     if parameters["force_overwrite"]:
#                         del first_group[second_ion_network.key]
#                     else:
#                         continue
#                 second_group = first_group.create_dataset(
#                     second_ion_network.key,
#                     data=first_ion_network.align_nodes(
#                         second_ion_network,
#                         parameters
#                     ),
#                     compression="lzf"
#                 )
#                 second_group.attrs["creation_time"] = time.asctime()
#                 dimension_overlap = first_ion_network.dimension_overlap(
#                     second_ion_network
#                 )
#                 for parameter_key, parameter_value in parameters.items():
#                     if parameter_key.startswith("max_edge_deviation"):
#                         if parameter_key[24:] not in dimension_overlap:
#                             continue
#                     second_group.attrs[parameter_key] = parameter_value
