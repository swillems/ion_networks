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
import network


class Evidence(object):
    """
    An evidence set containing positive and negative edge evidence, as well as
    node evidence.
    """

    def __init__(
        self,
        evidence_file_name=None,
        ion_network=None,
        parameters={},
        logger=logging.getLogger()
    ):
        # TODO: Docstring
        self.logger = logger
        if ion_network is not None:
            self.ion_network = ion_network
            self.file_name = os.path.join(
                ion_network.file_name_path,
                ion_network.file_name_base
            ) + ".evidence.hdf"
        elif evidence_file_name is not None:
            self.file_name = evidence_file_name
            self.ion_network = network.Network(
                os.path.join(
                    self.file_name_path,
                    self.file_name_base
                ) + ".inet.hdf"
            )
        else:
            raise ValueError("Evidence file needs a corresponding ion-network.")
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
                self.remove_evidence_from(other)
                other.remove_evidence_from(self)
            else:
                self.logger.info(
                    f"Evidence for {self.file_name} and {other.file_name} "
                    f"has already been collected."
                )
                return False
        return True

    def is_evidenced_with(self, other):
        # TODO: Docstring
        with h5py.File(self.file_name, "a") as evidence_file:
            if other.ion_network.file_name_base in evidence_file:
                return True
        return False

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
        aligned_ion_indices,
        self_edges,
        other_edges,
        parameters
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
            positive = pairwise_alignment * other_edges * pairwise_alignment_T
            positive = (positive + positive.T).multiply(self_edges)
            positive_mask = (self_edges.astype(np.int8) + positive).data == 2
            alignment_mask = pairwise_alignment.indptr[:-1] != pairwise_alignment.indptr[1:]
            left_node_indices, right_node_indices = self_edges.nonzero()
            negative_mask = alignment_mask[
                left_node_indices
            ] & alignment_mask[
                right_node_indices
            ] & ~positive_mask
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
                data=aligned_ion_indices,
                compression="lzf",
            )
            dimension_overlap = self.ion_network.dimension_overlap(
                other.ion_network
            )
            evidence_group.attrs["creation_time"] = time.asctime()
            for parameter_key, parameter_value in parameters.items():
                if parameter_key.startswith("max_alignment_absolute_error_"):
                    if parameter_key[24:] not in dimension_overlap:
                        continue
                evidence_group.attrs[parameter_key] = parameter_value

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
            other.ion_network.file_name_base,
            return_as_mask=False
        )
        other_ions = other.get_aligned_nodes_from_group(
            self.ion_network.file_name_base,
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
        with h5py.File(self.file_name, "r") as evidence_file:
            if group_name == "":
                ions = np.zeros(self.ion_network.node_count, dtype=np.int)
                for group_name in evidence_file:
                    ions[evidence_file[group_name]["aligned_nodes"][...]] += 1
            else:
                ions = evidence_file[group_name]["aligned_nodes"][...]
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
        with h5py.File(self.file_name, "r") as evidence_file:
            if group_name == "":
                edges = np.zeros(self.ion_network.edge_count, dtype=np.int)
                for group_name in evidence_file:
                    edges += evidence_file[group_name][mask_name][...]
            else:
                edges = evidence_file[group_name][mask_name][...]
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
        with h5py.File(self.file_name, "r") as evidence_file:
            ion_networks = list(evidence_file)
        return ion_networks

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

    @property
    def file_name_base(self):
        """
        Get the file name base of the evidence.

        Returns
        -------
        str
            The file name base of the evidence as a string.
        """
        return os.path.basename(self.file_name)[:-13]

    @property
    def file_name_path(self):
        """
        Get the file name path of the evidence.

        Returns
        -------
        str
            The file name path of the evidence as a string.
        """
        return os.path.dirname(self.file_name)
