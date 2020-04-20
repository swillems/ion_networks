#!python

# builtin
import time
# external
import numpy as np
import scipy
# local
import network
import utils


class Evidence(utils.HDF_MS_Run_File):
    """
    An evidence set containing positive and negative edge evidence, as well as
    node evidence.
    """

    def __init__(
        self,
        reference,
        new_file=False,
        is_read_only=True,
        logger=None
    ):
        # TODO: Docstring
        file_name = self.convert_reference_to_trimmed_file_name(reference)
        super().__init__(
            f"{file_name}.evidence.hdf",
            new_file,
            is_read_only,
            logger,
        )
        # self.ion_network = network.Network(reference)
        self.ion_network = network.Network(file_name + ".inet.hdf")
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
                pass
                # self.remove_evidence_from(other)
                # other.remove_evidence_from(self)
            else:
                self.logger.info(
                    f"Evidence for {self.file_name} and {other.file_name} "
                    f"has already been collected."
                )
                return False
        return True

    def is_evidenced_with(self, other):
        # TODO: Docstring
        return other.run_name in self.get_group_list

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
        if group_name == "":
            ions = np.zeros(self.ion_network.node_count, dtype=np.int)
            for group_name in self.get_group_list:
                ion_mask = self.get_dataset(
                    "aligned_nodes",
                    parent_group_name=group_name,
                )
                ions[ion_mask] += 1
        else:
            ions = self.get_dataset(
                "aligned_nodes",
                parent_group_name=group_name,
            )
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
        if group_name == "":
            edges = np.zeros(self.ion_network.edge_count, dtype=np.int)
            for group_name in self.get_group_list:
                edges += self.get_dataset(
                    mask_name,
                    parent_group_name=group_name,
                )
        else:
            edges = self.get_dataset(
                mask_name,
                parent_group_name=group_name,
            )
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
        return self.get_group_list

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
