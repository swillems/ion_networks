#!python

import logging
import time
import functools

import h5py
import scipy
import numpy as np


class Alignment(object):
    # TODO: Docstring

    def __init__(
        self,
        alignment_file_name,
        ion_networks=None,
        parameters={},
        logger=logging.getLogger('ion_network_log')
    ):
        """
        Loads an alignment. Alternatively, an alignment is created if list of
        ion-networks is provided.

        Parameters
        ----------
        alignment_file_name : str
            The file name of the alignment.
        ion_networks : None or iterable[ion-network]
            The ion-networks to align. If this is None, the alignment file
            is opened as read only.
        parameters : dict
            A dictionary with optional parameters for the alignment of
            ion-networks.
        logger : logging.logger
            The logger that indicates all progress
        """
        self.file_name = alignment_file_name
        self.logger = logger
        if ion_networks is not None:
            self.align_ion_networks(ion_networks, parameters)

    def align_ion_networks(self, ion_networks, parameters):
        """
        Pairwise align multiple ion-networks against each other.

        Parameters
        ----------
        ion_networks : iterable[ion_network]
            The ion-networks that will all be pairwise aligned ageainst each
            other.
        parameters : dict
            A dictionary with optional parameters for the alignment of
            ion-networks.
        """
        with h5py.File(
            self.file_name,
            parameters["file_mode"]
        ) as alignment_file:
            ion_networks = sorted(ion_networks)
            for index, first_ion_network in enumerate(
                ion_networks[:-1]
            ):
                if first_ion_network.key in alignment_file:
                    first_group = alignment_file[first_ion_network.key]
                else:
                    first_group = alignment_file.create_group(
                        first_ion_network.key
                    )
                for second_ion_network in ion_networks[index + 1:]:
                    if second_ion_network.key in first_group:
                        if parameters["force_overwrite"]:
                            del first_group[second_ion_network.key]
                        else:
                            continue
                    second_group = first_group.create_dataset(
                        second_ion_network.key,
                        data=first_ion_network.align_nodes(
                            second_ion_network,
                            parameters
                        )
                    )
                    second_group.attrs["creation_time"] = time.asctime()
                    for parameter_key, parameter_value in parameters.items():
                        second_group.attrs[parameter_key] = parameter_value

    @functools.lru_cache
    def get_alignment(
        self,
        first_ion_network,
        second_ion_network,
        return_as_scipy_csr=False
    ):
        """
        Get the pairwise alignment between two ion-networks.

        Parameters
        ----------
        first_ion_network : ion_network
            The first ion-networks of the alignment.
        second_ion_network : ion_network
            The second ion-networks of the alignment.
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
        swap = False
        if first_ion_network > second_ion_network:
            first_ion_network, second_ion_network = second_ion_network, first_ion_network
            swap = True
        with h5py.File(self.file_name, "r") as alignment_file:
            alignment = alignment_file[first_ion_network.key][
                second_ion_network.key
            ][...]
            if return_as_scipy_csr:
                alignment = scipy.sparse.csr_matrix(
                    (
                        np.ones(alignment.shape[0], dtype=np.bool),
                        (alignment[:, 0], alignment[:, 1])
                    ),
                    shape=(
                        first_ion_network.node_count,
                        second_ion_network.node_count
                    )
                )
                if swap:
                    alignment = alignment.T.tocsr()
            elif swap:
                alignment = alignment[:, ::-1]
        return alignment


if __name__ == "__main__":
    pass
