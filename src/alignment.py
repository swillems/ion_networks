#!python


import h5py
import logging
import time


class Alignment(object):

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
        ion_networks : None or list[ion-network]
            The ion-networks to align. If this is None, the alignment file
            is opened as read only.
        parameters : dict
            A dictionary with optional parameters for the creation of an
            ion-network.
        logger : logging.logger
            The logger that indicates all progress
        """
        self.file_name = alignment_file_name
        self.logger = logger
        if ion_networks is not None:
            self.align_ion_networks(ion_networks, parameters)

    def align_ion_networks(self, ion_networks, parameters):
        # TODO
        with h5py.File(
            self.file_name,
            parameters["alignment_file_write_mode"]
        ) as alignment_file:
            ion_networks = sorted(ion_networks)
            for index, first_ion_network in enumerate(
                ion_networks[:-1]
            ):
                for second_ion_network in ion_networks[index + 1:]:
                    if first_ion_network.key not in list(alignment_file.keys()):
                        alignment_group = alignment_file.create_group(
                            first_ion_network.key
                        )
                        alignment_group = alignment_group.create_group(
                            second_ion_network.key
                        )
                    elif second_ion_network.key not in list(
                        alignment_file[first_ion_network.key].keys()
                    ):
                        alignment_group = alignment_file[
                            first_ion_network.key
                        ].create_group(
                            second_ion_network.key
                        )
                    else:
                        alignment_group = alignment_file[
                            first_ion_network.key
                        ][
                            second_ion_network.key
                        ]
                    if parameters["alignment_file_force_overwrite_nodes"]:
                        align_nodes = True
                    else:
                        try:
                            alignment_group["nodes"]["first_indices"]
                            alignment_group["nodes"]["second_indices"]
                            align_nodes = False
                        except KeyError:
                            align_nodes = True
                    if align_nodes:
                        if "nodes" in list(alignment_group.keys()):
                            del alignment_group["nodes"]
                        alignment_group.create_group("nodes")
                        first_indices, second_indices = first_ion_network.align_nodes(
                            second_ion_network,
                            parameters
                        )
                        alignment_group["nodes"].create_dataset(
                            "first_indices",
                            data=first_indices
                        )
                        alignment_group["nodes"].create_dataset(
                            "second_indices",
                            data=second_indices
                        )
                        alignment_group["nodes"].attrs[
                            "creation_time"
                        ] = time.asctime()
                    if parameters["alignment_file_force_overwrite_edges"]:
                        align_edges = True
                    else:
                        try:
                            alignment_group["edges"]["first_indices"]
                            alignment_group["edges"]["second_indices"]
                            align_edges = False
                        except KeyError:
                            align_edges = True
                    if align_edges:
                        if "edges" in list(alignment_group.keys()):
                            del alignment_group["edges"]
                        alignment_group.create_group("edges")
                        first_indices, second_indices = first_ion_network.align_edges(
                            second_ion_network,
                            parameters
                        )
                        alignment_group["edges"].create_dataset(
                            "first_indices",
                            data=first_indices
                        )
                        alignment_group["edges"].create_dataset(
                            "second_indices",
                            data=second_indices
                        )
                    if align_nodes or align_edges:
                        if "parameters" in list(alignment_group.keys()):
                            del alignment_group["parameters"]
                        alignment_group["parameters"].attrs[
                            "creation_time"
                        ] = time.asctime()
                        for parameter_key, parameter_value in parameters.items():
                            alignment_group["parameters"].attrs[
                                parameter_key
                            ] = parameter_value

    def get_alignment(self, first_ion_network, second_ion_network):
        # TODO
        if first_ion_network > second_ion_network:
            first_ion_network, second_ion_network = second_ion_network, first_ion_network
            swap = True
        else:
            swap = False
        with h5py.File(self.file_name, "r") as alignment_file:
            alignment_group = alignment_file[
                first_ion_network.key
            ][
                second_ion_network.key
            ]
            nodes_first_indices = alignment_group["nodes/first_indices"][...]
            nodes_second_indices = alignment_group["nodes/second_indices"][...]
            edges_first_indices = alignment_group["edges/first_indices"][...]
            edges_second_indices = alignment_group["edges/second_indices"][...]
        if swap:
            return nodes_second_indices, nodes_first_indices, edges_second_indices, edges_first_indices
        else:
            return nodes_first_indices, nodes_second_indices, edges_first_indices, edges_second_indices


if __name__ == "__main__":
    pass
