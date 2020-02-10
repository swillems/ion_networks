#!python

import os

import network
import evidence
import gui as inet_gui
import utils


def convert_data_to_csv(
    input_path,
    output_directory,
    data_type,
    parameter_file_name,
    log_file_name
):
    """
    Convert centroided MSMS data to a unified csv that can be read as an
    ion-network.

    Parameters
    ----------
    input_path : iterable[str]
        An iterable with file and/or folder names.
    output_directory : str or None
        If provided, all new files will be saved in this directory.
    data_type : str
        The data type of the input files. Options are:
            'dda'
            'sonar'
            'hdmse'
            'swimdia'
    parameter_file_name : str or None
        If provided, parameters will be read from this file.
    log_file_name : str or None
        If provided, all logs will be written to this file.
    """
    parameters = utils.read_parameters_from_json_file(file_name=parameter_file_name)
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        input_file_names = utils.get_file_names_with_extension(
            input_path,
            extension=utils.DATA_TYPE_FILE_EXTENSIONS[data_type]
        )
        file_count = len(input_file_names)
        logger.info(f"Found {file_count} files to process.")
        for input_file_name in sorted(input_file_names):
            data = utils.read_data_from_file(data_type, input_file_name, logger)
            file_name_base = os.path.splitext(
                os.path.basename(input_file_name)
            )[0]
            output_file_name = os.path.join(
                output_directory,
                f"{file_name_base}.inet.csv"
            )
            utils.write_data_to_csv_file(data, output_file_name, logger)


def create_ion_network(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name
):
    """
    Create ion-networks from unified csv files.

    Parameters
    ----------
    input_path : iterable[str]
        An iterable with file and/or folder names.
    output_directory : str or None
        If provided, all new files will be saved in this directory.
    parameter_file_name : str or None
        If provided, parameters will be read from this file.
    log_file_name : str or None
        If provided, all logs will be written to this file.
    """
    parameters = utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="create"
    )
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        input_file_names = utils.get_file_names_with_extension(
            input_path,
            ".inet.csv"
        )
        file_count = len(input_file_names)
        logger.info(f"Found {file_count} .inet.csv files to process.")
        for csv_file_name in input_file_names:
            local_file_name = os.path.basename(csv_file_name)
            if output_directory is None:
                output_path = os.path.dirname(csv_file_name)
            else:
                output_path = output_directory
            ion_network_file_name = os.path.join(
                output_path,
                f"{local_file_name[:-9]}.inet.hdf"
            )
            network.Network(
                network_file_name=ion_network_file_name,
                centroided_csv_file_name=csv_file_name,
                parameters=parameters,
                logger=logger
            )


def evidence_ion_network(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name
):
    """
    Evidence ion-networks with each other.

    Parameters
    ----------
    input_path : iterable[str]
        An iterable with file and/or folder names.
    output_directory : str or None
        If provided, all new files will be saved in this directory.
    parameter_file_name : str or None
        If provided, parameters will be read from this file.
    log_file_name : str or None
        If provided, all logs will be written to this file.
    """
    parameters = utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="evidence"
    )
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        input_file_names = utils.get_file_names_with_extension(
            input_path,
            ".inet.hdf"
        )
        file_count = len(input_file_names)
        logger.info(f"Found {file_count} .inet.hdf files to process.")
        ion_networks = [
            network.Network(file_name) for file_name in input_file_names
        ]
        for ion_network in ion_networks:
            local_file_name = os.path.basename(ion_network.file_name)
            if output_directory is None:
                output_path = os.path.dirname(ion_network.file_name)
            else:
                output_path = output_directory
            evidence_file_name = os.path.join(
                output_path,
                f"{local_file_name[:-9]}.evidence.hdf"
            )
            evidence.Evidence(
                evidence_file_name=evidence_file_name,
                ion_network=ion_network,
                # alignment=align_ion_networks(ion_networks, parameters),
                evidence_ion_networks=ion_networks,
                parameters=parameters,
                logger=logger
            )


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


def show_ion_network(
    ion_network_file_name,
    evidence_file_name,
    parameter_file_name,
    log_file_name
):
    # TODO: Docstring
    parameters = utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="show"
    )
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        inet_gui.GUI(
            network.Network(ion_network_file_name),
            evidence.Evidence(evidence_file_name),
            logger
        )


def start_gui():
    # TODO: implement
    raise NotImplementedError
