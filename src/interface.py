#!python

# builtin
import os
# local
import network
import evidence
import gui
import utils
import user_interface


def convert_data_formats_to_csvs(
    input_path,
    data_type,
    output_directory,
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
    if parameter_file_name is None:
        parameter_file_name = ""
    if output_directory is None:
        output_directory = ""
    if log_file_name is None:
        log_file_name = ""
    parameters = utils.read_parameters_from_json_file(
        file_name=parameter_file_name
    )
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        # logger.info("Running command: ")
        # logger.info("")
        if output_directory != "":
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
            if output_directory == "":
                output_path = os.path.dirname(input_file_name)
            else:
                output_path = output_directory
            output_file_name = os.path.join(
                output_path,
                f"{file_name_base}.inet.csv"
            )
            utils.write_data_to_csv_file(data, output_file_name, logger)


def create_ion_networks(
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
    if parameter_file_name is None:
        parameter_file_name = ""
    if output_directory is None:
        output_directory = ""
    if log_file_name is None:
        log_file_name = ""
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
            if output_directory == "":
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


def evidence_ion_networks(
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
    if parameter_file_name is None:
        parameter_file_name = ""
    if output_directory is None:
        output_directory = ""
    if log_file_name is None:
        log_file_name = ""
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
        evidence_files = []
        for ion_network in ion_networks:
            local_file_name = os.path.basename(ion_network.file_name)
            if output_directory == "":
                output_path = os.path.dirname(ion_network.file_name)
            else:
                output_path = output_directory
            evidence_file_name = os.path.join(
                output_path,
                f"{local_file_name[:-9]}.evidence.hdf"
            )
            evidence_files.append(
                evidence.Evidence(
                    evidence_file_name=evidence_file_name,
                    ion_network=ion_network,
                    parameters=parameters,
                    logger=logger
                )
            )
        for index, evidence_file in enumerate(evidence_files[:-1]):
            for secondary_evidence_file in evidence_files[index + 1:]:
                evidence_file.mutual_collect_evidence_from(
                    secondary_evidence_file,
                    parameters=parameters,
                    logger=logger
                )


def show_ion_network(
    ion_network_file_name,
    evidence_file_name,
    parameter_file_name,
    log_file_name
):
    # TODO: Docstring
    # TODO: Implement
    if parameter_file_name is None:
        parameter_file_name = ""
    if log_file_name is None:
        log_file_name = ""
    parameters = utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
    )
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        inet = network.Network(ion_network_file_name)
        evi = evidence.Evidence(
            evidence_file_name=evidence_file_name,
            ion_network=inet,
            parameters=parameters,
            logger=logger
        )
        gui.GUI(
            inet,
            evi,
            logger
        )


def start_gui():
    # TODO: Docstring
    with user_interface.GUI() as gui:
        gui.run()
