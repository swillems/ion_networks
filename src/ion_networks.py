#!python

import logging
import sys
import os

import click

from alignment import Alignment
from network import Network
from evidence import Evidence
import conversion
import parameter_io


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    context_settings=CONTEXT_SETTINGS,
    help="Analysis of LC-[...]-MSMS data with ion-networks."
)
def cli():
    pass


@click.command(
    "convert",
    help="Convert input files."
)
@click.option(
    "--input_path",
    "-i",
    help="An [input.*] file with centroided ion peaks that needs to be "
        "converted to an [input.csv] file before ion-network creation. "
        "This flag can be set multiple times to create multiple files. "
        "Equally, multiple directories can be provided to create multiple "
        "[input.csv] for all [input.*] files contained in these directories.",
    multiple=True,
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    "--output_directory",
    "-o",
    help="For each [input.*] file, an [input.csv] file is created. "
        "WARNING: This overrides already existing files without confirmation.",
    required=True,
    type=click.Path(file_okay=False)
)
@click.option(
    '--data_type',
    '-d',
    help="The data type of the [input.*] file. If this is DDA, an .mgf file "
    "that was centroided with ms-convert is expected with the field "
    "RTINSECONDS as the LC dimension. "
    "For HDMSE, SONAR and SWIM-DIA, a .csv generated with Waters' Apex3d is "
    "expected, typically generated as follows 'Apex3D64.exe -pRawDirName "
    "sample.raw -outputDirName peak_picked_sample_folder -lockMassZ2 785.8426 "
    "-lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 "
    "-leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 "
    "-bEnableCentroids 0'."
    ,
    required=True,
    type=click.Choice(
        ['DDA', 'HDMSE', "SONAR", "SWIMDIA"],
        case_sensitive=False
    )
)
@click.option(
    "--parameter_file",
    "-p",
    "parameter_file_name",
    help="A [parameters.json] file with optional parameters.",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--log_file",
    "-l",
    "log_file_name",
    help="Save the log to a [log.txt] file.",
    type=click.Path(dir_okay=False)
)
def convert(
    input_path,
    output_directory,
    data_type,
    parameter_file_name,
    log_file_name
):
    parameters = parameter_io.read(
        file_name=parameter_file_name,
        default="create"
    )
    with open_logger(log_file_name, parameters=parameters) as logger:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        if data_type == "dda":
            input_files = get_files_with_extension(input_path, extension=".mgf")
        elif data_type == "sonar":
            input_files = get_files_with_extension(input_path, extension=".csv")
        elif data_type == "hdmse":
            input_files = get_files_with_extension(input_path, extension=".csv")
        elif data_type == "swimdia":
            input_files = get_files_with_extension(input_path, extension=".csv")
        file_count = len(input_files)
        logger.info(f"Found {file_count} files to process.")
        for full_file_name in sorted(input_files):
            if data_type == "dda":
                data = conversion.read_mgf(full_file_name, logger)
            elif data_type == "sonar":
                data = conversion.read_sonar(full_file_name, logger)
            elif data_type == "hdmse":
                data = conversion.read_hdmse(full_file_name, logger)
            elif data_type == "swimdia":
                data = conversion.read_swim_dia(full_file_name, logger)
            base_file_name = os.path.basename(full_file_name)
            out_file_name = os.path.join(output_directory, base_file_name)
            conversion.write(data, out_file_name, logger)


@click.command(
    "create",
    help="Create ion-networks."
)
@click.option(
    "--input_path",
    "-i",
    help="An [input.csv] file with centroided ion peaks. "
        "The columns RT, MZ2 and LOGINT need to be present. "
        "All other columns are automatically interpreted as precursor "
        "selection/separation dimensions. "
        "This flag can be set multiple times to create multiple ion-networks. "
        "Equally, multiple directories can be provided to create "
        "ion-networks for all [input.csv] files contained in these directories.",
    multiple=True,
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    "--output_directory",
    "-o",
    help="For each [input.csv] file, an ion-networks is created as "
        "[input.inet.hdf]. If no output_directory is set, these are saved in the "
        "original directory as the [input.csv] files. "
        "WARNING: This overrides already existing files without confirmation.",
    type=click.Path(file_okay=False)
)
@click.option(
    "--parameter_file",
    "-p",
    "parameter_file_name",
    help="A [parameters.json] file with optional parameters.",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--log_file",
    "-l",
    "log_file_name",
    help="Save the log to a [log.txt] file.",
    type=click.Path(dir_okay=False)
)
def create(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name
):
    parameters = parameter_io.read(
        file_name=parameter_file_name,
        default="create"
    )
    with open_logger(log_file_name, parameters=parameters) as logger:
        input_csv_files = get_files_with_extension(input_path, ".csv")
        file_count = len(input_csv_files)
        logger.info(f"Found {file_count} .csv files to process.")
        for csv_file_name in input_csv_files:
            local_path = os.path.dirname(csv_file_name)
            local_file_name = os.path.basename(csv_file_name)
            if output_directory is None:
                ion_network_file_name = os.path.join(
                    local_path,
                    f"{local_file_name[:-4]}.inet.hdf"
                )
            else:
                ion_network_file_name = os.path.join(
                    output_directory,
                    f"{local_file_name[:-4]}.inet.hdf"
                )
            Network(
                network_file_name=ion_network_file_name,
                centroided_csv_file_name=csv_file_name,
                parameters=parameters,
                logger=logger
            )


@click.command(
    "align",
    help="Align ion-networks."
)
@click.option(
    "--input_path",
    "-i",
    help="The [input.inet.hdf] ion-network files to align. "
        "This flag can be set multiple times to align multiple "
        "ion-networks pairwise."
        "Equally, multiple directories can be provided to align all "
        "[input.inet.hdf] files contained in these directories.",
    required=True,
    multiple=True,
    type=click.Path(exists=True)
)
@click.option(
    "--output_file",
    "-o",
    "alignment_file_name",
    help="A new [alignment.hdf] file with all pairwise alignments. "
        "WARNING: This overrides already existing alignments "
        "without confirmation unless 'force_overwrite' "
        "is set to false in a [parameters.json] file.",
    required=True,
    type=click.Path(dir_okay=False)
)
@click.option(
    "--parameter_file",
    "-p",
    "parameter_file_name",
    help="A [parameters.json] file with optional parameters.",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--log_file",
    "-l",
    "log_file_name",
    help="Save the log to a [log.txt] file.",
    type=click.Path(dir_okay=False)
)
def align(
    input_path,
    alignment_file_name,
    parameter_file_name,
    log_file_name
):
    parameters = parameter_io.read(
        file_name=parameter_file_name,
        default="align"
    )
    with open_logger(log_file_name, parameters=parameters) as logger:
        ion_network_names = get_files_with_extension(input_path, ".inet.hdf")
        ion_networks = {Network(file_name) for file_name in ion_network_names}
        file_count = len(ion_networks)
        logger.info(f"Found {file_count} .inet.hdf files to process.")
        Alignment(
            alignment_file_name,
            ion_networks=ion_networks,
            parameters=parameters,
            logger=logger
        )


@click.command(
    "evidence",
    help="Evidence ion-networks."
)
@click.option(
    "--input_path",
    "-i",
    help="The ion-network file (.inet.hdf) to evidence."
        "This flag can be set multiple times to evidence multiple "
        "ion-networks against each other."
        "Alternatively, directories can be provided to evidence all "
        "[.inet.hdf] files contained in these directories.",
    required=True,
    multiple=True,
    type=click.Path(exists=True)
)
@click.option(
    "--alignment_file",
    "-a",
    help="The alignment file (.hdf) from where to get the evidence. "
        "Only those ion-networks with actual alignments will be used as "
        "evidence for eachother.",
    required=True,
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--output_directory",
    "-o",
    help="For each [input.inet.hdf] file, an evidence is created as "
        "[input.evidence.hdf]. If no output_path is set, these are saved in the "
        "original directory as the [input.inet.hdf] files. "
        "WARNING: This overrides already existing files without confirmation.",
    type=click.Path(file_okay=False)
)
@click.option(
    "--parameter_file",
    "-p",
    "parameter_file_name",
    help="A [parameters.json] file with optional parameters.",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--log_file",
    "-l",
    "log_file_name",
    help="Save the log to a [log.txt] file.",
    type=click.Path(dir_okay=False)
)
def evidence(
    input_path,
    alignment_file_name,
    output_directory,
    parameter_file_name,
    log_file_name
):
    parameters = parameter_io.read(
        file_name=parameter_file_name,
        default="evidence"
    )
    with open_logger(log_file_name, parameters=parameters) as logger:
        alignment = Alignment(alignment_file_name)
        ion_network_names = get_files_with_extension(input_path, ".inet.hdf")
        ion_networks = sorted(
            {
                Network(file_name) for file_name in ion_network_names
            }
        )
        file_count = len(ion_networks)
        logger.info(f"Found {file_count} .inet.hdf files to process.")
        for ion_network in ion_networks:
            local_path = os.path.dirname(ion_network.file_name)
            local_file_name = os.path.basename(ion_network.file_name)
            if output_directory is None:
                evidence_file_name = os.path.join(
                    local_path,
                    f"{local_file_name[:-9]}.evidence.hdf"
                )
            else:
                evidence_file_name = os.path.join(
                    output_directory,
                    f"{local_file_name[:-9]}.evidence.hdf"
                )
            Evidence(
                evidence_file_name=evidence_file_name,
                ion_network=ion_network,
                alignment=alignment,
                evidence_ion_networks=ion_networks,
                parameters=parameters,
                logger=logger
            )


@click.command(
    "show",
    help="Show ion-networks."
)
@click.option(
    "--ion_network_file",
    "-i",
    help="The ion-network file (.inet.hdf) to show."
        "This flag can be set multiple times to evidence multiple "
        "ion-networks.",
    required=True,
    multiple=True,
    type=click.File('r')
)
@click.option(
    "--evidence_file",
    "-e",
    help="The evidence file (.evidence.hdf).",
    required=True,
    type=click.File('r')
)
@click.option(
    "--parameter_file",
    "-p",
    "parameter_file_name",
    help="A [parameters.json] file with optional parameters.",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--log_file",
    "-l",
    "log_file_name",
    help="Save the log to a [log.txt] file.",
    type=click.Path(dir_okay=False)
)
def show(
    ion_network_file,
    evidence_file,
    parameter_file_name,
    log_file_name
):
    # TODO: implement
    raise NotImplementedError
    parameters = parameter_io.read(
        file_name=parameter_file_name,
        default="create"
    )
    with open_logger(log_file_name, parameters=parameters) as logger:
        pass


@click.command(
    "gui",
    help="Graphical user interface for ion-networks."
)
def gui():
    # TODO: implement
    raise NotImplementedError


class open_logger(object):
    # TODO: Docstring

    def __init__(
        self,
        log_file_name,
        log_level=logging.DEBUG,
        overwrite=False,
        parameters={}
    ):
        """
        Create a logger to track all progress.

        Parameters
        ----------
        log_file_name : str
            If a log_file_name is provided, the current log is appended to this
            file.
        log_level : int
            The level at which log messages are returned. by default this is
            logging.DEBUG.
        overwrite : bool
            If overwrite is True, the current log is not appended to the file name
            but overwrites it instead.
        """
        logger = logging.getLogger('ion_network_log')
        formatter = logging.Formatter(
            '%(asctime)s > %(message)s'
        )
        logger.setLevel(log_level)
        console_handler = logging.StreamHandler(stream=sys.stdout)
        console_handler.setLevel(log_level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        if log_file_name is None:
            log_file_name = parameters["log_file_name"]
        if log_file_name is not None:
            directory = os.path.dirname(log_file_name)
            if not os.path.exists(directory):
                os.makedirs(directory)
            if overwrite:
                mode = "w"
            else:
                mode = "a"
            file_handler = logging.FileHandler(
                log_file_name,
                mode=mode
            )
            file_handler.setLevel(log_level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            parameters["log_file_name"] = log_file_name
        else:
            parameters["log_file_name"] = ""
        logger.info("=" * 50)
        logger.info("Executing command: ion_networks.py " + " ".join(sys.argv[1:]))
        logger.info("")
        self.logger = logger
        self.log_file_name = log_file_name

    def __enter__(self):
        return self.logger

    def __exit__(self, type, value, traceback):
        if type is not None:
            self.logger.exception("Errors occurred, execution incomplete!")
            handlers = self.logger.handlers[:]
            for handler in handlers:
                handler.close()
                self.logger.removeHandler(handler)
            sys.exit()
        else:
            self.logger.info("Successfully finished execution.")


def get_files_with_extension(input_path, extension=None):
    """
    Get all files with a specific extension from a list of files and folders.

    Parameters
    ----------
    input_path : iterable[str]
        An iterable with files or folders from which all files with a specific
        extension need to be selected.
    extension : str
        The extension of the files of interest.

    Returns
    -------
    list
        A sorted list with unique file names with the specific extension.
    """
    input_files = set()
    for current_path in input_path:
        if os.path.isfile(current_path):
            if (extension is None) or current_path.endswith(extension):
                input_files.add(current_path)
        elif os.path.isdir(current_path):
            for current_file_name in os.listdir(current_path):
                if (extension is None) or current_file_name.endswith(extension):
                    file_name = os.path.join(
                        current_path,
                        current_file_name
                    )
                    input_files.add(file_name)
    return sorted(input_files)


if __name__ == "__main__":
    cli.add_command(convert)
    cli.add_command(create)
    cli.add_command(align)
    cli.add_command(evidence)
    cli.add_command(show)
    cli.add_command(gui)
    cli()
