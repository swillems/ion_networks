#!python

import click
import logging
import sys
import os
import h5py

import network
import parameter_io


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    context_settings=CONTEXT_SETTINGS,
    help="Analysis of LC-[...]-MSMS data with ion-networks."
)
def cli():
    pass


@click.command(
    "create",
    help="Create ion-networks."
)
@click.option(
    "--input_path",
    "-i",
    help="An input file (.csv) with centroided ion peaks. "
        "The columns RT, MZ2 and LOGINT need to be present. "
        "All other columns are automatically interpreted as precursor "
        "selection/separation dimensions. "
        "This flag can be set multiple times to create multiple ion-networks. "
        "Alternatively, directories can be provided to create "
        "ion-networks for all .csv files contained in these directories.",
    multiple=True,
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    "--output_path",
    "-o",
    help="For each [input.csv] file, an ion-networks is created as "
        "[input.hdf]. If no output_path is set, these are saved in the "
        "original directory as the [input.csv] file. "
        "WARNING: This overrides already existing files without confirmation.",
    type=click.Path(file_okay=False)
)
@click.option(
    "--parameter_file",
    "-p",
    "parameter_file_name",
    help="A file with optional parameters (.json).",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--log_file",
    "-l",
    "log_file_name",
    help="Save the log to a file (.txt).",
    type=click.Path(dir_okay=False)
)
def create(
    input_path,
    output_path,
    parameter_file_name,
    log_file_name
):
    parameters = parameter_io.read(
        file_name=parameter_file_name,
        default="create"
    )
    with open_logger(log_file_name, parameters=parameters) as logger:
        input_csv_files = set()
        for current_path in input_path:
            if os.path.isfile(current_path):
                if current_path.endswith(".csv"):
                    input_csv_files.add(current_path)
            elif os.path.isdir(current_path):
                for current_file_name in os.listdir(current_path):
                    if current_file_name.endswith(".csv"):
                        file_name = os.path.join(
                            current_path,
                            current_file_name
                        )
                        input_csv_files.add(file_name)
        input_csv_files = sorted(input_csv_files)
        file_count = len(input_csv_files)
        logger.info(f"Found {file_count} .csv files to process.")
        for csv_file_name in input_csv_files:
            local_path = os.path.dirname(csv_file_name)
            local_file_name = os.path.basename(csv_file_name)
            if output_path is None:
                ion_network_file_name = os.path.join(
                    local_path,
                    f"{local_file_name[:-4]}.hdf"
                )
            else:
                ion_network_file_name = os.path.join(
                    output_path,
                    f"{local_file_name[:-4]}.hdf"
                )
            network.Network(
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
    help="The ion-network files (.hdf) to align. "
        "This flag can be set multiple times to align multiple "
        "ion-networks pairwise.",
    required=True,
    multiple=True,
    type=click.Path(exists=True)
)
@click.option(
    "--alignment_file",
    "-a",
    "alignment_file_name",
    help="A new alignment file (.hdf) with all pairwise alignments. "
        "WARNING: This overrides already existing files "
        "without confirmation.",
    required=True,
    type=click.Path(dir_okay=False)
)
# @click.option(
#     "--write_mode",
#     "-w",
#     help="Alignment file write mode.",
#     # type=click.File('r')
# )
@click.option(
    "--parameter_file",
    "-p",
    "parameter_file_name",
    help="A file with optional parameters (.json).",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--log_file",
    "-l",
    "log_file_name",
    help="Save the log to a file (.txt).",
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
        input_hdf_files = set()
        for current_path in input_path:
            if os.path.isfile(current_path):
                if current_path.endswith(".hdf"):
                    input_hdf_files.add(current_path)
            elif os.path.isdir(current_path):
                for current_file_name in os.listdir(current_path):
                    if current_file_name.endswith(".hdf"):
                        file_name = os.path.join(
                            current_path,
                            current_file_name
                        )
                        input_hdf_files.add(file_name)
        input_hdf_files = sorted(input_hdf_files)
        file_count = len(input_hdf_files)
        logger.info(f"Found {file_count} .hdf files to process.")
        with h5py.File(
            alignment_file_name,
            parameters["alignment_file_write_mode"]
        ) as alignment_file:
            for index, first_ion_network_file in enumerate(ion_network_file[:-1]):
                first_ion_network_file_name = first_ion_network_file
                first_ion_network = network.Network(
                    network_file_name=first_ion_network_file_name,
                    parameters=parameters,
                    logger=logger
                )
                for second_ion_network_file in ion_network_file[index + 1:]:
                    second_ion_network_file_name = second_ion_network_file
                    second_ion_network = network.Network(
                        network_file_name=second_ion_network_file_name,
                        parameters=parameters,
                        logger=logger
                    )
                    first_ion_network.align(
                        second_ion_network,
                        alignment_file,
                        parameters=parameters,
                    )


@click.command(
    "evidence",
    help="Evidence ion-networks."
)
@click.option(
    "--ion_network_file",
    "-i",
    help="The ion-network file (.hdf) to evidence."
        "This flag can be set multiple times to evidence multiple "
        "ion-networks.",
    required=True,
    multiple=True,
    type=click.File('r')
)
@click.option(
    "--alignment_file",
    "-a",
    help="The alignment file (.hdf) from where to get the evidence. "
        "If a single ion-network was provided, evidence is drawn from all "
        "ion-networks present in the alignment file. If multiple ion-networks "
        "are provided that are present in the alignment file, only those will "
        "be used as evidence for eachother.",
    required=True,
    type=click.File('r')
)
@click.option(
    "--evidence_file",
    "-e",
    help="A new evidence file (.hdf) for the ion-network. If not set, "
        "an '[ion_network_file].evidence.hdf' file will be created per "
        "ion-network. WARNING: This overrides already existing files without "
        "confirmation.",
    multiple=True,
    type=click.File('w')
)
@click.option(
    "--parameter_file",
    "-p",
    "parameter_file_name",
    help="A file with optional parameters (.json).",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--log_file",
    "-l",
    "log_file_name",
    help="Save the log to a file (.txt).",
    type=click.Path(dir_okay=False)
)
def evidence(
    ion_network_file,
    alignment_file,
    evidence_file,
    parameter_file_name,
    log_file_name
):
    parameters = parameter_io.read(
        file_name=parameter_file_name,
        default="evidence"
    )
    with open_logger(log_file_name, parameters=parameters) as logger:
        alignment_file_name = alignment_file
        for index, first_ion_network_file in enumerate(ion_network_file):
            first_ion_network_file_name = first_ion_network_file
            ion_network = network.Network(
                network_file_name=first_ion_network_file_name,
                parameters=parameters,
                logger=logger
            )
            if len(evidence_file) > index:
                evidence_file_name = evidence_file[index]
                ion_network.evidence(
                    alignment_file_name,
                    evidence_file_name=evidence_file_name,
                    parameters=parameters
                )
            else:
                ion_network.evidence(
                    alignment_file_name,
                    parameters=parameters
                )


@click.command(
    "show",
    help="Show ion-networks."
)
@click.option(
    "--ion_network_file",
    "-i",
    help="The ion-network file (.hdf) to show."
        "This flag can be set multiple times to evidence multiple "
        "ion-networks.",
    required=True,
    multiple=True,
    type=click.File('r')
)
@click.option(
    "--evidence_file",
    "-e",
    help="The evidence file (.hdf).",
    required=True,
    type=click.File('r')
)
@click.option(
    "--parameter_file",
    "-p",
    "parameter_file_name",
    help="A file with optional parameters (.json).",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--log_file",
    "-l",
    "log_file_name",
    help="Save the log to a file (.txt).",
    type=click.Path(dir_okay=False)
)
def show(
    ion_network_file,
    evidence_file,
    parameter_file_name,
    log_file_name
):
    # TODO
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
    # TODO
    pass


class open_logger(object):

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
        # TODO
        return self.logger

    def __exit__(self, type, value, traceback):
        # TODO
        if type is not None:
            self.logger.exception("Errors occurred, execution incomplete!")
            handlers = self.logger.handlers[:]
            for handler in handlers:
                handler.close()
                self.logger.removeHandler(handler)
            sys.exit()
        else:
            self.logger.info("Succesfully finished execution.")


if __name__ == "__main__":
    cli.add_command(create)
    cli.add_command(align)
    cli.add_command(evidence)
    cli.add_command(show)
    cli.add_command(gui)
    cli()
