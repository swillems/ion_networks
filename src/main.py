#!python

import click
import logging
import sys
import os

import network
import parameter_io


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    context_settings=CONTEXT_SETTINGS,
    help="Analysis of MSMS data with ion-networks."
)
def cli():
    pass


@click.command(
    "create",
    help="Create ion-networks."
)
@click.option(
    "--centroided_file",
    "-c",
    help="The centroided input file (.csv) with centroided ion peaks. "
        "The column INTENSITY and MZ2 need to be present. All other columns "
        "are automically interpreted as precursor selection/separation "
        "dimensions. This flag can be set multiple times to create multiple "
        "ion-networks.",
    required=True,
    multiple=True,
    type=click.File('r')
)
@click.option(
    "--ion_network_file",
    "-i",
    help="A new ion-network file (.hdf). This flag can be set multiple times "
        "to save multiple ion-networks in the same order. If not set, "
        "a new '[raw_file].hdf' file will be created per raw file. WARNING: "
        "This overrides already existing files without confirmation.",
    multiple=True,
    type=click.File('w')
)
# @click.option(
#     "--input_directory",
#     "-d",
#     help="The centroided input file (.csv) with centroided ion peaks. "
#         "The column INTENSITY needs to be present. All other columns are "
#         "automically interpreted as selection/separation dimensions for"
#         "This flag can be set multiple times to create multiple ion-networks.",
#     required=True,
#     is_eager=True,
#     type=click.Path(exists=True)
# )
# https://stackoverflow.com/questions/55584012/python-click-dependent-options-on-another-option
# Use defaults to avoid required?
@click.option(
    "--parameter_file",
    "-p",
    help="A parameter file (.json).",
    type=click.File('r')
)
@click.option(
    "--log_file",
    "-l",
    help="Save the log to a file (.txt).",
    type=click.File('w')
)
def create(
    centroided_file,
    ion_network_file,
    parameter_file,
    log_file
):
    parameters = parameter_io.read(
        file_name=get_file_name(parameter_file),
        default="create"
    )
    with open_logger(get_file_name(log_file)) as logger:
        for index, centroided_csv_file in enumerate(centroided_file):
            centroided_csv_file_name = get_file_name(centroided_csv_file)
            if len(ion_network_file) > index:
                ion_network_file_name = get_file_name(ion_network_file[index])
            else:
                ion_network_file_name = f"{centroided_csv_file_name}.hdf"
            network.Network(
                network_file_name=ion_network_file_name,
                centroided_csv_file_name=centroided_csv_file_name,
                parameters=parameters,
                logger=logger
            )


@click.command(
    "align",
    help="Align ion-networks."
)
@click.option(
    "--ion_network_file",
    "-i",
    help="The ion-network file (.hdf) to align. "
        "This flag can be set multiple times to align multiple "
        "ion-networks pairwise.",
    required=True,
    multiple=True,
    type=click.File('r')
)
@click.option(
    "--alignment_file",
    "-a",
    help="A new alignment file (.hdf) with all pairwise alignments. If not set, "
        "an 'alignment.hdf' file will be created in directory of the first "
        "[ion_network_file]. WARNING: This overrides already existing files "
        "without confirmation.",
    type=click.File('w')
)
@click.option(
    "--parameter_file",
    "-p",
    help="A parameters file (.json).",
    type=click.File('r')
)
@click.option(
    "--log_file",
    "-l",
    help="Save the log to a file (.txt).",
    type=click.File('w')
)
def align(
    ion_network_file,
    alignment_file,
    parameter_file,
    log_file
):
    parameters = parameter_io.read(
        file_name=get_file_name(parameter_file),
        default="create"
    )
    with open_logger(get_file_name(log_file)) as logger:
        alignment_file_name = get_file_name(alignment_file)
        for index, first_ion_network_file in enumerate(ion_network_file[:-1]):
            first_ion_network_file_name = get_file_name(first_ion_network_file)
            first_ion_network = network.Network(
                network_file_name=first_ion_network_file_name,
                parameters=parameters,
                logger=logger
            )
            for second_ion_network_file in ion_network_file[index + 1:]:
                second_ion_network_file_name = get_file_name(second_ion_network_file)
                second_ion_network = network.Network(
                    network_file_name=second_ion_network_file_name,
                    parameters=parameters,
                    logger=logger
                )
                first_ion_network.align(
                    second_ion_network,
                    alignment_file_name=alignment_file_name,
                    parameters=parameters,
                    logger=logger
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
    help="A parameters file (.json).",
    type=click.File('r')
)
@click.option(
    "--log_file",
    "-l",
    help="Save the log to a file (.txt).",
    type=click.File('w')
)
def evidence(
    ion_network_file,
    alignment_file,
    evidence_file,
    parameter_file,
    log_file
):
    parameters = parameter_io.read(
        file_name=get_file_name(parameter_file),
        default="create"
    )
    with open_logger(get_file_name(log_file)) as logger:
        alignment_file_name = get_file_name(alignment_file)
        for index, first_ion_network_file in enumerate(ion_network_file):
            first_ion_network_file_name = get_file_name(first_ion_network_file)
            ion_network = network.Network(
                network_file_name=first_ion_network_file_name,
                parameters=parameters,
                logger=logger
            )
            if len(evidence_file) > index:
                evidence_file_name = get_file_name(evidence_file[index])
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
    help="A parameters file (.json).",
    type=click.File('r')
)
@click.option(
    "--log_file",
    "-l",
    help="Save the log to a file (.txt).",
    type=click.File('w')
)
def show(
    ion_network_file,
    evidence_file,
    parameter_file,
    log_file
):
    # TODO
    parameters = parameter_io.read(
        file_name=get_file_name(parameter_file),
        default="create"
    )
    with open_logger(get_file_name(log_file)) as logger:
        pass


@click.command(
    "gui",
    help="Graphical user interface for ion-networks."
)
def gui():
    # TODO
    pass


def get_file_name(file):
    """
    Return the file name and closes the file.

    Parameters
    ----------
    file : type
        The (un)opened file.

    Returns
    -------
    str
        The file name. If no proper file was provided, None is returned
    """
    try:
        file_name = file.name
        file.close()
    except AttributeError:
        file_name = None
    return file_name


class open_logger(object):

    def __init__(self, log_file_name, log_level=logging.DEBUG, overwrite=False):
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
        logger = logging.getLogger('network_log')
        formatter = logging.Formatter(
            '%(asctime)s > %(message)s'
        )
        logger.setLevel(log_level)
        console_handler = logging.StreamHandler(stream=sys.stdout)
        console_handler.setLevel(log_level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
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
        logger.info("=" * 50)
        logger.info("Executing from the command-line:")
        logger.info(" ".join(sys.argv))
        logger.info("")
        self.logger = logger
        self.log_file_name = log_file_name

    def __enter__(self):
        return self.logger

    def __exit__(self, type, value, traceback):
        if type is not None:
            self.logger.exception("Errors occurred, execution incomplete!")
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
