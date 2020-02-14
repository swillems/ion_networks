#!python

# external
import click
# local
import interface


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    context_settings=CONTEXT_SETTINGS,
    help="Analysis of LC-[...]-MSMS data with ion-networks."
)
def cli():
    pass


@click.command(
    "convert",
    help="Convert [input.*] files with centroided ions to unified "
        "[input.inet.csv] csv files.",
    short_help="Convert various input formats to unified input."
)
@click.option(
    "--input_path",
    "-i",
    help="An [input.*] file with centroided ion peaks that needs to be "
        "converted to a unified [input.inet.csv] file. "
        "Individual files can be provided, as well as folders."
        "This flag can be set multiple times.",
    multiple=True,
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    "--output_directory",
    "-o",
    help="For each [input.*] file, an [input.inet.csv] file is created. "
        "If no output directory is provided, each [input.inet.csv] file is "
        "placed in the same folder as its corresponding [input.*] file. "
        "WARNING: This overrides already existing files without confirmation.",
    type=click.Path(file_okay=False)
)
@click.option(
    '--data_type',
    '-d',
    help="The data type of the [input.*] file. If this is DDA, an .mgf file "
    "that was centroided with ms-convert is expected with the field "
    "RTINSECONDS as LC coordinate. "
    "For HDMSE, SONAR and SWIM-DIA, a .csv generated with Waters' Apex3d is "
    "expected, typically generated as follows 'Apex3D64.exe -pRawDirName "
    "sample.raw -outputDirName peak_picked_sample_folder -lockMassZ2 785.8426 "
    "-lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 "
    "-leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 "
    "-bEnableCentroids 0'.",
    required=True,
    type=click.Choice(
        ['DDA', 'HDMSE', "SONAR", "SWIMDIA"],
        case_sensitive=True
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
    help="Save the log to a [log.txt] file. "
        "This log file can also be supplied through a [parameters.json] file. "
        "If the log file already exists, the new log data is appended.",
    type=click.Path(dir_okay=False)
)
def convert(
    input_path,
    output_directory,
    data_type,
    parameter_file_name,
    log_file_name
):
    interface.convert_data_formats_to_csvs(
        input_path,
        output_directory,
        data_type,
        parameter_file_name,
        log_file_name
    )


@click.command(
    "create",
    help="Create [input.inet.hdf] ion-network files from unified "
        "[input.inet.csv] files.",
    short_help="Create ion-networks from unified input."
)
@click.option(
    "--input_path",
    "-i",
    help="A unified [input.inet.csv] file with centroided ion peaks. "
        "The columns PRECURSOR_RT, FRAGMENT_MZ and FRAGMENT_LOGINT always "
        "need to be present. "
        "All other columns (e.g. PRECURSOR_MZ or PRECURSOR_DT) are "
        "automatically interpreted as precursor dimensions. "
        "Individual files can be provided, as well as folders."
        "This flag can be set multiple times.",
    multiple=True,
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    "--output_directory",
    "-o",
    help="For each [input.inet.csv] file, an [input.inet.hdf] ion-network "
        "file is created. "
        "If no output directory is provided, each [input.inet.hdf] file is "
        "placed in the same folder as its corresponding [input.inet.csv] file. "
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
    help="Save the log to a [log.txt] file. "
        "This log file can also be supplied through a [parameters.json] file. "
        "If the log file already exists, the new log data is appended.",
    type=click.Path(dir_okay=False)
)
def create(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name
):
    interface.create_ion_networks(
        input_path,
        output_directory,
        parameter_file_name,
        log_file_name
    )


@click.command(
    "evidence",
    help="Collect pairwise evidence for [input.inet.hdf] ion-network files as "
        "[input.evidence.hdf] evidence files.",
    short_help="Collect evidence for ion-networks."
)
@click.option(
    "--input_path",
    "-i",
    help="An [input.inet.hdf] ion-network file."
        "Individual files can be provided, as well as folders."
        "This flag can be set multiple times.",
    required=True,
    multiple=True,
    type=click.Path(exists=True)
)
@click.option(
    "--output_directory",
    "-o",
    help="For each [input.inet.hdf] file, an [input.evidence.hdf] evidence "
        "file is created. "
        "If no output directory is provided, each [input.evidence.hdf] file is "
        "placed in the same folder as its corresponding [input.inet.hdf] file. "
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
    help="Save the log to a [log.txt] file. "
        "This log file can also be supplied through a [parameters.json] file. "
        "If the log file already exists, the new log data is appended.",
    type=click.Path(dir_okay=False)
)
def evidence(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name
):
    interface.evidence_ion_networks(
        input_path,
        output_directory,
        parameter_file_name,
        log_file_name
    )


# TODO: Implement
@click.command(
    "show",
    help="Show and browse ion-networks and their evidence.",
    short_help="Show and browse ion-networks."
)
@click.option(
    "--input_path",
    "-i",
    "ion_network_file_name",
    help="The ion-network file (.inet.hdf) to show.",
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    "--evidence_file",
    "-e",
    "evidence_file_name",
    help="The corresponding evidence file (.evidence.hdf).",
    type=click.Path(exists=True)
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
    help="Save the log to a [log.txt] file. "
        "This log file can also be supplied through a [parameters.json] file. "
        "If the log file already exists, the new log data is appended.",
    type=click.Path(dir_okay=False)
)
def show(
    ion_network_file_name,
    evidence_file_name,
    parameter_file_name,
    log_file_name
):
    interface.show_ion_network(
        ion_network_file_name,
        evidence_file_name,
        parameter_file_name,
        log_file_name
    )


@click.command(
    "gui",
    help="Graphical user interface to analyse ion-networks.",
)
def gui():
    interface.start_gui()


if __name__ == "__main__":
    cli.add_command(convert)
    cli.add_command(create)
    cli.add_command(evidence)
    cli.add_command(show)
    cli.add_command(gui)
    cli()
