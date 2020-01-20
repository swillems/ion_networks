#!python

import click

import network


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    context_settings=CONTEXT_SETTINGS,
    help="Analysis of LC-IMS-MSMS data with ion-networks."
)
def cli():
    pass


@click.command(
    "create",
    help="Create ion-networks.",
    # invoke_without_command=False
)
@click.option(
    "--raw_file",
    "-r",
    help="The raw input file (.csv or .hdf) with centroided ion peaks. "
        "This flag can be set multiple times to create multiple ion-networks.",
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
@click.option(
    "--parameters",
    "-p",
    help="A parameter file (.json).",
    type=click.File('r')
)
def create(
    raw_file,
    ion_network_file,
    parameters
):
    for index, raw_file_name in enumerate(raw_file):
        if len(ion_network_file) > index:
            network.Network(
                raw_file_name=raw_file_name,
                network_file_name=ion_network_file[index],
                parameters=parameters
            )
        else:
            network.Network(
                raw_file_name=raw_file_name,
                parameters=parameters
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
    "--parameters",
    "-p",
    help="A parameters file (.json).",
    type=click.File('r')
)
def align(
    ion_network_file,
    alignment_file,
    parameters
):
    for index, first_ion_network_file in enumerate(ion_network_file[:-1]):
        first_ion_network = network.Network(
            network_file_name=first_ion_network_file,
            parameters=parameters
        )
        for second_ion_network_file in ion_network_file[index + 1:]:
            second_ion_network = network.Network(
                network_file_name=second_ion_network_file,
                parameters=parameters
            )
            first_ion_network.align(
                second_ion_network,
                alignment_file_name=alignment_file,
                parameters=parameters
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
    "--parameters",
    "-p",
    help="A parameters file (.json).",
    type=click.File('r')
)
def evidence(
    ion_network_file,
    alignment_file,
    evidence_file,
    parameters
):
    for index, ion_network_file in enumerate(ion_network_file):
        ion_network = network.Network(
            network_file_name=ion_network_file,
            parameters=parameters
        )
        if len(evidence_file) > index:
            ion_network.evidence(
                alignment_file,
                evidence_file_name=evidence_file[index],
                parameters=parameters
            )
        else:
            ion_network.evidence(
                alignment_file,
                parameters=parameters
            )


@click.command(
    "show",
    help="Show ion-networks."
)
@click.option(
    "--ion_network_file",
    "-i",
    help="The ion-network file (.hdf) to show.",
    required=True,
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
    "--parameters",
    "-p",
    help="A parameters file (.json).",
    type=click.File('r')
)
def show(
    ion_network_file,
    evidence_file,
    parameters
):
    pass


@click.command(
    "gui",
    help="Graphical user interface for ion-networks."
)
def gui():
    pass


if __name__ == "__main__":
    cli.add_command(create)
    cli.add_command(align)
    cli.add_command(evidence)
    cli.add_command(show)
    cli.add_command(gui)
    cli()
