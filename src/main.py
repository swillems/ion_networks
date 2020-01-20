#!python

import click
import logging


@click.group()
def cli():
    """
    Analysis of LC-IMS-MSMS data with ion-networks.

    """
    pass


@click.command(
    'create',
    short_help='Create ion-networks'
)
def create():
    pass


@click.command(
    'align',
    short_help='Align ion-networks'
)
def align():
    pass


@click.command(
    'evidence',
    short_help='Evidence ion-networks'
)
def evidence():
    pass


@click.command(
    'show',
    short_help='Show ion-networks'
)
def show():
    pass


@click.command(
    'gui',
    short_help='GUI for ion-networks'
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
