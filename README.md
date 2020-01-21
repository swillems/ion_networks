# Ion-networks
Analysis of LC-IMS-MSMS data with ion-networks.

## Installation
This repository requires [conda](https://conda.io/projects/conda/en/latest/index.html). After installation of conda, run the following commands in a terminal to download this repository in the desired location:

```bash
cd /location/where/to/install
git clone https://github.com/swillems/ion_networks.git
conda env create -f ion_networks/install/environment.yml
```

## Usage

In a terminal, run the following command:

```bash
cd /location/where/to/install/ion_networks/src
conda activate ion_networks
python main.py
```

Options are:
```bash
Usage: main.py [OPTIONS] COMMAND [ARGS]...

  Analysis of LC-IMS-MSMS data with ion-networks.

Options:
  -h, --help  Show this message and exit.

Commands:
  align     Align ion-networks.
  create    Create ion-networks.
  evidence  Evidence ion-networks.
  gui       Graphical user interface for ion-networks.
  show      Show ion-networks.
```

Command options are:
```bash
Usage: main.py create [OPTIONS]

  Create ion-networks.

Options:
  -r, --raw_file FILENAME         The raw input file (.csv or .hdf) with
                                  centroided ion peaks. This flag can be set
                                  multiple times to create multiple ion-
                                  networks.  [required]
  -i, --ion_network_file FILENAME
                                  A new ion-network file (.hdf). This flag can
                                  be set multiple times to save multiple ion-
                                  networks in the same order. If not set, a
                                  new '[raw_file].hdf' file will be created
                                  per raw file. WARNING: This overrides
                                  already existing files without confirmation.
  -p, --parameters FILENAME       A parameter file (.json).
  -h, --help                      Show this message and exit.
```

```bash
Usage: main.py align [OPTIONS]

  Align ion-networks.

Options:
  -i, --ion_network_file FILENAME
                                  The ion-network file (.hdf) to align. This
                                  flag can be set multiple times to align
                                  multiple ion-networks pairwise.  [required]
  -a, --alignment_file FILENAME   A new alignment file (.hdf) with all
                                  pairwise alignments. If not set, an
                                  'alignment.hdf' file will be created in
                                  directory of the first [ion_network_file].
                                  WARNING: This overrides already existing
                                  files without confirmation.
  -p, --parameters FILENAME       A parameters file (.json).
  -h, --help                      Show this message and exit.
```

```bash
Usage: main.py evidence [OPTIONS]

  Evidence ion-networks.

Options:
  -i, --ion_network_file FILENAME
                                  The ion-network file (.hdf) to evidence.This
                                  flag can be set multiple times to evidence
                                  multiple ion-networks.  [required]
  -a, --alignment_file FILENAME   The alignment file (.hdf) from where to get
                                  the evidence. If a single ion-network was
                                  provided, evidence is drawn from all ion-
                                  networks present in the alignment file. If
                                  multiple ion-networks are provided that are
                                  present in the alignment file, only those
                                  will be used as evidence for eachother.
                                  [required]
  -e, --evidence_file FILENAME    A new evidence file (.hdf) for the ion-
                                  network. If not set, an
                                  '[ion_network_file].evidence.hdf' file will
                                  be created per ion-network. WARNING: This
                                  overrides already existing files without
                                  confirmation.
  -p, --parameters FILENAME       A parameters file (.json).
  -h, --help                      Show this message and exit.
```
