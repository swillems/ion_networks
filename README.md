# Ion-networks
Analysis of LC-[...]-MSMS data with ion-networks.

## Installation
The ion-networks repository was developed on a Ubuntu 18.04 machine with the [python](https://docs.python.org/3.8/) language and standard python packages. It is likely to function on other systems as well, but this has not been verified.

### Windows 10
For Windows users, the recommended approach is to install a Windows subsystem for Linux (WSL) with Ubuntu 18.04 by following [these steps](https://docs.microsoft.com/en-us/windows/wsl/install-win10). After the WSL has been installed and a user account has been created, close the WSL and download and install [MobaXterm v11.1](https://mobaxterm.mobatek.net/download-home-edition.html) or higher. Then, follow the installation steps for Ubuntu 18.04 within the MobaXterm WSL.

### Ubuntu 18.04
This repository requires python to be run within a [conda](https://conda.io/projects/conda/en/latest/index.html) environment. If not installed, the following commands will perform a full installation of Anaconda3, followed by the installation of the ion-networks repository by [downloading](https://github.com/swillems/ion_networks/tree/master/install/install.sh) and running the install script at the desired location:

```bash
cd /desired/installation/path/
wget https://raw.githubusercontent.com/swillems/ion_networks/master/install/install.sh && bash install.sh && rm install.sh
# wget --auth-no-challenge -O install.sh https://raw.githubusercontent.com/swillems/ion_networks/master/install/install.sh?token=XXX && bash install.sh && rm install.sh
```

If conda is already installed or if the ion-networks repository needs to be updated, [download](https://github.com/swillems/ion_networks/tree/master/install/update.sh) and/or run the ```update.sh``` bash script contained in the [install](https://github.com/swillems/ion_networks/tree/master/install) folder of the ion-networks repository.

## Usage
Five basic modules have been implemented for analysis of ion-networks:

1. Convert
2. Create
3. Evidence
4. Show

These modules can be run either with a GUI or through a CLI.

### Windows 10
For Windows users, a (desktop) shortcut can be created in Windows that automatically runs the ion-networks GUI within the WSL and MobaXterm. This can be done by opening MobaXterm and pressing the *session* button on the top left. Select the rightmost tab *WSL* and set the Linux distribution to Ubuntu in the *Basic WSL settings* tab. Click the *Advanced WSL settings* tab and copy ```ion_networks.py gui``` to the *Execute the following commands at startup* window. Finally, click the *Bookmark settings* tab and change the *Session name* to e.g. *ion_network_gui*. Click the *Create a desktop shortcut to this session* button and select both options *Hide terminal on startup* and *Close MobaXterm on exit* before pressing *OK* in this popup. Confirm the session settings with *OK*. A pop-up with the GUI running should have appeared in your taskbar, allowing you to test the installation. For subsequent use, double-clicking the Windows desktop icon suffices to run the ion-networks GUI.

### Ubuntu 18.04
The ion-networks software can be run within a terminal with the command ```ion_networks.py``` (this alias is set by default during installation). Possible commands are:

```bash
Usage: ion_networks.py [OPTIONS] COMMAND [ARGS]...

  Analysis of LC-[...]-MSMS data with ion-networks.

Options:
  -h, --help  Show this message and exit.

Commands:
  convert   Convert various input formats to unified input.
  create    Create ion-networks from unified input.
  evidence  Collect evidence for ion-networks.
  gui       Graphical user interface to analyse ion-networks.
  show      Show and browse ion-networks.
```

Each command then comes with its own help function through ```ion_networks.py COMMAND -h```.

Typically, a workflow looks as follows:

```bash
ion_networks.py convert -i project_folder/centroided_data/ -o project_folder/ion_networks/ -hdmse -l project_folder/ion_networks/log.txt
ion_networks.py create -i project_folder/ion_networks/ -l project_folder/ion_networks/log.txt
ion_networks.py evidence -i project_folder/ion_networks/ -l project_folder/ion_networks/log.txt
ion_networks.py show -i project_folder/ion_networks/ -l project_folder/ion_networks/log.txt
```

Alternatively, a GUI can be used by running the command ```ion_networks.py gui```.
