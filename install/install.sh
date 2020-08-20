#!sh

if ! hash conda 2>/dev/null; then
  echo "Conda not found."
  echo "Downloading conda."
  if [ "$(uname)" == "Darwin" ]; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    echo "Installing conda."
    sh Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
    echo "Updating conda."
    eval "$(~/miniconda3/bin/conda shell.bash hook)"
    conda update -n root conda -y
    echo "Initializing conda."
    conda init
    conda config --set auto_activate_base false
    echo "Cleaning up conda installation."
    source ~/.bashrc
    rm Miniconda3-latest-Linux-x86_64.sh
  elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    echo "Installing conda."
    sh Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
    echo "Updating conda."
    eval "$(~/miniconda3/bin/conda shell.bash hook)"
    conda update -n root conda -y
    echo "Initializing conda."
    conda init
    conda config --set auto_activate_base false
    echo "Cleaning up conda installation."
    source ~/.bashrc
    rm Miniconda3-latest-MacOSX-x86_64.sh
  else
    echo "Detected unknown OS"
  fi
else
  echo "Conda is already installed."
fi

if ! hash ion_networks.py 2>/dev/null; then
  DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
  if ! [[ "$DIR" == */ion_networks/install ]]; then
    echo "Downloading latest ion-networks repository from GitHub."
    git clone https://github.com/swillems/ion_networks.git
  else
    echo "Source already downloaded"
  fi
  echo "Installing ion-networks."
  conda env create --file ion_networks/install/environment.yml
  eval "$(conda shell.bash hook)"
  pip_command="$(conda activate ion_networks; which pip)"
  if [ "$(uname)" == "Darwin" ]; then
    echo "Detected Mac OS X"
    git clone https://github.com/compomics/ms2pip_c.git
    sed -i .bak '/-fno-var-tracking-assignments/d' ms2pip_c/setup.py
    "${pip_command}" install cython
    "${pip_command}" install './ms2pip_c'
  elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    "${pip_command}" install ms2pip
  else
    echo "Detected unknown OS"
  fi
  "${pip_command}" install './ion_networks[ms2pip]'
  ion_networks_bin="$(conda activate ion_networks; which ion_networks)"
  if [ -n "$ZSH_VERSION" ]; then
     echo "Adding ion-networks alias to ~/.zshrc."
     echo "alias ion_networks='"${ion_networks_bin}"'" >> ~/.zshrc
  elif [ -n "$BASH_VERSION" ]; then
     echo "Adding ion-networks alias to ~/.bashrc."
     echo "alias ion_networks='"${ion_networks_bin}"'" >> ~/.bashrc
  else
     echo "Unknown shell."
  fi
else
  echo "Ion-networks are already installed."
  echo "Update with command 'bash ion_networks/install/update.sh'."
fi
