#!bash

if ! hash conda 2>/dev/null; then
  echo "Conda not found."
  echo "Downloading conda."
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  echo "Installing conda."
  bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
  echo "Updating conda."
  eval "$(~/miniconda3/bin/conda shell.bash hook)"
  conda update -n root conda -y
  echo "Initializing conda."
  conda init
  conda config --set auto_activate_base false
  echo "Cleaning up conda installation."
  source ~/.bashrc
  rm Miniconda3-latest-Linux-x86_64.sh
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
  ion_networks_command="$(conda activate ion_networks; which python)"
  cd ion_networks
  ion_networks_command "/ion_networks/setup.py" install
  cd ..
  if [ -n "$ZSH_VERSION" ]; then
     echo "Adding ion-networks.py alias to ~/.zshrc."
     echo "alias ion_networks.py='"${ion_networks_command}" "$(pwd)"/ion_networks/src/ion_networks.py'" >> ~/.zshrc
  elif [ -n "$BASH_VERSION" ]; then
     echo "Adding ion-networks.py alias to ~/.bashrc."
     echo "alias ion_networks.py='"${ion_networks_command}" "$(pwd)"/ion_networks/src/ion_networks.py'" >> ~/.bashrc
  else
     echo "Unknown shell."
  fi
  # export ion_networks.py='"${ion_networks_command}" "$(pwd)"/ion_networks/src/ion_networks.py'
else
  echo "Ion-networks are already installed."
  echo "Update with command 'bash ion_networks/install/update.sh'."
fi
