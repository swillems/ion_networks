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
  echo "Downloading latest ion-networks repository from GitHub."
  git clone https://github.com/swillems/ion_networks.git
  echo "Installing ion-networks."
  conda env create --file ion_networks/install/environment.yml
  echo "Adding ion-networks.py alias to ~/.bashrc."
  ion_networks_command="$(conda activate; which python)"
  echo "alias ion_networks.py='"${ion_networks_command}" "$(pwd)"/ion_networks/src/ion_networks.py'" >> ~/.bashrc
  # export ion_networks.py='"${ion_networks_command}" "$(pwd)"/ion_networks/src/ion_networks.py'
else
  echo "Ion-networks are already installed."
  echo "Update with command 'bash ion_networks/install/update.sh'."
fi
