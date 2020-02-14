#!bash

if ! hash conda 2>/dev/null; then
  echo "Conda not found."
  echo "Downloading conda."
  wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
  # wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  echo "Installing conda."
  # bash Anaconda3-2019.10-Linux-x86_64.sh -b -p ~/anaconda3
  bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
  echo "Updating conda."
  conda update -n root conda
  eval "$(~/anaconda3/bin/conda bash hook)"
  conda init
  conda config --set auto_activate_base false
  conda activate
  conda install -c conda-forge nb_conda_kernels
  conda deactivate
  # echo 'function conda_notebooks() { conda activate && jupyter_notebook && conda deactivate; }' >> ~/.bashrc
  # apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
else
  echo "Conda already installed."
fi

if ! hash ion_networks.py 2>/dev/null; then
  echo "Downloading latest ion-networks repository from github."
  git clone https://github.com/swillems/ion_networks.git
  echo "Installing ion-networks."
  conda env create --file ion_networks/install/environment.yml
  # sed -i '/function ion_networks.py() { conda activate ion_networks && python .* "$@" && conda deactivate; }/d' ~/.bashrc
  echo "Adding ion-networks.py command to ~/.bashrc."
  # echo 'function ion_networks.py() { conda activate ion_networks && python '$(pwd)'/ion_networks/src/ion_networks.py "$@" && conda deactivate; }' >> ~/.bashrc
  conda activate ion_networks
  ion_networks_command="$(which python)"
  conda deactivate
  echo "alias ion_networks.py='"${ion_networks_command}" "$(pwd)"/ion_networks/src/ion_networks.py'" >> ~/.bashrc
  source ~/.bashrc
else
  echo "Ion-networks are already installed."
  echo "Update with command: bash ion_networks install/update.sh"
fi
