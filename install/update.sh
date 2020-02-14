#!bash

if ! hash ion_networks.py 2>/dev/null; then
  DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
  cd $DIR
  cd ..
  echo "Ion-networks are already installed, updating from GitHub."
  git stash
  git pull
  conda env update --file install/environment.yml
else
  echo "Ion-networks are not installed yet."
  echo "Try running 'source install/install.sh'."
  # echo "Downloading ion-networks from github."
  # git clone https://github.com/swillems/ion_networks.git
  # echo "Installing ion-networks."
  # conda env create --file ion_networks/install/environment.yml
  # # sed -i '/function ion_networks.py() { conda activate ion_networks && python .* "$@" && conda deactivate; }/d' ~/.bashrc
  # echo 'function ion_networks.py() { conda activate ion_networks && python '$(pwd)'/ion_networks/src/ion_networks.py "$@" && conda deactivate; }' >> ~/.bashrc
fi
