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
fi
