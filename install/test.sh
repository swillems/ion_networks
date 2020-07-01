#!bash

echo "Testing installation"
mkdir ion_networks/projects
tar xf ion_networks/install/test_data.tar.xz --directory projects
ion_networks.py create -i ion_networks/projects/test_data -l ion_networks/projects/test_data
ion_networks.py evidence -i ion_networks/projects/test_data -l ion_networks/projects/test_data
ion_networks.py show
