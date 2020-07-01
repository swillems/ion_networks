#!bash

echo "Testing installation"
mkdir ion_networks/projects
tar xf install/test_data.tar.xz --directory projects
ion_networks.py create -i projects/test_data -l projects/test_data
ion_networks.py evidence -i projects/test_data -l projects/test_data
ion_networks.py show
