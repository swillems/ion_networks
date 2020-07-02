#!bash

echo "Testing installation"
mkdir ion_networks/projects
tar xf ion_networks/install/test_data.tar.xz --directory ion_networks/projects
eval "$(conda shell.bash hook)"
ion_networks_command="$(conda activate ion_networks; which python)"
ion_networks='"${ion_networks_command}" "$(pwd)"/ion_networks/src/ion_networks.py'
eval "$ion_networks" create -i ion_networks/projects/test_data -l ion_networks/projects/test_data
eval "$ion_networks" evidence -i ion_networks/projects/test_data -l ion_networks/projects/test_data
cd ion_networks/projects
eval "$ion_networks" show
cd ../..
