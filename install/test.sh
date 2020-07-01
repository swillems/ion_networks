#!bash

echo "Testing installation"
mkdir ion_networks/projects
tar xf ion_networks/install/test_data.tar.xz --directory ion_networks/projects
eval "$(conda shell.bash hook)"
ion_networks_command="$(conda activate ion_networks; which python)"
ion_networks='"${ion_networks_command}" "$(pwd)"/ion_networks/src/ion_networks.py'
"$ion_networks" create -i ion_networks/projects/test_data -l ion_networks/projects/test_data
"$ion_networks" evidence -i ion_networks/projects/test_data -l ion_networks/projects/test_data
"$ion_networks" show
