#!bash

echo "Testing installation"
mkdir ion_networks/projects
tar xf ion_networks/install/test_data.tar.xz --directory ion_networks/projects
eval "$(conda shell.bash hook)"
ion_networks_bin="$(conda activate ion_networks; which ion_networks)"
eval "$ion_networks_bin" create -i ion_networks/projects/test_data -l ion_networks/projects/test_data
eval "$ion_networks_bin" evidence -i ion_networks/projects/test_data -l ion_networks/projects/test_data
eval "$ion_networks_bin" show
