#!bash

conda remove --name ion_networks --all -y
sed '/^alias ion_networks/d' -i ~/.bashrc
