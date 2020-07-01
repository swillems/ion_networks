#!bash

conda remove --name ion_networks --all
sed '/^alias ion_networks/d' -i ~/.bashrc
