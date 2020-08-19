#!bash

conda remove --name ion_networks --all -y
if test -f "~/.zshrc"; then
   sed '/^alias ion_networks/d' -i ~/.zshrc
fi
if test -f "~/.bashrc"; then
   sed '/^alias ion_networks/d' -i ~/.bashrc
fi
