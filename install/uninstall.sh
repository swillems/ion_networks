#!bash

conda remove --name ion_networks --all -y
if [ -n "$ZSH_VERSION" ]; then
   sed '/^alias ion_networks/d' -i ~/.zshrc
elif [ -n "$BASH_VERSION" ]; then
   sed '/^alias ion_networks/d' -i ~/.bashrc
else
   echo "Unknown shell."
fi
