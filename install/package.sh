python setup.py sdist
pip install twine
twine upload dist/*

conda create -n ion_networks_installer python=3.6
conda activate ion_networks_installer
pip install pypiwin32
pip install pyinstaller
pip install ion_networks
pip install matplotlib==3.2.2
pyinstaller package_ion_networks.py --noconfirm --add-data ~/miniconda3/envs/ion_networks/lib/python3.6/site-packages/ion_networks/lib:ion_networks/lib --windowed --onefile -n ion_networks
conda deactivate
