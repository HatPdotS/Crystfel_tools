Load anaconda and copy and run the script below

The script takes to arguments:

1. the path to where you want to create the conda environment
2. the path where you want to clone the python source code into

Example scripts will be setup inside the diectory you run this from in a folder example scripts

'''

#!/bin/bash


CONDAPATH=$1 # name of the conda environemt

INSTALLPATH=$2 # directory where custom packages will be installed

conda create -p $CONDAPATH # creates the conda environment

conda activate $CONDANAME # acitvate conda environment

conda install python=3.11

git clone https://github.com/HatPdotS/Crystfel_tools.git $INSTALLPATH/crystfel_tools # clones crystfel tool from peters github 
git clone https://github.com/HatPdotS/Slurm-tools.git  $INSTALLPATH/slurm_tools # clones slurm tools from peters github

pip install pandas h5py numpy tabulate # installs other dependencies

pip install -e $INSTALLPATH/slurm_tools # installs slurmtools
pip install -e $INSTALLPATH/crystfel_tools # installs crystfeltools
cp -r $INSTALLPATH/crystfel_tools/example_scripts .
'''