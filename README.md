# MD_Scripts
Various scripts for molecular dynamics simulations.

When running MD simulations one often has to do several small tasks that get repetitive and can be trivially automated.
This repository will store some scripts that I write during my work that can perhaps be used by others easily. 

## Scripts
Here I will briefly review each of the scripts as well as where to find them.
#### General
* box_calculator.py will return the side length of a cubic box of a given density of atoms that you gave input. It will
do the same for a spherical box

#### LAMMPS
* LAMMPS_Install.sh is bash script that allows one to install LAMMPS using cmake with most of its packages. Some of the user
packages are missing but can be very easily added. If you have any troubles make an issue and I can add it for you. Note that
this script is designed to be run on a cluster with a package manager. If you are installing on a local system, the module 
load commands at the start of the script may not be necessary, or can be converted into paths for compilers and other libraries. 

#### VASP
* generate_vasp_potcar.py is a python script that takes several potcar files and returns a single, concatenated file. This will
avoid any errors that might come out of an OS interpreting the format of the files poorly.
* generate_vasp_poscar.py is a python script that will collect an xyz coordinate file and build a single poscar file with
given information. If using this along with the potcar generator, remember to add the files in the correct order

## Running the scripts
These scripts are written in several languages depending on what fit the best at the time. This is a general overview of 
how they can be run.

#### Python scripts
The python scripts all assume python3 and can be run with
```
python script.py
```
For information about what the script does, one can usually run
```
python script.py -h
```
A notable exception is the box calculator which uses inputs rather than argparse. In this case, by simply running the script, all
necessary input will be asked for.

#### Bash scripts
Bash scripts can be called either with the bash command:
```
bash script.sh
```
or they can be made into an executable with
```
chmod +x script.sh
```
which can then be run with the usual
```
./script.sh
```

