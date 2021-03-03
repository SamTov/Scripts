#!/bin/bash

# A script to download and install the LAMMPS simulation package
# with most of the standard addons. None of the additions require
# additional downloads.

# ################## #
# Import BEE modules #
# ################## #

module purge
module load gcc/8.4.0
module load openmpi/4.0.5
module load openblas/0.3.13
module load python/3.8.7
module load py-pip/20.2
module load cmake/3.19.2
module load fftw/3.3.9

# ######################## #
# Define package variables #
# ######################## #

# Standard Packages
asphere="-D PKG_ASPHERE=on"
body="-D PKG_BODY=on"
class2="-D PKG_CLASS2=on"
colloid="-D PKG_COLLOID=on"
compress="-D PKG_COMPRESS=on"
coreshell="-D PKG_CORESHELL=on"
dipole="-D PKG_DIPOLE=on"
gou="-D PKG_GOU=on"
granular="-D PKG_GRANULAR=on"
kim="-D PKG_KIM=on"
kokkos="-D PKG_KOKKOS=on"
kspace="-D PKG_KSPACE=on"
latte="-D PKG_LATTE=on"
manybody="-D PKG_MANYBODY=on"
mc="-D PKG_MC=on"
message="-D PKG_MESSAGE=on"
misc="-D PKG_MISC=on"
mliap="-D PKG_MLIAP=on"
molecule="-D PKG_MOLECULE=on"
mpiio="-D PKG_MPIIO=on"
mscg="-D PKG_MSCG=on"
opt="-D PKG_OPT=on"
peri="-D PKG_PERI=on"
poems="-D PKG_POEMS=on"
python="-D BUILD_SHARED_LIBS=on -D LAMMPS_EXCEPTIONS=on -D PKG_PYTHON=on"
qeq="-D PKG_QEQ=on"
replica="-D PKG_REPLICA=on"
rigid="-D PKG_RIGID=on"
shock="-D PKG_SHOCK=on"
snap="-D PKG_SNAP=on"
spin="-D PKG_SPIN=on"
srd="-D PKG_SRD=on"
voronoi="-D PKG_VORONOI=on"

# User Packages
drude="-D PKG_USER-DRUDE=on"

pkg_list=(
	${asphere}
	${body}
	${class2}
	${colloid}
	${compress}
	${coreshell}
	${dipole}
	${gou}
	${granular}
	${kim}
	${kokkos}
	${kspace}
	${latte}
	${manybody}
	${mc}
	${message}
	${misc}
	${mliap}
	${molecule}
	${mpiio}
	${opt}
	${peri}
	${poems}
	${python}
	${qeq}
	${replica}
	${rigid}
	${shock}
	${snap}
	${spin}
	${srd}
	${voronoi}
	${drude}
)

# #################### #
# Function definitions #
# #################### #
function build_package_string () {
	pkg_string=""
	for item in "${pkg_list[@]}"; do
		pkg_string="${pkg_string} ${item}"
	done
}
build_package_string

# ################## #
# Begin installation #
# ################## #

echo Beginning installation of LAMMPS

# Clone the repository
git clone -b unstable https://github.com/lammps/lammps.git lammps ; cd lammps

mkdir build ; cd build

cmake -C ../cmake/presets/minimal.cmake ${pkg_string} ../cmake
cmake --build .
cmake --install .

echo Build finished


