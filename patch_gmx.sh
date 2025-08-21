#!/bin/bash
PTCH_ROOT=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
DEST_DIR="gromacs/gmxana"

NEMODULES=''
for CPP in $(ls "${PTCH_ROOT}"/*.cpp)
do
    CPP0=$(basename "$CPP")
    MODULES+="${CPP0/.cpp} "
done


if [[ -z $1 ]]; then
echo "This script will patch your GROMACS source directory to add following gmx commands:"
echo "$MODULES"
echo "To use the patch provide the path to GORMACS \"src\" directory."
echo "Rebuild GROMACS after applying patch."

else
    GMX_ROOT=$1
    echo "Patching files to make new modules as functional gmx commands:"
    for P in $(ls ${PTCH_ROOT}/*.patch)
    do
        P0=$(basename $P)
        F0=${P0/.patch}
        F=$(find ${GMX_ROOT} -name ${F0})
        if [[ (-z ${F}) || ( $(echo ${F} | wc -l ) -gt 1 ) ]]; then
            echo "ERROR: file ${F0} not foung in ${GMX_ROOT} folder or more than one copy was found!"
            echo "${GMX_ROOT} should be pointing to the GROMACS/src directory"
        else
            patch $F $P
        fi
    done

    for M in ${MODULES}
    do
	echo "Copying ${M}.cpp to ${GMX_ROOT}/${DEST_DIR}" 
        cp ${PTCH_ROOT}/${M}.* ${GMX_ROOT}/${DEST_DIR}
    done

fi
