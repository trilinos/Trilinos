#!/bin/bash

module purge

case $TARGET_COMPILER in
    gnu)
        module load gnu/4.7.0
        module load openmpi-gnu/1.4
        ;;
    pgi)
        module load pgi/12.5
        module load openmpi-pgi/1.4
        ;;
    intel | *)
        module load intel/12.1
        module load openmpi-intel/1.4
        ;;
esac

module list
