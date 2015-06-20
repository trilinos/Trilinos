#!/bin/bash

module purge

case $TARGET_COMPILER in
    gnu)
        module load gnu/4.7.1
        module load openmpi-gnu/1.6
        ;;
    pgi)
        module load compilers/pgi-12.8
        module load mpi/openmpi-1.4.3_oobpr_pgi-12.8-0
        ;;
    intel | *)
        module load compilers/intel-11.1-f064-c064
        module load mpi/openmpi-1.4.3_oobpr_intel-11.1-f064-c064
        ;;
esac

module list


