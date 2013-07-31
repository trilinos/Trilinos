#!/bin/bash

module purge

case $TARGET_COMPILER in
    gnu)
        module load compilers/gcc-4.6.2
        module load mpi/openmpi-1.4.3_oobpr_gcc-4.6.2
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

module load libraries/papi-3.7.2

module list


