#!/bin/bash
module purge

module use /home/projects/x86-64-haswell/modulefiles
module use /home/projects/sparc/modules
module use /home/projects/x86-64-knl/modulefiles
module load sparc-dev/intel-17.0.0-knl

#module load devpack/openmpi/1.10.0/intel/16.1.056/cuda/none
#module load intel/compilers/16.1.056
#module load openmpi/1.10.0/intel/16.1.056
#module load hdf5/1.8.15/openmpi/1.10.0/intel/16.1.056
#module load boost/1.59.0/openmpi/1.10.0/intel/16.1.056
#module load netcdf-exo/4.3.3.1/openmpi/1.10.0/intel/16.1.056
#module load cgns/3.2.1/openmpi/1.10.0/intel/16.1.056
#module load parmetis/4.0.3/openmpi/1.10.0/intel/16.1.056
#module load superlu-dist/4.3.0/openmpi/1.10.0/intel/16.1.056
#module load pnetcdf/1.6.1/openmpi/1.10.0/intel/16.1.056
#module load totalview/2016.01.06
#module load cmake/3.4.3
#module load intel/vtune/2016.1.1
#module load intel/inspector/2016.1.1
#module load git/20150310
#module load zlib/1.2.8
#module load intel/advisor/2016.1.20
#module load seacas/serial/20160328

unset ATTB_ENV

export PATH=$PATH:$HOME/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/bin

export NETCDF_BASE_DIR=$NETCDF_ROOT
export HDF_BASE_DIR=$HDF5_ROOT
export BOOST_BASE_DIR=$BOOST_ROOT
export LAPACK_ROOT=${OPENBLAS_ROOT}/lib
export BLAS_ROOT=${OPENBLAS_ROOT}/lib
export PATH=$SEACAS_ROOT:$PATH
export MPI_HOME=$MPI_ROOT
