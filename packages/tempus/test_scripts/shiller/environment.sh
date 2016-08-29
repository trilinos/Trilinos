#!/bin/bash
module purge
module load devpack/openmpi/1.10.0/gcc/4.8.4/cuda/7.5.18

unset ATTB_ENV

export PATH=$PATH:$HOME/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/bin

export NETCDF_BASE_DIR=$NETCDF_ROOT
export HDF_BASE_DIR=$HDF5_ROOT
export BOOST_BASE_DIR=$BOOST_ROOT
