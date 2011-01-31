#!/bin/bash

#
# This is the script that I use to do checkin testing with GCC and OpenMPI
# installed with the install-gcc.py and install-openmpi.py scripts.
#

EXTRA_ARGS=$@

TRILINOS_TOOLSET_BASE=/usr/local/trilinos-toolset

echo "
-DTrilinos_EXTRA_LINK_FLAGS:STRING='-Wl,-rpath,$TRILINOS_TOOLSET_BASE/lib64'
" > COMMON.config

echo "
-DMPI_BASE_DIR:PATH=$TRILINOS_TOOLSET_BASE
" > MPI_DEBUG.config

echo "
-DCMAKE_CXX_COMPILER:PATH=$TRILINOS_TOOLSET_BASE/bin/g++
-DCMAKE_C_COMPILER:PATH=$TRILINOS_TOOLSET_BASE/bin/gcc
" > SERIAL_RELEASE.config

../../../Trilinos/checkin-test.py \
-j12 \
--ctest-timeout=180 \
$EXTRA_ARGS  
