#!/bin/bash

# Used to test Trilinos on u233.ornl.gov

# NOTE: To use this, you must first prepend /opt/trilinos-toolset/bin
# to your path to find eg and cmake!

EXTRA_ARGS=$@

TRILINOS_TOOLSET_BASE=/opt/gcc-4.5.1/trilinos-toolset

echo "
-DTrilinos_EXTRA_LINK_FLAGS:STRING='-Wl,-rpath,$TRILINOS_TOOLSET_BASE/lib64'
-DTPL_BLAS_LIBRARIES=/usr/lib64/libblas.so.3
-DTPL_LAPACK_LIBRARIES=/usr/lib64/liblapack.so.3
" > COMMON.config

echo "
-DMPI_BASE_DIR:PATH=$TRILINOS_TOOLSET_BASE
" > MPI_DEBUG.config

echo "
-DCMAKE_CXX_COMPILER:PATH=$TRILINOS_TOOLSET_BASE/bin/g++
-DCMAKE_C_COMPILER:PATH=$TRILINOS_TOOLSET_BASE/bin/gcc
" > SERIAL_RELEASE.config

../../../Trilinos/checkin-test.py \
-j16 \
--ctest-timeout=180 \
$EXTRA_ARGS  

# NOTE: By default we use 16 processes which is 1/2 of the 32
# processes on this machine.  This way two people can build and test
# Trilinos without taxing the machine too much.
