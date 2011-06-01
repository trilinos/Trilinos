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

echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_DEBUG_SYMBOLS:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DTPL_ENABLE_MPI:BOOL=ON
-DTPL_ENABLE_Boost:BOOL=ON
" > MPI_DEBUG_BOOST.config

../../../Trilinos/checkin-test.py \
-j12 \
--ctest-timeout=180 \
--ctest-options="-E '(Ifpack_BlockCheby_MPI_4|Moertel_test1_MPI_2)'" \
$EXTRA_ARGS  

# The above two tests Ifpack_BlockCheby_MPI_4 and Moertel_test1_MPI_2 are
# disabled becuses they have checked STL errors (see bugs 5203 and 5204).
