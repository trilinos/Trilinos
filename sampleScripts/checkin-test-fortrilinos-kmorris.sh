#!/bin/bash

#
# This is the script that I used to checkin to Trilinos on
# trilinos-test2.sandia.gov.  You can copy this script and adapt it to
# your own use on this machine.
#
# NOTE: You need to add '/home/trilinos/bin' to your path before you
# run this script.
#

EXTRA_ARGS=$@

#
# Set up basic environment options
#

#echo "-DBUILD_SHARED_LIBS:BOOL=ON" > COMMON.config

echo "
-DMPI_BASE_DIR:PATH=/usr/local/openmpi/nagfor
-DMPI_BIN_DIR:PATH=/usr/local/openmpi/nagfor/bin
-DMPI_INCLUDE_PATH:PATH=/usr/local/openmpi/nagfor/include
-DMPI_Fortran_COMPILER:FILEPATH=/usr/local/openmpi/nagfor/bin/mpif90
-DMPI_CXX_COMPILER:FILEPATH=/usr/local/openmpi/nagfor/bin/mpicxx
-DMPI_C_COMPILER:FILEPATH=/usr/local/openmpi/nagfor/bin/mpicc
"  > MPI_DEBUG.config

echo "
-DCMAKE_CXX_COMPILER:FILEPATH=g++
-DCMAKE_C_COMPILER:FILEPATH=gcc
-DCMAKE_Fortran_COMPILER:FILEPATH=nagfor
" > SERIAL_RELEASE.config

echo " 
-D CMAKE_Fortran_FLAGS:STRING='-f2003 -g -C=all -wmismatch=mpi_scatterv,mpi_alltoallv,mpi_gatherv,mpi_allgatherv,mpi_bcast '
-D TPL_ENABLE_MPI:BOOL=ON 
-D MPI_BASE_DIR:PATH=/usr/local/openmpi/nagfor
-D MPI_BIN_DIR:PATH=/usr/local/openmpi/nagfor/bin
-D MPI_INCLUDE_PATH:PATH=/usr/local/openmpi/nagfor/include
-D MPI_USE_COMPILER_WRAPPERS:BOOL=ON 
-D MPI_Fortran_COMPILER:FILEPATH=/usr/local/openmpi/nagfor/bin/mpif90
-D MPI_CXX_COMPILER:FILEPATH=/usr/local/openmpi/nagfor/bin/mpicxx
-D MPI_C_COMPILER:FILEPATH=/usr/local/openmpi/nagfor/bin/mpicc
-D HAVE_GCC_ABI_DEMANGLE:BOOL=ON 
-D ForTrilinos_ENABLE_OBJECT_ORIENTED:BOOL=ON 
" > ForTrilinos_MPI.config

echo "
-DCMAKE_CXX_COMPILER:FILEPATH=g++
-DCMAKE_C_COMPILER:FILEPATH=gcc
-DCMAKE_Fortran_COMPILER:FILEPATH=nagfor
-DCMAKE_Fortran_FLAGS:STRING='-f2003 -g -C=all'
-DForTrilinos_ENABLE_OBJECT_ORIENTED:BOOL=ON 
" > ForTrilinos_SERIAL.config

#
# Run the checkin-test.py script with more arguments
#

../../Trilinos/checkin-test.py \
--do-all \
--make-options="-j6" \
--ctest-options="-j8" \
--ctest-timeout=600 \
--enable-packages=ForTrilinos \
--enable-all-packages=off \
--enable-fwd-packages \
--no-eg-git-version-check \
--extra-builds=ForTrilinos_MPI,ForTrilinos_Serial \
$EXTRA_ARGS  
