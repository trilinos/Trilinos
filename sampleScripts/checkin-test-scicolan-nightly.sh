#!/bin/bash

#
# This is the script that I used to checkin to Trilinos on my scico-lan
# machine.  You can copy this script and adapt it to your own machine.
#

#
# Allow command-line arguments to pass through to cmake configure!
#

EXTRA_ARGS=$@

#
# Set build options
#

echo "
-DTrilinos_ENABLE_Fortran:BOOL=OFF
-D TPL_ENABLE_Boost:BOOL=ON
-D Boost_INCLUDE_DIRS:PATH=/usr/netpub/boost_1_40_0
-D TPL_ENABLE_Netcdf:BOOL=ON
-D Netcdf_INCLUDE_DIRS:PATH=/projects/seacas/current/include
-D Netcdf_LIBRARY_DIRS:PATH=/projects/seacas/current/lib
-D Trios_ENABLE_XDMF:BOOL=ON
-D TPL_ENABLE_XMDF:BOOL=ON
-D XDMF_INCLUDE_DIRS:PATH=/projects/seacas/current/include
-D XDMF_LIBRARY_DIRS:PATH=/projects/seacas/current/lib
-D TPL_ENABLE_HDF5:BOOL=ON
-D HDF5_INCLUDE_DIRS:PATH=/projects/seacas/current/include
-D HDF5_LIBRARY_DIRS:PATH=/projects/seacas/current/lib
-D EpetraExt_USING_HDF5:BOOL=OFF
-D CMAKE_BUILD_DEBUG:BOOL=ON
" > COMMON.config

echo "
-D MPI_BASE_DIR:PATH=/usr/netpub/mpi/OpenMPI/1.4/64Bit/gnu
" > MPI_DEBUG.config


echo "
-DCMAKE_C_COMPILER:PATH=/usr/bin/gcc
-DCMAKE_CXX_COMPILER:PATH=/usr/bin/g++
" > SERIAL_RELEASE.config

#
# Run the standard checkin testing script with my specializations
#

../Trilinos/checkin-test.py \
--no-eg-git-version-check \
--make-options=-j8 \
--ctest-options=-j4 \
--ctest-timeout=180 \
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS  
