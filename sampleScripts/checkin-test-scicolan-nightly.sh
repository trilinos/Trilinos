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
-D Boost_INCLUDE_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/code/TPLs_src/boost
-D Boost_LIBRARY_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/build/SierraBuild.trilinos_dev.devsierra/boost/1.46.1/gcc-4.4.4/release/address-model-64/mpi-openmpi/runtime-link-shared
-D TPL_ENABLE_gtest:BOOL=ON
-D gtest_INCLUDE_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/code/TPLs_src/gtest/include
-D gtest_LIBRARY_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/build/SierraBuild.trilinos_dev.devsierra/gtest/1.5.0/gcc-4.4.4/release/address-model-64/mpi-openmpi/runtime-link-shared/
-D TPL_ENABLE_Netcdf:BOOL=ON
-D TPL_Netcdf_INCLUDE_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/code/TPLs_src/netcdf/src/include
-D Netcdf_LIBRARY_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/build/SierraBuild.trilinos_dev.devsierra/netcdf/4.0-snl3/gcc-4.4.4/release/address-model-64/mpi-openmpi/runtime-link-shared
-D Trios_ENABLE_XDMF:BOOL=OFF
-D TPL_ENABLE_XMDF:BOOL=OFF
-D TPL_ENABLE_HDF5:BOOL=OFF
-D Claps_ENABLE_SClaps:BOOL=ON
" > COMMON.config
#-D CMAKE_BUILD_DEBUG:BOOL=ON
#-D Netcdf_INCLUDE_DIRS:PATH=/projects/seacas/current/include
#-D TPL_ENABLE_BoostLib:BOOL=ON
#-D BoostLib_INCLUDE_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/code/TPLs_src/boost
#-D BoostLib_LIBRARY_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/build/SierraBuild.trilinos_dev.devsierra/boost/1.46.1/gcc-4.4.4/release/address-model-64/mpi-openmpi/runtime-link-shared

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
--no-rebase \
$EXTRA_ARGS  
