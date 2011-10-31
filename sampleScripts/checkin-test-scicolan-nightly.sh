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
" > COMMON.config

echo "
" > MPI_DEBUG.config


echo "
" > SERIAL_RELEASE.config

echo "
-D Trilinos_ENABLE_TriKota:BOOL=OFF
-D Trilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON
-D TPL_ENABLE_Boost:BOOL=ON
-D Boost_INCLUDE_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/code/TPLs_src/boost
-D Boost_LIBRARY_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/build/SierraBuild.trilinos_dev.devsierra/boost/1.46.1/gcc-4.4.4/release/address-model-64/mpi-openmpi/runtime-link-shared
-D TPL_ENABLE_gtest:BOOL=ON
-D gtest_INCLUDE_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/code/TPLs_src/gtest/include
-D gtest_LIBRARY_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/build/SierraBuild.trilinos_dev.devsierra/gtest/1.5.0/gcc-4.4.4/release/address-model-64/mpi-openmpi/runtime-link-shared/
-D TPL_ENABLE_Netcdf:BOOL=ON
-D TPL_Netcdf_INCLUDE_DIRS:PATH=/projects/seacas/linux_rhel5/current/include
-D Netcdf_LIBRARY_DIRS:PATH=/projects/seacas/linux_rhel5/current/lib
-D SEACAS_ENABLE_XDMF:BOOL=OFF
-D Trilinos_ENABLE_SEACASMat2exo:BOOL=OFF
-D TPL_ENABLE_XMDF:BOOL=OFF
-D TPL_ENABLE_HDF5:BOOL=OFF
" > MPI_DEBUG_STK.config

echo "
-D Trilinos_ENABLE_TriKota:BOOL=OFF
-D Trilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON
-D TPL_ENABLE_Boost:BOOL=ON
-D Boost_INCLUDE_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/code/TPLs_src/boost
-D Boost_LIBRARY_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/build/SierraBuild.trilinos_dev.devsierra/boost/1.46.1/gcc-4.4.4/release/address-model-64/mpi-openmpi/runtime-link-shared
-D TPL_ENABLE_gtest:BOOL=ON
-D gtest_INCLUDE_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/code/TPLs_src/gtest/include
-D gtest_LIBRARY_DIRS:PATH=/scratch/stana/PROJECTS/NightlyTesting/SierraTridev/build/SierraBuild.trilinos_dev.devsierra/gtest/1.5.0/gcc-4.4.4/release/address-model-64/mpi-openmpi/runtime-link-shared/
-D TPL_ENABLE_Netcdf:BOOL=ON
-D TPL_Netcdf_INCLUDE_DIRS:PATH=/projects/seacas/linux_rhel5/current/include
-D Netcdf_LIBRARY_DIRS:PATH=/projects/seacas/linux_rhel5/current/lib
-D SEACAS_ENABLE_XDMF:BOOL=OFF
-D Trilinos_ENABLE_SEACASMat2exo:BOOL=OFF
-D TPL_ENABLE_XMDF:BOOL=OFF
-D TPL_ENABLE_HDF5:BOOL=OFF
" > SERIAL_RELEASE_STK.config


#
# Run the standard checkin testing script with my specializations
#

../Trilinos/checkin-test.py \
--no-eg-git-version-check \
--ss-extra-builds=MPI_DEBUG_STK,SERIAL_RELEASE_STK \
--make-options=-j8 \
--ctest-options=-j4 \
--ctest-timeout=180 \
$EXTRA_ARGS  

