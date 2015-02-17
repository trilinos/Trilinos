#!/bin/bash

# Used to test Trilinos on any of the ORNL CASL Fissile 4 machines
#
# This script requires that the VERA dev env be loaded by sourcing the script:
#
#  . /projects/vera/load_dev_env.[sh,csh]
#
# You can source this script either in your shell startup script
# (e.g. .bash_profile) or you can source it manually whenever you need to set
# up to build VERA software.
#
# NOTE: This script should not be directly modifed by typical CASL
# developers except, perhaps to add new extra builds.

EXTRA_ARGS=$@

# The default location for this directory tree is:
#
#  Trilinos.base
#    Trilinos    (your Trilinos soruce tree)
#    BUILDS
#      CHECKIN   (where you run this script from)
#
if [ "$TRILINOS_BASE_DIR" == "" ] ; then
  TRILINOS_BASE_DIR=../..
fi

TRILINOS_BASE_DIR_ABS=$(readlink -f $TRILINOS_BASE_DIR)

DRIVERS_BASE_DIR="$TRILINOS_BASE_DIR_ABS/Trilinos/cmake/ctest/drivers/pu241"

# Packages in Trilinos to disable (mostly for auotmated CI server)
DISABLE_PACKAGES=CTrilinos,ForTrilinos,PyTrilinos,Didasko,Mesquite,Phdmesh,Pliris,Claps,Amesos2,STK,FEApp,TriKota,Optika

# Check to make sure that the env has been loaded correctly
if [ "$LOADED_VERA_DEV_ENV" != "gcc461" ] ; then
  echo "Error, must source /projects/vera/gcc-4.6.1/load_dev_env.[sh,csh] before running checkin-test-vera.sh!"
  exit 1
fi

echo "
-DTrilinos_ENABLE_CXX11=OFF
-DTrilinos_EXCLUDE_PACKAGES=CTrilinos
-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
" > COMMON.config

#
# Built-in Primary Stable (PS) builds (DO NOT MODIFY)
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.6.1-mpi-debug-ps-options.cmake'
" > MPI_DEBUG.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.6.1-serial-release-ps-options.cmake'
" > SERIAL_RELEASE.config

#
# Standard Secondary Stable (SS) builds (DO NOT MODIFY)
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.6.1-mpi-debug-ss-options.cmake,$DRIVERS_BASE_DIR/trilinos-tpls-gcc.4.6.1.cmake'
" > MPI_DEBUG_ST.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.6.1-serial-release-ss-options.cmake,$DRIVERS_BASE_DIR/trilinos-tpls-gcc.4.6.1.cmake'
" > SERIAL_RELEASE_ST.config

#
# Extra builds
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.6.1-mpi-debug-ps-options.cmake'
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=OFF
" > MPI_RELEASE.config

#
# Invocation
#

$TRILINOS_BASE_DIR/Trilinos/checkin-test.py \
-j16 \
--ctest-timeout=180 \
--st-extra-builds=MPI_DEBUG_ST,SERIAL_RELEASE_ST \
--disable-packages=$DISABLE_PACKAGES \
--skip-case-no-email \
--ctest-options="-E '(Piro_AnalysisDriver|Stokhos_Linear2D_Diffusion_GMRES_KLR|Panzer_STK_ResponseLibraryTest|MueLu_|Amesos2_|Rythmos_ImplicitRK_UnitTest_MPI_1|SEACASExodus_exodus_unit_tests|Intrepid_test_Discretization_Basis_HGRAD_TRI_Cn_FEM_Test_02_MPI_1|Intrepid_test_Discretization_Basis_HDIV_TET_In_FEM_Test_02_MPI_1|Intrepid_test_Discretization_Basis_HGRAD_TET_Cn_FEM_Test_02_MPI_1|Sundance_BesselTest2D_MPI_1|ThyraTpetraAdapters_TpetraThyraWrappersUnitTests_serial|Ifpack2_RILUKSingleProcessUnitTests)'" \
$EXTRA_ARGS

# NOTE: By default we use 16 processes which is 1/2 of the 32 processes on a
# fissile 4 machine.  This way two people can build and test without taxing
# the machine too much.

# Here is the template for excluding tests that can be aded to
# checkin-test.py call above
#
#   --ctest-options="-E '(Test1|Test2)'" \
