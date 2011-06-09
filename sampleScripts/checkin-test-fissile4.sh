#!/bin/bash

# Used to test Trilinos on any of the ORNL fissle 4 machines
# (e.g. u233, u235, pu239, and pu241).

# NOTE: To use this, you must first prepend /opt/trilinos-toolset/bin
# to your path to find eg and cmake!

# NOTE: This script automatically picks up any CASL VRI related extra
# repos and adds them to --extra-repos.  If you want to override that,
# you can just pass in --extra-repos=??? to drop off extra repos or
# select the set that you want.

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

# Must be in the correct order (see ExtraExternalRepositories.casl_vri.cmake)
EXTRA_REPOS_FULL_LIST="StartCCMExt DeCARTExt CASLBOA CASLRAVE LIMEExt PSSDriversExt"
DRIVERS_BASE_DIR="$TRILINOS_BASE_DIR_ABS/Trilinos/cmake/ctest/drivers/pu241"

# Pakage in Trilinos to disable (mostly for auotmated CI server)
DISABLE_PACKAGES=PyTrilinos,TriKota,Optika,Sundance,Stokhos,STK

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=$DRIVERS_BASE_DIR/gcc-4.5.1-mpi-options.cmake
" > MPI_DEBUG.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=$DRIVERS_BASE_DIR/gcc-4.5.1-serial-options.cmake
" > SERIAL_RELEASE.config

#
# Extra intel builds added with --extra-builds=INTEL_RELEASE,...
#

# note: the pvm dirs below can be removed when configure_options_files supports multiple fragments
echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=$DRIVERS_BASE_DIR/intel-12.0.4-options.cmake;$DRIVERS_BASE_DIR/SubmitToCaslDev.cmake;$DRIVERS_BASE_DIR/casl-vri-packages-coupled.cmake
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF
-DDART_TESTING_TIMEOUT:STRING=660
" > VERA_INTEL.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=$DRIVERS_BASE_DIR/intel-12.0.4-options.cmake;$DRIVERS_BASE_DIR/SubmitToCaslDev.cmake;$DRIVERS_BASE_DIR/casl-vri-packages-coupled.cmake
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF
-DDART_TESTING_TIMEOUT:STRING=660
" > VERA_INTEL_VERACOUPLINGS.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=$DRIVERS_BASE_DIR/intel-12.0.4-options.cmake;$DRIVERS_BASE_DIR/SubmitToCaslDev.cmake;$DRIVERS_BASE_DIR/casl-vri-packages-native.cmake
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF
-DTrilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=OFF
-DTrilinos_ALL_PACKAGES:BOOL=OFF
-DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=OFF
-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
-DTrilinos_ENABLE_TESTS:BOOL=ON
-DTrilinos_ENABLE_CASLRAVE:BOOL=ON
-DTrilinos_ENABLE_CASLBOA:BOOL=ON
-DDART_TESTING_TIMEOUT:STRING=660
" > VERA_INTEL_NATIVEBUILDS.config

#
# Load up the list of extra repos based on what is present:
#

EXTRA_REPOS=
for extra_repo in $EXTRA_REPOS_FULL_LIST; do
  #echo $extra_repo
  EXTRA_REPO_PATH=$TRILINOS_BASE_DIR/Trilinos/$extra_repo
  #echo $EXTRA_REPO_PATH
  if [ -d $EXTRA_REPO_PATH ]; then
    EXTRA_REPOS=$EXTRA_REPOS$extra_repo,
  fi
done
#echo "EXTRA_REPOS=$EXTRA_REPOS"

#
# Invocation
#

$TRILINOS_BASE_DIR/Trilinos/checkin-test.py \
--extra-repos=$EXTRA_REPOS \
-j16 \
--ctest-timeout=180 \
--disable-packages=$DISABLE_PACKAGES \
$EXTRA_ARGS  

# NOTE: By default we use 16 processes which is 1/2 of the 32
# processes on this machine.  This way two people can build and test
# Trilinos without taxing the machine too much.
