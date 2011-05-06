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
#    Trilinos
#    BUILDS
#      CHECKIN
#
if [ "$TRILINOS_BASE_DIR" == "" ] ; then
  TRILINOS_BASE_DIR=../..
fi

TRILINOS_TOOLSET_BASE=/opt/gcc-4.5.1/trilinos-toolset

EXTRA_REPOS_FULL_LIST="LIMEExt PSSDriversExt"

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

#
# Extra intel builds added with --extra-builds=INTEL_11064_SERIAL_DEBUG,...
#
# NOTE: You must do 'source /opt/casldev/env/casl_dev_env.sh' before
# using the intel builds.

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DBUILD_SHARED_LIBS:BOOL=ON
-DTPL_ENABLE_Boost:BOOL=ON
-DBoost_INCLUDE_DIRS:FILEPATH=/opt/tpls_src/boost_1_46_1
-DTPL_ENABLE_BinUtils:BOOL=ON
-DCMAKE_C_COMPILER:FILEPATH=/opt/intel/Compiler/11.1/064/bin/intel64/icc
-DCMAKE_CXX_COMPILER:FILEPATH=/opt/intel/Compiler/11.1/064/bin/intel64/icpc
-DCMAKE_Fortran_COMPILER:FILEPATH=/opt/intel/Compiler/11.1/064/bin/intel64/ifort
-DTPL_BLAS_LIBRARIES:STRING='-L${MKLROOT}/lib/em64t -lmkl_intel_lp64 -lmkl_blas95_lp64 -lmkl_core -lmkl_sequential'
-DTPL_LAPACK_LIBRARIES:STRING='-L${MKLROOT}/lib/em64t -lmkl_lapack95_lp64'
-DTrilinos_ENABLE_TESTS:BOOL=ON
-DDART_TESTING_TIMEOUT:STRING=180.0
-DTeuchos_ENABLE_STACKTRACE:BOOL=ON
" > INTEL_11064_SERIAL_DEBUG.config


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
$EXTRA_ARGS  

# NOTE: By default we use 16 processes which is 1/2 of the 32
# processes on this machine.  This way two people can build and test
# Trilinos without taxing the machine too much.
