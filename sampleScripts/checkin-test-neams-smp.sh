#!/bin/bash

# Simple script I used to drive checkin-test.py from within a CHECKIN subdir of the Trilinos
# source repo, which uses a CHECKIN subdir to store the output files.  This
# makes it easy to do commands like --pull with extra repos to see what the
# status of things are. 
#
# Glen Hansen, gahanse@sandia.gov
#
# Here is my directory structure:
#   Codes/
#       Trilinos/
#           packages/
#           checkin-test.py
#           checkin_message
#           build/
#               CHECKIN/
#                   checkin-test-neams-smp.sh
#               MPI_DEBUG/
#               SERIAL_DEBUG/
#
# Make sure to "cd ~/Codes/Trilinos/build/CHECKIN; ln -s ../../sampleScripts/checkin-test-neams-smp.sh ."
#
# To use this script just run it, for example, like (this builds and tests):
#
#  $ ./checkin-test-neams-smp.sh --do-all
#
# (This pushes after successful tests have been run)
#
#  $ ./checkin-test-neams-smp.sh --push

for word in "$@"
do
  if [ ${word} == "--help" ]; then
    echo "
To run this script (this builds and tests):

$ ./checkin-test-neams-smp.sh --do-all

(This pushes after successful tests have been run)

$ ./checkin-test-neams-smp.sh --push "
exit

  fi
done



#
# Set up configuration files
#
BASE=/repository/usr/local
BOOSTDIR=$BASE/boost/boost_1_51_0
MPI_BASE_DIR=$BASE/gcc/gcc-4.7.2
NETCDF=$BASE/parallel/netcdf-4.2.1.1
HDFDIR=$BASE/parallel/hdf5-1.8.9
SNETCDF=$BASE/serial/netcdf-4.2.1.1
SHDFDIR=$BASE/serial/hdf5-1.8.9
MKL_PATH=$BASE/intel
LABLAS_LIBRARY_DIRS="$MKL_PATH/mkl/lib/intel64;/usr/lib64"
LABLAS_LIBRARY_NAMES="mkl_intel_lp64;mkl_sequential;mkl_core;pthread"

echo "
-D BUILD_SHARED_LIBS:BOOL=ON
-D Boost_INCLUDE_DIRS:FILEPATH=\"$BOOSTDIR/include\"
-D Boost_LIBRARY_DIRS:FILEPATH=\"$BOOSTDIR/lib\"
-D TPL_ENABLE_Boost:BOOL=ON
-D TPL_ENABLE_Netcdf:BOOL=ON
-D TPL_ENABLE_HDF5:BOOL=ON
-D TPL_ENABLE_BLAS:STRING=ON 
-D TPL_ENABLE_LAPACK:STRING=ON 
-D BLAS_LIBRARY_DIRS:STRING=$LABLAS_LIBRARY_DIRS 
-D BLAS_LIBRARY_NAMES:STRING=$LABLAS_LIBRARY_NAMES 
-D LAPACK_LIBRARY_DIRS:STRING=$LABLAS_LIBRARY_DIRS 
-D LAPACK_LIBRARY_NAMES:STRING=$LABLAS_LIBRARY_NAMES 
" > COMMON.config

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
-DTPL_ENABLE_MPI:BOOL=ON 
-DMPI_BASE_DIR:PATH=$MPI_BASE_DIR
-D Netcdf_INCLUDE_DIRS:PATH=$NETCDF/include 
-D Netcdf_LIBRARY_DIRS:PATH=$NETCDF/lib 
-D HDF5_INCLUDE_DIRS:PATH=$HDFDIR/include 
-D HDF5_LIBRARY_DIRS:PATH=$HDFDIR/lib 
" > MPI_DEBUG.config

#echo "
#-DCMAKE_BUILD_TYPE:STRING=DEBUG
#-DTPL_ENABLE_MPI:BOOL=ON 
#-DMPI_BASE_DIR:PATH=$MPI_BASE_DIR
#-D Netcdf_INCLUDE_DIRS:PATH=$NETCDF/include 
#-D Netcdf_LIBRARY_DIRS:PATH=$NETCDF/lib 
#-D HDF5_INCLUDE_DIRS:PATH=$HDFDIR/include 
#-D HDF5_LIBRARY_DIRS:PATH=$HDFDIR/lib 
#" > MPI_DEBUG_SS.config

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
-D Netcdf_INCLUDE_DIRS:PATH=$SNETCDF/include 
-D Netcdf_LIBRARY_DIRS:PATH=$SNETCDF/lib 
-D HDF5_INCLUDE_DIRS:PATH=$SHDFDIR/include 
-D HDF5_LIBRARY_DIRS:PATH=$SHDFDIR/lib 
" > SERIAL_DEBUG.config

EXTRA_ARGS=$@

../../checkin-test.py \
--make-options="-j 16" \
--ctest-options="-j 16" \
--send-email-to="" \
--no-eg-git-version-check \
$EXTRA_ARGS
