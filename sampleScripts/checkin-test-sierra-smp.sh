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
#                   checkin-test-sierra-smp.sh
#               MPI_DEBUG/
#               SERIAL_DEBUG/
#
# Make sure to "cd ~/Codes/Trilinos/build/CHECKIN; ln -s ../../sampleScripts/checkin-test-sierra-smp.sh ."
#
# To use this script just run it, for example, like (this builds and tests):
#
#  $ ./checkin-test-sierra-smp.sh --do-all
#
# (This pushes after successful tests have been run)
#
#  $ ./checkin-test-sierra-smp.sh --push

for word in "$@"
do
  if [ ${word} == "--help" ]; then
    echo "
To run this script (this builds and tests):

$ ./checkin-test-sierra-smp.sh --do-all

(This pushes after successful tests have been run)

$ ./checkin-test-sierra-smp.sh --push "
exit

  fi
done



#
# Set up configuration files
#
BOOSTDIR=/ascldap/users/gahanse/local/boost_1_50_0
MPI_BASE_DIR=/sierra/Sntools/extras/mpi/RHEL5/openmpi-1.4.5/gcc-4.6.3-64Bit
NETCDF=/ascldap/users/gahanse/local/parallel
HDFDIR=/ascldap/users/gahanse/local/parallel
BLASDIR=/sierra/Sntools/extras/compilers/acml4.3.0/gfortran64/lib
FORTRAN_LIBS="-lgfortran -lboost_system"

echo "
-D BUILD_SHARED_LIBS:BOOL=ON
-D Boost_INCLUDE_DIRS:FILEPATH=\"$BOOSTDIR/include\"
-D Boost_LIBRARY_DIRS:FILEPATH=\"$BOOSTDIR/lib\"
-D TPL_ENABLE_Boost:BOOL=ON
-D TPL_ENABLE_Netcdf:BOOL=ON
-D TPL_Netcdf_INCLUDE_DIRS:PATH=$NETCDF/include
-D TPL_Netcdf_LIBRARY_DIRS:PATH=$NETCDF/lib
-D Trilinos_EXTRA_LINK_FLAGS:STRING=\"-lgfortran -lboost_system\"
-D BLAS_LIBRARY_DIRS:FILEPATH=$BLASDIR \
-D BLAS_LIBRARY_NAMES:STRING=\"acml\"
-D LAPACK_LIBRARY_DIRS:FILEPATH=$BLASDIR
-D LAPACK_LIBRARY_NAMES:STRING=\"acml\"
" > COMMON.config

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
-DTPL_ENABLE_MPI:BOOL=ON 
-DMPI_BASE_DIR:PATH=$MPI_BASE_DIR
" > MPI_DEBUG.config

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
-DTPL_ENABLE_MPI:BOOL=ON 
-DMPI_BASE_DIR:PATH=$MPI_BASE_DIR
" > MPI_DEBUG_SS.config

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
" > SERIAL_DEBUG.config

EXTRA_ARGS=$@

../../checkin-test.py \
--make-options="-j8" \
--ctest-options="-j8" \
--send-email-to="" \
--no-eg-git-version-check \
$EXTRA_ARGS
