#!/bin/bash

# Simple script I used to drive checkin-test.py from within a CHECKIN subdir of the Trilinos
# source repo, which uses a CHECKIN subdir to store the output files.  This
# makes it easy to do commands like --pull with extra repos to see what the
# status of things are. 
#
# Here is my directory structure:
#   Codes/
#       Trilinos/
#           packages/
#           checkin-test.py
#           checkin_message
#       TriBuild/
#           CHECKIN/
#               checkin-test-neams-smp.sh
#           MPI_DEBUG/
#           SERIAL_DEBUG/
#
# Make sure to "cd ~/Codes/TriBuild/CHECKIN; ln -s ../../Trilinos/sampleScripts/checkin-test-neams-smp.sh ."
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
BOOSTDIR=/ascldap/users/gahanse/local/boost_1_47_0
MPI_BASE_DIR=/usr/netpub/mpi/OpenMPI/1.4/64Bit/gnu
NETCDF=/ascldap/users/gahanse/local

echo "
-D BUILD_SHARED_LIBS:BOOL=ON
-D Boost_INCLUDE_DIRS:FILEPATH=\"$BOOSTDIR\"
-D Boost_LIBRARY_DIRS:FILEPATH=\"$BOOSTDIR/stage/lib\"
-D TPL_ENABLE_Boost:BOOL=ON
-D TPL_ENABLE_Netcdf:BOOL=ON
-D TPL_Netcdf_INCLUDE_DIRS:PATH=$NETCDF/include
-D TPL_Netcdf_LIBRARY_DIRS:PATH=$NETCDF/lib
" > COMMON.config

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
-DTPL_ENABLE_MPI:BOOL=ON 
-DMPI_BASE_DIR:PATH=$MPI_BASE_DIR
" > MPI_DEBUG.config

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
" > SERIAL_DEBUG.config

EXTRA_ARGS=$@

../../Trilinos/checkin-test.py \
--make-options="-j8" \
--ctest-options="-j8" \
--send-email-to="" \
--no-eg-git-version-check \
$EXTRA_ARGS
