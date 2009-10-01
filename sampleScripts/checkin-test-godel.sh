#!/bin/bash

#
# This is the script that I used to checkin to Trilinos on godel.sandia.gov.
# You can copy this script and adapt it to your own machine.
#

EXTRA_ARGS=$@

#
# Set up basic environment options
#

echo "-DBUILD_SHARED_LIBS:BOOL=ON" > COMMON.config

echo "-DMPI_BASE_DIR:PATH=/usr/lib64/openmpi/1.2.7-gcc" > MPI_DEBUG.config

echo "-DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++"     >  SERIAL_RELEASE.config
echo "-DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc"       >> SERIAL_RELEASE.config
echo "-DCMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/f77" >> SERIAL_RELEASE.config


#
# Run the checkin-test.py script with more arguments
#


/home/rabartl/PROJECTS/Trilinos.base/Trilinos/cmake/python/checkin-test.py \
--make-options="-j6" \
--ctest-options="-j6" \
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS  
