#!/bin/bash

#
# This is the script that I use for remote test/push on
# godel.sandia.gov pulling from gabriel.sandia.gov.  You can copy this
# script and adapt it to your own machine.
#
# To do a remote test/push by pulling from gabriel.sandia.gov, just pass in the options:
#
#   --do-all --push
#
# then just wait for the return email
#

EXTRA_ARGS=$@

#
# Set up basic environment options
#

echo "-DBUILD_SHARED_LIBS:BOOL=ON" > COMMON.config

echo "-DMPI_BASE_DIR:PATH=/usr/lib64/openmpi/1.2.7-gcc" > MPI_DEBUG.config

echo "
-DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++
-DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc
-DCMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/f77
" > SERIAL_RELEASE.config


#
# Run the checkin-test.py script with more arguments
#

../../Trilinos/checkin-test.py \
--make-options="-j6" \
--ctest-options="-j6" \
--ctest-timeout=180 \
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS

# NOTE: Even though godel has 8 cores, I will only use 6 of them so that I
# don't dominate the machine.

# Above, I left in the --commit-msg-header option in case I need to
# fix problems and add more commits.  I like editing a file before I
# commit.
