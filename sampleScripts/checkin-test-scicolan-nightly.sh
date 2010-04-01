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
" > COMMON.config

echo "
-DMPI_BASE_DIR:PATH=/home/sntools/extras/mpi/mpich-1.2.7p1-gcc-4.2.4-64Bit
-DMPI_EXEC_PRE_NUMPROCS_FLAGS:STRING=--all-local
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
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS  
