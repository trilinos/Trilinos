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

echo > MPI_DEBUG.config
echo "MPI_BASE_DIR:PATH=/home/sntools/extras/mpi/mpich-1.2.7p1-gcc-4.2.4-64Bit" >> MPI_DEBUG.config
echo "MPI_EXEC_PRE_NUMPROCS_FLAGS:STRING=--all-local" >> MPI_DEBUG.config

#
# Run the standard checkin testing script with my specializations
#

/var/scratch/rabartl/PROJECTS/Sierra.new/Aria_Trilinos/code/TPLs_src/Trilinos/dev/cmake/python/checkin-test.py \
--make-options="-j8" \
--ctest-options="-j8" \
--without-serial-release \
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS  
