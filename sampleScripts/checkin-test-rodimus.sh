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
-DTPL_BLAS_LIBRARIES:STRING=/usr/lib64/libblas.so.3
-DTPL_LAPACK_LIBRARIES:STRING=/usr/lib64/liblapack.so.3
" > COMMON.config

echo "
-DMPI_BASE_DIR:PATH=/home/dridzal/software/openmpi-1.4.3/install
" > MPI_DEBUG.config

echo "
-DCMAKE_CXX_COMPILER:FILEPATH=/home/dridzal/software/gcc/install_gcc-4.5.1/bin/g++
-DCMAKE_C_COMPILER:FILEPATH=/home/dridzal/software/gcc/install_gcc-4.5.1/bin/gcc
-DCMAKE_Fortran_COMPILER:FILEPATH=/home/dridzal/software/gcc/install_gcc-4.5.1/bin/gfortran
" > SERIAL_RELEASE.config

#
# Run the standard checkin testing script with my specializations
#

export LD_LIBRARY_PATH=/home/dridzal/software/openmpi-1.4.3/gcc-4.5.1/lib:/home/dridzal/software/gcc/install_gcc-4.5.1/lib64/:/home/dridzal/software/gcc/install_mpc/lib:/home/dridzal/software/gcc/install_mpfr/lib:/home/dridzal/software/gcc/install_gmp/lib:

../../checkin-test.py \
--make-options=-j12 \
--ctest-options=-j12 \
--ctest-timeout=180 \
$EXTRA_ARGS  
