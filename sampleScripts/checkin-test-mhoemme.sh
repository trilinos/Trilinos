#!/bin/bash
#
# Mark Hoemmen <mhoemme@sandia.gov>
# Check-in test script for my Sandia workstation sponge-bob.sandia.gov
#   (Linux (RHEL), x86_64, 8 cores).
#
# This script was based originally on checkin-test-gabriel.sh, and
# then heavily modified.  You are welcome to copy and modify this
# script for your own purposes, subject to the same software license
# under which Trilinos is released.
#
# Unusual features: 
#
# * Custom GCC and OpenMPI user-space installs.  PATH and
#   LD_LIBRARY_PATH are locally modified to ensure that the right C,
#   C++, Fortran, and MPI wrappers are invoked, both for the build and
#   for testing (the MPI_EXEC variable in particular).
#
# * Shared libraries are enabled; this accelerates the build.
#
# * AMD's vendor BLAS and LAPACK implementation (ACML) is used, rather
#   than the default system install.  This is a nonthreaded (but
#   thread-safe) version.  If you use a threaded version, you might
#   need to set the number of threads to 1 via an environment
#   variable, so that it doesn't contend for parallel resources when
#   invoked in parallel from multiple threads.  
#
# * I build with Boost enabled, using a user-space install of Boost.
#   You'll need the Boost libraries installed, not just the header
#   files.
#
# Notes:
#
# On occasion, it's useful to exclude a particular test from the
# check-in tests.  For example, the test might work everywhere else
# but on your computer (I've had this happen once).  You can do this
# by passing in the name of the test as one of the ctest options; you
# just have to be careful how you quote the option, otherwise the name
# of the test to exclude goes away.  For example, in order to exclude
# the following test:
#
# Teuchos_Dependencies_Serialization_test
#
# run this script as follows:
#
# ./checkin-test-mhoemme.sh --do-all --push --ctest-options=\"-E Teuchos_Dependencies_Serialization_test\"
#
# You need the backslashes before the double quotes; otherwise,
# ctest's options get set to "-E".
#
# You should specify the number of processes that cmake and ctest
# should use via a "-jN" argument to this script: e.g., -j9 to use 9
# processes.  It's typical practice to oversubscribe CPUs slightly for
# parallel builds, in order to overlap disk I/O and computation.


# If the first argument isn't already in the PATH, prepend it to the
# PATH (or append it if the second argument is "after").
concat_path () {
    if ! echo $PATH | /bin/egrep -q "(^|:)$1($|:)" ; then
	if [ "$2" = "after" ] ; then
	    PATH=$PATH:$1
	else
	    PATH=$1:$PATH
	fi
    fi
    export PATH
}

# If the first argument isn't already in LD_LIBRARY_PATH, prepend it
# to LD_LIBRARY_PATH (or append it if the second argument is "after").
concat_ld_library_path () {
    if ! echo $LD_LIBRARY_PATH | /bin/egrep -q "(^|:)$1($|:)" ; then
	if [ "$2" = "after" ] ; then
	    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$1
	else
	    LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
	fi
    fi
    export LD_LIBRARY_PATH
}

# Make sure that the build happens with the user-space install of GCC.
#
# The exports of PATH and LD_LIBRARY_PATH won't escape this script.
if [ -d $HOME/pkg/gcc-4.5.1/lib64 ]; then
    concat_ld_library_path $HOME/pkg/gcc-4.5.1/lib64
else
    concat_ld_library_path $HOME/pkg/gcc-4.5.1/lib
fi
concat_path $HOME/pkg/gcc-4.5.1/bin

# Make sure that the build happens with the correct MPI implementation.
concat_path $HOME/pkg/openmpi-1.4.3/bin
concat_ld_library_path $HOME/pkg/openmpi-1.4.3/lib

# Print out which GCC version will be used for the build.
echo "GCC and G++ version used for the build:"
echo ""
echo `which gcc`
echo `which g++`
echo ""

# You can pass any valid arguments for the original check-in script as
# extra arguments to this wrapper script.  We pick up the extra
# arguments here.
EXTRA_ARGS=$@

#
# Configuration options for all builds
#
echo "
-D BUILD_SHARED_LIBS:BOOL=ON
-D BLAS_LIBRARY_DIRS:FILEPATH=\"/home/mhoemme/pkg/acml4.4.0/gfortran64/lib\"
-D BLAS_LIBRARY_NAMES:STRING=\"acml\"
-D LAPACK_LIBRARY_DIRS:FILEPATH=\"/home/mhoemme/pkg/acml4.4.0/gfortran64/lib\"
-D LAPACK_LIBRARY_NAMES:STRING=\"acml\"
-D TPL_ENABLE_TBB:BOOL=OFF
-D TPL_ENABLE_Boost:BOOL=ON
-D Boost_LIBRARY_DIRS:FILEPATH=\"${HOME}/pkg/boost/lib\"
-D Boost_INCLUDE_DIRS:FILEPATH=\"${HOME}/pkg/boost/include\"
-D Kokkos_ENABLE_TSQR:BOOL=ON
-D Belos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON
" > COMMON.config

#
# Configuration options for MPI debug build
#
echo "
-D CMAKE_BUILD_TYPE:STRING=RELEASE
-D MPI_BASE_DIR:FILEPATH=\"${HOME}/pkg/openmpi-1.4.3\"
-D MPI_EXEC:FILEPATH=\"${HOME}/pkg/openmpi-1.4.3/bin/mpiexec\"
-D MPI_Fortran_COMPILER:FILEPATH=\"${HOME}/pkg/openmpi-1.4.3/bin/mpif90\"
-D MPI_CXX_COMPILER:FILEPATH=\"${HOME}/pkg/openmpi-1.4.3/bin/mpicxx\"
-D MPI_C_COMPILER:FILEPATH=\"${HOME}/pkg/openmpi-1.4.3/bin/mpicc\"
" > MPI_DEBUG.config

#
# Run the standard checkin testing script with my specializations
#
../../Trilinos/checkin-test.py \
--ctest-timeout=500 \
--no-eg-git-version-check \
$EXTRA_ARGS  

