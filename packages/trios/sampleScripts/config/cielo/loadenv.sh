#!/bin/bash

# This script gets sourced by do-configure.  But it can also be soruced in your shell like this:
# bash> TARGET_COMPILER=pgi source loadenv.sh
#
# The steps taken below can be summarized as:
# - check that TARGET_COMPILER is set.
# - check that TARGET_COMPILER is set to a known value.
# - cleanup the module list by unloading all known PrgEnv modules.
# - load the correct PrgEnv for the TARGET_COMPILER
# - load compiler neutral modules
# - perform any last minute compiler specific module swaps


# Load the appropriate bash environment
. /opt/modules/default/init/bash


if [ -z $TARGET_COMPILER ]; then
  echo "=============================="
  echo "=============================="
  echo "TARGET_COMPILER is not defined in your environment.  Please use one of the known values - pgi, gnu, cray or intel.  Aborting."
  echo "=============================="
  echo "=============================="
  return
fi

case ${TARGET_COMPILER} in
  pgi|gnu|cray|intel)
    # TARGET_COMPILER is recognized.  Carry on.
    ;;
  *)
    echo "=============================="
    echo "=============================="
    echo "'$TARGET_COMPILER' is not a recognized value for TARGET_COMPILER.  Please use one of the known values - pgi, gnu, cray or intel.  Aborting."
    echo "=============================="
    echo "=============================="
    return
    ;;
esac

# Unload any existing PrgEnv modules
module unload PrgEnv-pgi
module unload PrgEnv-gnu
module unload PrgEnv-cray
module unload PrgEnv-intel

module load PrgEnv-${TARGET_COMPILER}


if [ -z $VALGRIND_BUILD ]; then
  # These specific module versions are required (defaults don't work)
  echo "Loading fftw"
  module load fftw/3.2.2.1
  module load xt-papi
else
  echo "Unloading libsci and fftw"
  # Unload any problematic modules
  module unload xt-libsci
  module unload fftw
  module unload xt-papi
fi

# If you need to something special for your compiler.  Do it here.
case ${TARGET_COMPILER} in
  pgi)
    module swap pgi/10.9.0
    ;;
  gnu)
    ;;
  cray)
    ;;
  intel)
    ;;
esac


module list
