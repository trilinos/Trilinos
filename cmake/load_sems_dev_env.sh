# Load the standard Trilinos SEMS develoment environment
#
# This script is sourced as:
#
#   $ source load_sems_dev_env.sh \
#     [ <compiler-and-version>] [ <mpi-and-version> [<cmake-and-version>] ] ]
#
# where:
#
#   <compiler-and-version> is the SEMS module name for the compiler and
#   version, e.g. sems-gcc/4.8.4 (default).
#
#   <mpi-and-version> is the module name for the MPI implementation and
#   version, e.g. sems-openmpi/1.8.7 (default).
#
# Once sourced, this module also laods the SEMS modules for all of the TPLs
# and tools that Trilinos can use (see below for details).
#
# One can use the default compiler/version and OpenMPI/version with:
#
#   $ source load_sems_dev_env.sh
#
# or one can specify a different compiler/version but use the default
# OpenMPI/version using, for example:
#
#   $ source load_sems_dev_env.sh sems-clang/3.5.2
#
# (which uses the default OpenMPI version, see below).
#
# or one can specify a different compiler/version and mpi/version with:
#
#   $ source load_sems_dev_env.sh sems-clang/3.5.2 sems-openmpi/1.87
#
# or one can specify a different compiler/version, mpi/version, and
# cmake/version with:
#
#   $ source load_sems_dev_env.sh sems-clang/3.5.2 sems-openmpi/1.87 \
#     sems-cmake/3.5.2
#
# For verbose output, set the env var TRILINOS_SEMS_DEV_ENV_VERBOSE=1.
#
# To see the modules that are loaded, run:
#
#   $ module list
#
# To unload a loaded Trilinos SEMS Dev Env, do:
#
#   $ source unload_sems_dev_env.sh
#
# If a current Trilinos SEMS Dev Env is loaded, this script will automatically
# unload the curent version and load the new reqeusted version.
#
# NOTE: After sourcing this script one can safely unload the current modules
# and load different modules for Python and CMake since these are just tools
# that don't depend on the selected compiler or MPI.  Also, one can also
# unload the default Boost module and load a different version since Boost
# does not depend on the MPI implementation and no other loaded TPLs depend on
# Boost.
#
# NOTE: If the modules don't load correctly because of some conflict , then
# consider calling:
#
#   $ module purge
#
# to clean house and then source this script again.
#

# Get the base dir for the sourced script
called=$_
#[[ $called != $0 ]] && echo "Script is being sourced" || echo "Script is being run"
#echo "\$BASH_SOURCE ${BASH_SOURCE[@]}"
_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "_SCRIPT_DIR = '$_SCRIPT_DIR'"

#
# A) Get the SEMS Dev Env to load
#

sems_compiler_and_version_load=$1
if [ "$sems_compiler_and_version_load" == "" ] ; then
  sems_compiler_and_version_load=sems-gcc/4.8.4
fi
#echo "sems_compiler_and_version_load = $sems_compiler_and_version_load"

sems_mpi_and_version_load=$2
if [ "$sems_mpi_and_version_load" == "" ] ; then
  sems_mpi_and_version_load=sems-openmpi/1.8.7
fi
#echo "sems_mpi_and_version_load = $sems_mpi_and_version_load"

sems_cmake_and_version_load=$3
if [ "$sems_cmake_and_version_load" == "" ] ; then
  sems_cmake_and_version_load=sems-cmake/3.5.2
fi
#echo "sems_cmake_and_version_load = $sems_cmake_and_version_load"

TRILINOS_SEMS_DEV_ENV_TO_LOAD="$sems_compiler_and_version_load $sems_mpi_and_version_load $sems_cmake_and_version_load"

#
# B) If a SEMS Dev Env is already loaded, then unload the current one and load
# the new on.
#

if [ "$TRILINOS_SEMS_DEV_ENV_LOADED" != "" ] ; then
  if [ "$TRILINOS_SEMS_DEV_ENV_LOADED" == "$TRILINOS_SEMS_DEV_ENV_TO_LOAD" ] ; then
    if [ "${TRILINOS_SEMS_DEV_ENV_VERBOSE}" == "1" ] ; then
      echo "Already loaded the right Trilinos SEMS Dev Env" \
        "'$TRILINOS_SEMS_DEV_ENV_TO_LOAD' so just return!"
    fi
    return 0
  fi
  if [ "$TRILINOS_SEMS_DEV_ENV_LOADED" != "$TRILINOS_SEMS_DEV_ENV_TO_LOAD" ] ; then
    . $_SCRIPT_DIR/unload_sems_dev_env.sh
  fi
fi

if [ "${TRILINOS_SEMS_DEV_ENV_VERBOSE}" == "1" ] ; then
  echo "Loading Trilinos SEMS Dev Env =" \
    "'$TRILINOS_SEMS_DEV_ENV_TO_LOAD'!"
fi

#
# C) Load the modules
#

module load sems-env # In case user did 'module purge'
module load sems-python/2.7.9 
module load $sems_cmake_and_version_load
module load $sems_compiler_and_version_load
module load $sems_mpi_and_version_load
module load sems-boost/1.55.0/base 
module load sems-zlib/1.2.8/base 
module load sems-hdf5/1.8.12/parallel 
module load sems-netcdf/4.3.2/parallel 
module load sems-parmetis/4.0.3/parallel 
module load sems-scotch/6.0.3/parallel 

if [ "${TRILINOS_SEMS_DEV_ENV_VERBOSE}" == "1" ] ; then
  module list
fi

#
# D) Remember the loaded SEMS Dev Env
#

export TRILINOS_SEMS_DEV_ENV_LOADED="$TRILINOS_SEMS_DEV_ENV_TO_LOAD"
