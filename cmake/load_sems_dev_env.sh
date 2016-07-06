#
# Load the standard Trilinos SEMS develoment environment
#
# USAGE:
#
#   $ source load_sems_dev_env.sh \
#       [<compiler-and-version>] [<mpi-and-version>] [<cmake-and-version>]
#
# where:
#
#   <compiler-and-version> is the SEMS module name for the compiler and
#   version, e.g. sems-gcc/4.8.4.  To use the default just pass in "default"
#   (see the default listed below).
#
#   <mpi-and-version> is the SEMS module name for the MPI implementation and
#   version, e.g. sems-openmpi/1.8.7.  To use the default just pass in
#   "default" (see the default listed below).
#
#   <cmake-and-version> is the SEMS module name for the CMake and its version,
#   e.g. sems-cmake/3.5.2.  To use the default just pass in "default" (see the
#   default listed below).
#
# Once sourced, this module also laods the SEMS modules for all of the TPLs
# and tools that Trilinos can use (see below for for the list of TPLs that get
# loaded).
#
# One can use the defaults with just:
#
#   $ source load_sems_dev_env.sh
#
# or one can specify a different compiler/version with:
#
#   $ source load_sems_dev_env.sh sems-clang/3.5.2
#
# or one can specify a different mpi/version with:
#
#   $ source load_sems_dev_env.sh default sems-openmpi/1.8.7
#
# or one can specify a different cmake/version with:
#
#   $ source load_sems_dev_env.sh default default sems-cmake/3.3.2
#
# If a different Trilinos SEMS Dev Env is alreadey loaded, this script will
# automatically unload the curent version and load the new reqeusted version.
#
# Set the env var:
#
#   $ export TRILINOS_SEMS_DEV_ENV_VERBOSE=1
#
# before sourcing this script in order to see verbose output (otherwise, it
# produces no output).
#
# To see the modules that are loaded, run:
#
#   $ module list
#
# To unload a loaded Trilinos SEMS Dev Env, do:
#
#   $ source unload_sems_dev_env.sh
#
# NOTE: After sourcing this script one can safely unload the current modules
# and load different modules for Python and CMake since these are just tools
# that don't depend on the selected compiler or MPI.  Also, one can unload the
# default Boost module and load a different version since Boost does not
# depend on the MPI implementation and no other loaded TPLs depend on Boost.
#
# NOTE: If the modules don't load correctly because of some conflict or some
# other issue, then one can clean house by calling:
#
#   $ module purge
#   $ unset TRILINOS_SEMS_DEV_ENV_LOADED
#
# Then one can source this script again in order to load the desired Trilinos
# SEMS Dev Env.
#

# Module defaults
sems_compiler_and_version_default=sems-gcc/4.8.4
sems_mpi_and_version_default=sems-openmpi/1.8.7
sems_cmake_and_version_default=sems-cmake/3.5.2

# Other defaults
sems_python_and_version_default=sems-python/2.7.9
sems_boost_and_version_default=sems-boost/1.55.0/base
# NOTE: If you change these versions, you also have to update the same list in
# the unload_sems_dev_env.sh module as well!

# Get the base dir for the sourced script
called=$_
#[[ $called != $0 ]] && echo "Script is being sourced" || echo "Script is being run"
#echo "\$BASH_SOURCE ${BASH_SOURCE[@]}"
_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "_SCRIPT_DIR = '$_SCRIPT_DIR'"

#
# A) Get the SEMS Dev Env to load and set defaults
#

sems_compiler_and_version_load=$1
if [ "$sems_compiler_and_version_load" == "" ] || [ "$sems_compiler_and_version_load" == "default" ] ; then
  sems_compiler_and_version_load=$sems_compiler_and_version_default
fi
#echo "sems_compiler_and_version_load = $sems_compiler_and_version_load"

sems_mpi_and_version_load=$2
if [ "$sems_mpi_and_version_load" == "" ] || [ "$sems_mpi_and_version_load" == "default" ] ; then
  sems_mpi_and_version_load=$sems_mpi_and_version_default
fi
#echo "sems_mpi_and_version_load = $sems_mpi_and_version_load"

sems_cmake_and_version_load=$3
if [ "$sems_cmake_and_version_load" == "" ] || [ "$sems_cmake_and_version_load" == "default" ] ; then
  sems_cmake_and_version_load=$sems_cmake_and_version_default
fi
#echo "sems_cmake_and_version_load = $sems_cmake_and_version_load"

TRILINOS_SEMS_DEV_ENV_TO_LOAD="$sems_compiler_and_version_load $sems_mpi_and_version_load $sems_cmake_and_version_load"

#
# B) If a SEMS Dev Env is already loaded, then unload the current one and load
# the new one.
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
# C) Load the modules (in the correct order)
#

module load sems-env # In case user did 'module purge'
module load $sems_python_and_version_default
module load $sems_cmake_and_version_load
module load $sems_compiler_and_version_load
module load $sems_mpi_and_version_load
module load $sems_boost_and_version_default
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
