#
# Purge the modules and load the standard Trilinos SEMS development
# environment
#
# USAGE:
#
#   $ source load_sems_dev_env.sh \
#       [<compiler-and-version>] [<mpi-and-version>] [<cmake-and-version>]
#
# where:
#
#   <compiler-and-version> is the SEMS module name for the compiler and
#   version, e.g. 'sems-gcc/4.8.4'.  To use the default just pass in "default"
#   (see the default listed below).
#
#   <mpi-and-version> is the SEMS module name for the MPI implementation and
#   version, e.g. 'sems-openmpi/1.6.5'.  To use the default just pass in
#   "default" (see the default listed below).
#
#   <cmake-and-version> is the SEMS module name for the CMake and its version,
#   e.g. 'sems-cmake/3.5.2'.  To use the default just pass in "default" (see
#   the default listed below).
#
# Once sourced, this script also loads the SEMS modules for all of the TPLs
# and tools that Trilinos can use that are provided by SEMS (see below for the
# list of TPLs that get loaded).
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
# This script will purge the current modules before loading the requested env
# modules.
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
# and load different modules for CMake, git, and Python since these are just
# tools that don't depend on the selected compiler or MPI.
#

# Get the base dir for the sourced script
called=$_
#[[ $called != $0 ]] && echo "Script is being sourced" || echo "Script is being run"
#echo "\$BASH_SOURCE ${BASH_SOURCE[@]}"
_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "_SCRIPT_DIR = '$_SCRIPT_DIR'"

if [ "$4" != "" ] ; then
  echo "ERROR, the source script 'load_sems_dev_env.sh' does not take more than 3 args! (Remove '$4' ...)"
  return 1
fi

source $_SCRIPT_DIR/std/sems/get_default_modules.sh

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
# B) Purge the current set of modules
#

module purge

#
# C) Load the modules (in the correct order)
#

module load sems-env
module load $sems_python_and_version_default
module load $sems_cmake_and_version_load
module load $sems_git_and_version_default

# The SEMS Intel modules point to an unsupported version of GCC.
# until this is fixed, the workaround is below.
# Please see https://github.com/trilinos/Trilinos/issues/2142
# for updates regarding the right solution.
if [[ $sems_compiler_and_version_load == "sems-intel/"* ]]; then
  module load sems-gcc/4.8.4
fi

module load $sems_compiler_and_version_load
module load $sems_mpi_and_version_load
module load $sems_boost_and_version_default
module load $sems_zlib_and_version_default
module load $sems_hdf5_and_version_default
module load $sems_netcdf_and_version_default
module load $sems_parmetis_and_version_default
module load $sems_scotch_and_version_default
module load $sems_superlu_and_version_default

if [ "${TRILINOS_SEMS_DEV_ENV_VERBOSE}" == "1" ] ; then
  module list
fi
