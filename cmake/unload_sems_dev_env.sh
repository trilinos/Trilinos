#
# Unload the current Trilinos SEMS development environment
#
# USAGE:
#
#   $ source unload_sems_dev_env.sh
#
# If no Trilinos SEMS Dev Env is currently loaded, then this is a do nothing
# script.
#
# NOTE: If the modules don't unload correctly (e.g. because the
# sems-clang/<version> module loaded a GCC module), then consider calling:
#
#   $ module purge
#   $ unset TRILINOS_SEMS_DEV_ENV_LOADED
#
# to clean house (but this will also unload the sems-env module which will
# result in none of the other sems-xxx modules being visible).
#

# Get the base dir for the sourced script
called=$_
#[[ $called != $0 ]] && echo "Script is being sourced" || echo "Script is being run"
#echo "\$BASH_SOURCE ${BASH_SOURCE[@]}"
_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "_SCRIPT_DIR = '$_SCRIPT_DIR'"

source $_SCRIPT_DIR/std/sems/get_default_modules.sh

#
# A) Exit if no Trilinos SEMS Dev Env is currently loaded
#

if [ "$TRILINOS_SEMS_DEV_ENV_LOADED" == "" ] ; then
  if [ "${TRILINOS_SEMS_DEV_ENV_VERBOSE}" == "1" ] ; then
    echo "TRILINOS_SEMS_DEV_ENV_LOADED='', nothing to unload!"
  fi
  return 0
fi

#
# B) Get the loaded dev env info to unload the correct modules
#

idx=0
for str in $TRILINOS_SEMS_DEV_ENV_LOADED; do
  #echo $str
  sems_dev_env_array[$idx]=$str
  idx=$(expr $idx+1)
done
sems_compiler_and_version_unload=${sems_dev_env_array[0]}
sems_mpi_and_version_unload=${sems_dev_env_array[1]}
sems_cmake_and_version_unload=${sems_dev_env_array[2]}
#echo "sems_compiler_and_version_unload = $sems_compiler_and_version_unload"
#echo "sems_mpi_and_version_unload = $sems_mpi_and_version_unload"
if [ "$TRILINOS_SEMS_DEV_ENV_VERBOSE" == "1" ] ; then
  echo "Unloading TRILINOS_SEMS_DEV_ENV_LOADED =" \
    "'$sems_compiler_and_version_unload" \
    "$sems_mpi_and_version_unload" \
    "$sems_cmake_and_version_unload'!"
fi

#
# C) Unload the modules (in the reverse order)
#

module unload $sems_superlu_and_version_default
#module unload $sems_scotch_and_version_default
module unload $sems_parmetis_and_version_default
module unload $sems_netcdf_and_version_default
module unload $sems_hdf5_and_version_default
module unload $sems_zlib_and_version_default
module unload $sems_boost_and_version_default
module unload $sems_mpi_and_version_unload
module unload $sems_compiler_and_version_unload
module unload $sems_cmake_and_version_unload
module unload $sems_python_and_version_default

if [ "${TRILINOS_SEMS_DEV_ENV_VERBOSE}" == "1" ] ; then
  module list
fi

#
# D) Wipe out the rememberd Trilinos SEMS Dev Env
#

export TRILINOS_SEMS_DEV_ENV_LOADED=
