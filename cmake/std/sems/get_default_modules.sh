PLATFORM_TYPE=`uname`

if [ "$PLATFORM_TYPE" = "Linux" ] ; then
  sems_compiler_and_version_default=sems-gcc/4.8.4
elif [ "$PLATFORM_TYPE" = "Darwin" ] ; then
  sems_compiler_and_version_default=sems-gcc/5.3.0
else
  echo "ERROR, unknown platform type '$PLATFORM_TYPE'!"
  exit 1
fi

sems_mpi_and_version_default=sems-openmpi/1.10.1
sems_python_and_version_default=sems-python/2.7.9
sems_cmake_and_version_default=atdm-cmake/3.11.1
sems_git_and_version_default=sems-git/2.10.1
sems_boost_and_version_default=sems-boost/1.63.0/base
sems_zlib_and_version_default=sems-zlib/1.2.8/base 
sems_hdf5_and_version_default=sems-hdf5/1.8.12/parallel 
sems_netcdf_and_version_default=sems-netcdf/4.4.1/exo_parallel
sems_parmetis_and_version_default=sems-parmetis/4.0.3/parallel 
sems_scotch_and_version_default=sems-scotch/6.0.3/nopthread_64bit_parallel   
sems_superlu_and_version_default=sems-superlu/4.3/base

# NOTE: The above defaults are what are used for the standard CI dev env so
# changing these will change the standard CI env.
