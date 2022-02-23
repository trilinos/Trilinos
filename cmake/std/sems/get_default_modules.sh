PLATFORM_TYPE=`uname`

if [ "$PLATFORM_TYPE" = "Linux" ] ; then
  sems_compiler_and_version_default=sems-archive-gcc/7.2.0
elif [ "$PLATFORM_TYPE" = "Darwin" ] ; then
  sems_compiler_and_version_default=sems-archive-gcc/5.3.0
else
  echo "ERROR, unknown platform type '$PLATFORM_TYPE'!"
  exit 1
fi

sems_mpi_and_version_default=sems-archive-openmpi/1.10.1
sems_python_and_version_default=sems-archive-python/2.7.9
sems_cmake_and_version_default=sems-archive-cmake/3.17.1
sems_ninja_and_version_default=sems-archive-ninja_fortran/1.8.2
sems_git_and_version_default=sems-archive-git/2.10.1
sems_boost_and_version_default=sems-archive-boost/1.63.0/base
sems_zlib_and_version_default=sems-archive-zlib/1.2.8/base
sems_hdf5_and_version_default=sems-archive-hdf5/1.10.6/parallel
sems_netcdf_and_version_default=sems-archive-netcdf/4.7.3/parallel
sems_parmetis_and_version_default=sems-archive-parmetis/4.0.3/parallel
sems_scotch_and_version_default=sems-archive-scotch/6.0.3/nopthread_64bit_parallel
sems_superlu_and_version_default=sems-archive-superlu/4.3/base

# NOTE: The above defaults are what are used for the standard CI dev env so
# changing these will change the standard CI env.
