#
# TPLs and other compiler-specific settings for Trilinos GCC 4.6.1 builds
#
# NOTE: This file only needs to get included when doing an actual configure,
# not during package dependency adjustment such as in the CI and Nightly
# CTest/CDash drivers.
#

SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;hdf5_cpp" CACHE STRING "")
SET(HDF5_LIBRARY_DIRS /projects/gcc-4.6.1/tpls/hdf5-1.8.7/serial/lib CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS /projects/gcc-4.6.1/tpls/hdf5-1.8.7/serial/include CACHE FILEPATH "")
SET(Netcdf_INCLUDE_DIRS /opt/gcc-4.6.1/tpls/netcdf-4.2/include CACHE FILEPATH "")
SET(Netcdf_LIBRARY_DIRS /opt/gcc-4.6.1/tpls/netcdf-4.2/lib CACHE FILEPATH "")
