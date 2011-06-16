#!/bin/bash
rm -rf CMakeCache.txt
cmake \
      -D CMAKE_INSTALL_PREFIX="/home/rppawlo/install_seacas" \
      -D Trilinos_EXTRA_REPOSITORIES="Panzer" \
      -D Trilinos_ENABLE_DEBUG=OFF \
      -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
      -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
      -D Trilinos_ENABLE_TESTS:BOOL=OFF \
      -D Trilinos_ENABLE_Pamgen:BOOL=OFF \
      -D Trilinos_ENABLE_SEACAS:BOOL=ON \
      -D SEACAS_ENABLE_EXODUS:BOOL=ON \
      -D SEACAS_ENABLE_APREPRO:BOOL=OFF \
      -D SEACAS_ENABLE_APPLICATIONS:BOOL=ON \
      -D TPL_ENABLE_MPI:BOOL=OFF \
      -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_INCLUDE_DIRS:FILEPATH="/home/rppawlo/install_netcdf/include" \
      -D Netcdf_LIBRARY_DIRS:FILEPATH="/home/rppawlo/install_netcdf/lib" \
      -D CMAKE_CXX_COMPILER:FILEPATH="g++" \
      -D CMAKE_C_COMPILER:FILEPATH="gcc" \
      -D Trilinos_EXTRA_LINK_FLAGS:STRING="-L/usr/lib -lm" \
      -D Trilinos_ENABLE_Fortran:BOOL=OFF \
      -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
      -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
      -D CMAKE_SKIP_RULE_DEPENDENCY=ON \
      -D Trilinos_ENABLE_STRONG_CXX_COMPILE_WARNINGS=OFF \
      -D Trilinos_ENABLE_STRONG_C_COMPILE_WARNINGS=OFF \
      -D Trilinos_ENABLE_SHADOW_WARNINGS=OFF \
       ../Trilinos
