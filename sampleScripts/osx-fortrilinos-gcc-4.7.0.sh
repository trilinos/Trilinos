#!/bin/sh
TRILINOS_PATH=/Users/rouson/Trilinos.base/Trilinos
EXTRA_ARGS=$@

rm -f CMakeCache.txt
 
cmake \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D MPI_BASE_DIR:PATH="/usr/local/openmpi/GCC47" \
  -D MPI_BIN_DIR:PATH="/usr/local/openmpi/GCC47/bin" \
  -D MPI_USE_COMPILER_WRAPPERS:BOOL=ON \
  -D MPI_Fortran_COMPILER:FILEPATH="/usr/local/openmpi/GCC47/bin/mpif90" \
  -D MPI_CXX_COMPILER:FILEPATH="/usr/local/openmpi/GCC47/bin/mpicxx" \
  -D MPI_C_COMPILER:FILEPATH="/usr/local/openmpi/GCC47/bin/mpicc" \
  -D CMAKE_Fortran_FLAGS:STRING=" " \
  -D HAVE_GCC_ABI_DEMANGLE:BOOL=ON \
  -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
  -D Trilinos_ENABLE_CTrilinos:BOOL=ON\
  -D Trilinos_ENABLE_ForTrilinos:BOOL=ON\
  -D Trilinos_ENABLE_AztecOO:BOOL=ON\
  -D ForTrilinos_ENABLE_AztecOO:BOOL=ON \
  -D ForTrilinos_ENABLE_TESTS:BOOL=ON \
  -D ForTrilinos_ENABLE_OBJECT_ORIENTED:BOOL=ON \
  -D ForTrilinos_ENABLE_EXAMPLES:BOOL=ON \
  -D ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS:BOOL=ON \
  -D ForTrilinos_DISABLE_FINAL_SUBROUTINES:BOOL=ON \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
$EXTRA_ARGS \
$TRILINOS_PATH
