#!/bin/sh
TRILINOS_PATH=/Users/rouson/working-trilinos/Trilinos
EXTRA_ARGS=$@

rm -f CMakeCache.txt
#
# nagfor.pl is a wrapper workaround for a problem with preprocessor.
# When the vendor fixes the bug you can revert to nagfor.
# Works with cmake-2.8.3-g67a4d or later
#
cmake \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG \
  -D CMAKE_CXX_COMPILER:FILEPATH=g++ \
  -D CMAKE_C_COMPILER:FILEPATH=gcc \
  -D CMAKE_Fortran_COMPILER:FILEPATH="$TRILINOS_PATH/packages/ForTrilinos/nagfor.pl" \
  -D CMAKE_Fortran_FLAGS:STRING="-f2003 -g -C=all" \
  -D HAVE_GCC_ABI_DEMANGLE:BOOL=ON \
  -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
  -D Trilinos_ENABLE_DEFAULT_PACKAGES:BOOL=OFF \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
  -D Trilinos_ENABLE_CTrilinos:BOOL=ON\
  -D Trilinos_ENABLE_ForTrilinos:BOOL=ON\
  -D Trilinos_ENABLE_AztecOO:BOOL=ON\
  -D ForTrilinos_ENABLE_AztecOO:BOOL=ON \
  -D ForTrilinos_ENABLE_TESTS:BOOL=ON \
  -D ForTrilinos_ENABLE_EXTENDED:BOOL=ON \
  -D ForTrilinos_ENABLE_EXAMPLES:BOOL=ON \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
$EXTRA_ARGS \
$TRILINOS_PATH



