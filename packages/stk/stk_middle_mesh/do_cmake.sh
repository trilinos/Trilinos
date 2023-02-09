#!/bin/bash

BUILD_SUFFIX=dbg  # opt or dbg or asan
# Note: if building with Trilinos, compilers should be MPI wrappers

#TRILINOS_BASE is an environment variable
#LINK_SUFFIX=shared
#TRILINOS=${TRILINOS_BASE}_serial_${LINK_SUFFIX}_${BUILD_SUFFIX}
#_legacy_trilinos=${TRILINOS_BASE}_${LINK_SUFFIX}_${BUILD_SUFFIX}
#if [[ ! -d ${TRILINOS} ]]
#then
#  echo -e "REMARK: Trilinos directory\n'${TRILINOS}'\ndoes not exist; using legacy-style directory:\n'${_legacy_trilinos}'"
#  TRILINOS=${_legacy_trilinos}
#fi

echo "TRILINOS_ROOT = $TRILINOS_ROOT"
# -lgfortran is required by OpenBlas (adding this flag to its CMake invocation didn't work)
if [[ $BUILD_SUFFIX == "dbg" ]]
then
  cmake -D CMAKE_C_COMPILER=`which mpicc` \
        -D CMAKE_CXX_COMPILER=`which mpicxx` \
        -D CMAKE_CXX_FLAGS="-Wall -g -O3 -fno-fast-math" \
        -D CMAKE_BUILD_TYPE=Debug \
        -D CMAKE_EXE_LINKER_FLAGS="-lgfortran" \
        -D CBLAS_INCLUDE_DIRS="$OPENBLAS_ROOT"/include/openblas \
        -D CBLAS_LIBRARY_DIRS="$OPENBLAS_ROOT"/lib64 \
        -D CBLAS_LIBRARIES="libopenblas.a" \
        -D Trilinos_ROOT=${TRILINOS_ROOT} \
        -D CMAKE_EXPORT_COMPILE_COMMANDS=1 \
        -D ENABLE_STK=ON \
        -G Ninja \
        ..
elif [[ $BUILD_SUFFIX == "asan" ]]
then
  cmake -D CMAKE_CXX_COMPILER=`which mpicxx` \
        -D CMAKE_CXX_FLAGS="-Wall -g -fsanitize=address -fno-omit-frame-pointer -O0" \
        -D CMAKE_BUILD_TYPE=Debug \
        -D CMAKE_EXE_LINKER_FLAGS="-lgfortran" \
        -D CBLAS_INCLUDE_DIRS="$OPENBLAS_ROOT"/include/openblas \
        -D CBLAS_LIBRARY_DIRS="$OPENBLAS_ROOT"/lib64 \
        -D CBLAS_LIBRARIES="libopenblas.a" \
        -D Trilinos_ROOT=${TRILINOS_ROOT} \
        -D CMAKE_EXPORT_COMPILE_COMMANDS=1 \
        -D ENABLE_STK=OFF \
        -G Ninja \
        ..
else
  cmake -D CMAKE_CXX_COMPILER=`which mpicxx` \
        -D CMAKE_CXX_FLAGS="-Wall -O3" \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_EXE_LINKER_FLAGS="-lgfortran" \
        -D CBLAS_INCLUDE_DIRS="$OPENBLAS_ROOT"/include/openblas \
        -D CBLAS_LIBRARY_DIRS="$OPENBLAS_ROOT"/lib64 \
        -D CBLAS_LIBRARIES="libopenblas.a" \
        -D Trilinos_ROOT=${TRILINOS_ROOT} \
        -D CMAKE_EXPORT_COMPILE_COMMANDS=1 \
        -D ENABLE_STK=OFF \
        -G Ninja \
        ..
fi

#cmake -D CMAKE_CXX_FLAGS="-Wall -O3 -lgfortran" -D CMAKE_BUILD_TYPE=Release ..
