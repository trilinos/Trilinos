#!/bin/bash

source /projects/sems/modulefiles/utils/sems-modules-init.sh

module load aue/cmake/3.26.2
module load aue/gcc/10.3.0
module load aue/ninja/1.11.1
module load cde/v3/boost/1.79.0-gcc-10.3.0-openmpi-4.1.2
module load apps/cuda-11.7

# change this to point to your Trilinos source
TRILINOS_HOME=/scratch/kxren/Trilinos
cd $TRILINOS_HOME
#GIT_SHA=`git rev-parse HEAD | cut -c 1-8`
GIT_SHA=`git rev-parse HEAD`
# jumps back to previous directory
cd - > /dev/null 2>&1

# creates a "unique" log
CMAKE_LOG_FILE=`mktemp cmake.XXXXX.log`

export OMPI_CXX=${TRILINOS_HOME}/packages/kokkos/bin/nvcc_wrapper
BUILD_TYPE="RELWITHDEBINFO"


rm -rf CMakeFiles CMakeCache.txt

module -w1 list
echo "trilinos source: $TRILINOS_HOME" >> $CMAKE_LOG_FILE
echo "trilinos build : $PWD" >> $CMAKE_LOG_FILE
echo "build_type     : $BUILD_TYPE" >> $CMAKE_LOG_FILE
echo "git sha        : $GIT_SHA" >> $CMAKE_LOG_FILE

#sleep 5

ARGS=(
  -GNinja
  -D Trilinos_ENABLE_TESTS=OFF
  -D Trilinos_ENABLE_EXAMPLES=OFF
  -D CMAKE_BUILD_TYPE:STRING="$BUILD_TYPE"

  -D Trilinos_ENABLE_Kokkos=ON
   -D Kokkos_ARCH_VOLTA70=ON
   -D Kokkos_ENABLE_CUDA=ON
   -D Kokkos_ENABLE_CUDA_UVM=OFF
   -D Kokkos_ENABLE_CUDA_LAMBDA=ON
   -D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=OFF

  #-D Trilinos_ENABLE_MueLu=ON
  #  -D MueLu_ENABLE_TESTS=ON
  #  -D MueLu_ENABLE_EXAMPLES=ON

  -D Trilinos_ENABLE_Ifpack2=OFF
    -D Ifpack2_ENABLE_TESTS=ON
    -D Ifpack2_ENABLE_EXAMPLES=ON

  -D Trilinos_ENABLE_Tpetra=ON
   -D Tpetra_INST_SERIAL=OFF
   -D Tpetra_INST_CUDA=ON
   -D Tpetra_ENABLE_TESTS=ON
   -D Tpetra_ENABLE_EXAMPLES=ON

  -D Trilinos_ENABLE_ShyLU=OFF
    -D Trilinos_ENABLE_ShyLU_DD=OFF
    -D Trilinos_ENABLE_ShyLU_Node=OFF
    -D Trilinos_ENABLE_ShyLU_NodeTacho=OFF
  -D Trilinos_ENABLE_Sacado=OFF
  -D Trilinos_ENABLE_Shards=OFF
  -D Trilinos_ENABLE_Pamgen=OFF
  -D Trilinos_ENABLE_Anasazi=OFF
  -D Trilinos_ENABLE_Teko=OFF
  -D Trilinos_ENABLE_Intrepid2=OFF

  -D TPL_ENABLE_MPI=ON

  -D CMAKE_C_COMPILER:PATH=`which mpicc`
  -D CMAKE_CXX_COMPILER:PATH=`which mpicxx`

)

echo ""
echo "" >> $CMAKE_LOG_FILE
echo "${ARGS[@]}" | sed "s/-D /-D/g" | tr ' ' '\n' >> $CMAKE_LOG_FILE
echo "" >> $CMAKE_LOG_FILE
echo ""

#cmake --debug-trycompile "${ARGS[@]}" ${TRILINOS_HOME} 2>&1 | tee -a $CMAKE_LOG_FILE
cmake "${ARGS[@]}" ${TRILINOS_HOME} 2>&1 | tee -a $CMAKE_LOG_FILE
