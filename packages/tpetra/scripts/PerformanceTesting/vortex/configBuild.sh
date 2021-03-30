#!/bin/bash

source $TRILINOS_SRC/cmake/std/atdm/load-env.sh cuda-10.1.243_gnu-7.3.1_spmpi-2019.06.24-release

echo "Switching to develop branch."
cd $TRILINOS_SRC
git fetch origin
git reset --hard origin/develop

cd $WORKSPACE/build

rm -rf CMakeFiles
rm -f CMakeCache.txt

#NOTE: the options Trilinos_CUDA_NUM_GPUS and Trilinos_CUDA_SLOTS_PER_GPU override the ATDM
#defaults. This is required to make ctest use all GPUs of 4 nodes.

cmake \
-DCMAKE_BUILD_TYPE=Release \
-DTrilinos_TEST_CATEGORIES=PERFORMANCE \
-DTrilinos_ENABLE_TESTS=ON \
-DTrilinos_ENABLE_EXAMPLES=ON \
-DTrilinos_ENABLE_Tpetra=ON \
-DTrilinos_ENABLE_MueLu=ON \
-DTrilinos_ENABLE_PanzerMiniEM=ON \
-DTrilinos_ENABLE_Percept=OFF \
-DTpetra_INST_SERIAL=ON \
-DTpetra_INST_INT_INT=ON \
-DXpetra_ENABLE_Epetra=ON \
-DMueLu_ENABLE_Epetra=ON \
-DTPL_ENABLE_HDF5=ON \
-DTPL_ENABLE_Netcdf=ON \
-DTPL_ENABLE_Matio=OFF \
-DTPL_ENABLE_X11=OFF \
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=$TRILINOS_SRC/cmake/std/atdm/ATDMDevEnv.cmake \
-DTPL_ENABLE_MPI=ON \
-DKokkos_ENABLE_SERIAL=ON \
-DMPI_EXEC_MAX_NUMPROCS=16 \
-DTrilinos_AUTOGENERATE_TEST_RESOURCE_FILE=OFF \
-DTrilinos_ENABLE_PyTrilinos=OFF \
$TRILINOS_SRC

make -j36
