#!/bin/bash

export WATCHR_BUILD_NAME="Stria Serial"
source $TRILINOS_SRC/cmake/std/atdm/load-env.sh stria-arm-release

#Update Trilinos (should always be on develop branch)
cd $TRILINOS_SRC
git fetch origin
git reset --hard origin/develop
export TRILINOS_GIT_SHA=`git rev-parse HEAD`

cd $WORKSPACE/build

rm -rf CMakeFiles
rm -f CMakeCache.txt

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
-DMPI_EXEC_MAX_NUMPROCS=36 \
-DTrilinos_ENABLE_PyTrilinos=OFF \
$TRILINOS_SRC

make -j36

ctest || true

