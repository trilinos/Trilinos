#!/bin/bash
# Requires cmake version > 3.12
# Paths to source
KOKKOS_PATH="${HOME}/Kokkos/kokkos" #path to kokkos source
KOKKOSKERNELS_PATH="../.."          #path to kokkos-kernels top directory

# Compiler - must be passed to kokkos and kokkos-kernels configurations
CXX=${KOKKOS_PATH}/bin/nvcc_wrapper #Options: icpc #g++ #clang++
CXXFLAGS="-Wall -pedantic -Werror -O3 -g -Wshadow -Wsign-compare -Wignored-qualifiers -Wempty-body -Wclobbered -Wuninitialized"

# Configure Kokkos (Unit Tests OFF) - Makefile located in kokkos-build
cmake -Bkokkos-build -DCMAKE_CXX_COMPILER=${CXX} -DKokkos_ARCH_PASCAL60=ON -DKokkos_ARCH_POWER8=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -DCMAKE_CXX_FLAGS="${CXXFLAGS}" -DCMAKE_INSTALL_PREFIX="${PWD}/kokkos-install" -DKokkos_ENABLE_TESTS=OFF ${KOKKOS_PATH}

# Build and Install Kokkos - install lib at ${PWD}/kokkos-install
cmake --build kokkos-build -j 8 --target install


# Configure KokkosKernels (Unit Tests OFF) - Makefile located in kokkoskernels-build
cmake -Bkokkoskernels-build -DCMAKE_CXX_COMPILER=${CXX} -DKokkos_ROOT="${PWD}/kokkos-install" -DKokkosKernels_INST_DOUBLE=ON -DKokkosKernels_INST_COMPLEX_DOUBLE=ON -DKokkosKernels_INST_ORDINAL_INT=ON -DKokkosKernels_INST_ORDINAL_INT64_T=ON -DKokkosKernels_INST_OFFSET_INT=ON -DKokkosKernels_INST_OFFSET_SIZE_T=ON -DKokkosKernels_INST_LAYOUTLEFT=ON -DKokkosKernels_ADD_DEFAULT_ETI=ON -DCMAKE_INSTALL_PREFIX="${PWD}/kokkoskernels-install" -DKokkosKernels_ENABLE_TESTS=OFF -DKokkosKernels_ENABLE_TPL_CUBLAS=OFF ${KOKKOSKERNELS_PATH}

# Build and Install KokkosKernels - install lib at ${PWD}/kokkoskernels-install
cmake --build kokkoskernels-build -j 8 --target install
