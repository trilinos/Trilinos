#!/bin/bash
#SBATCH -A cfd116
#SBATCH -J configureAndBuildTrilinos
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END

TRILINOS_SRC="$HOME/spock/sources/trilinos/Trilinos-muelu"
TRILINOS_BUILD="$HOME/spock/builds/performance"
BUILD_TYPE="RELEASE"

#proxies to allow git operations from compute nodes
export all_proxy=socks://proxy.ccs.ornl.gov:3128/
export ftp_proxy=ftp://proxy.ccs.ornl.gov:3128/
export http_proxy=http://proxy.ccs.ornl.gov:3128/
export https_proxy=http://proxy.ccs.ornl.gov:3128/
export no_proxy='localhost,127.0.0.0/8,*.ccs.ornl.gov,*.ncrc.gov',*.olcf.ornl.gov

source "${TRILINOS_SRC}/packages/tpetra/scripts/PerformanceTesting/spock/load_modules.sh"
module list

cd ${TRILINOS_SRC}
git fetch
git checkout develop
git merge

cd $TRILINOS_BUILD

echo "SOURCE = $TRILINOS_SRC"
echo "BUILD  = $PWD"
echo "BUILD_TYPE = $BUILD_TYPE"

rm -rf CMakeCache.txt CMakeFiles

ARGS=(
  -GNinja
  -D Trilinos_ENABLE_Fortran=OFF
  -D Trilinos_TEST_CATEGORIES=PERFORMANCE
  -D Trilinos_ENABLE_FORTRAN=OFF
  -D TPL_ENABLE_HDF5=OFF
  -D   HDF5_INCLUDE_DIRS:PATH="${CRAY_HDF5_PARALLEL_DIR}/include"
  -D   HDF5_LIBRARY_DIRS:PATH="${CRAY_HDF5_PARALLEL_DIR}/lib"
  -D TPL_ENABLE_Netcdf=OFF
  -D TPL_ENABLE_Matio=OFF
  -D TPL_ENABLE_X11=OFF
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON
  -D MPI_EXEC_MAX_NUMPROCS=16
  -D Trilinos_AUTOGENERATE_TEST_RESOURCE_FILE=OFF
  -D Trilinos_ENABLE_PyTrilinos=OFF
  -D CMAKE_BUILD_TYPE:STRING="$BUILD_TYPE"
  -D CMAKE_CXX_COMPILER:PATH="/opt/rocm-4.5.0/bin/hipcc"
  -D CMAKE_CXX_FLAGS:STRING="-I${CRAY_MPICH_DIR}/include"
  -D BUILD_SHARED_LIBS:BOOL=OFF
  -D Trilinos_EXTRA_LINK_FLAGS:STRING="-L${CRAY_MPICH_DIR}/lib -lmpi"
  -D Trilinos_ENABLE_TESTS:BOOL=OFF
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
  -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF
  -D Trilinos_ALLOW_NO_PACKAGES:BOOL=OFF
  -D Trilinos_ENABLE_OpenMP:BOOL=OFF
  -D Trilinos_ENABLE_COMPLEX_DOUBLE:BOOL=OFF

  -D Trilinos_ENABLE_Amesos:BOOL=ON
  -D Trilinos_ENABLE_Amesos2:BOOL=ON
  -D   Amesos2_ENABLE_SuperLU:BOOL=OFF
  -D   Amesos2_ENABLE_KLU2:BOOL=ON
  -D   Amesos2_ENABLE_Basker:BOOL=ON
  -D Trilinos_ENABLE_AztecOO:BOOL=ON
  -D Trilinos_ENABLE_Belos:BOOL=ON
  -D Trilinos_ENABLE_Epetra:BOOL=ON
  -D Trilinos_ENABLE_Galeri:BOOL=ON
  -D Trilinos_ENABLE_Ifpack:BOOL=ON
  -D Trilinos_ENABLE_Ifpack2:BOOL=ON
    -D Ifpack2_ENABLE_TESTS:BOOL=ON
  -D Trilinos_ENABLE_Kokkos:BOOL=ON
  -D   Kokkos_ARCH_ZEN2:BOOL=ON
  -D   Kokkos_ARCH_VEGA908:BOOL=ON
  -D   Kokkos_ENABLE_SERIAL:BOOL=ON
  -D   Kokkos_ENABLE_HIP:BOOL=ON
  -D Trilinos_ENABLE_KokkosKernels:BOOL=ON
  -D   KokkosKernels_ENABLE_SUPERNODAL_SPTRSV:BOOL=OFF
  -D   KokkosKernels_ENABLE_TPL_SUPERLU:BOOL=OFF
  -D Trilinos_ENABLE_ML:BOOL=ON
  -D Trilinos_ENABLE_MueLu:BOOL=ON
  -D   MueLu_ENABLE_Kokkos_Refactor:BOOL=ON
  -D   MueLu_ENABLE_TESTS:BOOL=ON
  -D   MueLu_ENABLE_Epetra:BOOL=OFF
  -D   MueLu_ENABLE_EXAMPLES:BOOL=ON
  -D Trilinos_ENABLE_Panzer:BOOL=ON
    -D Trilinos_ENABLE_PanzerMiniEM:BOOL=ON
    -D PanzerMiniEM_ENABLE_EXAMPLES:BOOL=ON
    -D PanzerMiniEM_ENABLE_TESTS:BOOL=ON
  -D Trilinos_ENABLE_Percept:BOOL=OFF
  -D Trilinos_ENABLE_Teko:BOOL=ON
  -D Trilinos_ENABLE_Tpetra:BOOL=ON
  -D   Tpetra_ENABLE_TESTS:BOOL=ON
  -D   Tpetra_ENABLE_EXAMPLES:BOOL=ON
  -D   Tpetra_ENABLE_CUDA:BOOL=OFF
  -D   Tpetra_INST_SERIAL:BOOL=ON
  -D   Tpetra_INST_OPENMP:BOOL=OFF
  -D   Tpetra_INST_HIP:BOOL=ON
  -D   Tpetra_INST_DOUBLE:BOOL=ON
  -D Trilinos_ENABLE_Teuchos:BOOL=ON
  -D   Teuchos_KOKKOS_PROFILING:BOOL=ON
  -D Trilinos_ENABLE_Xpetra:BOOL=ON
  -D   Xpetra_ENABLE_Kokkos_Refactor:BOOL=ON
  -D   Xpetra_ENABLE_Epetra:BOOL=OFF
  -D Trilinos_ENABLE_Zoltan2:BOOL=ON
  -D TPL_ENABLE_BLAS:BOOL=ON
  -D   BLAS_LIBRARY_DIRS:PATH="/ccs/home/jjhu/spock/install/lib64"
  -D TPL_ENABLE_Boost:BOOL=OFF
  -D   Boost_INCLUDE_DIRS:PATH="${OLCF_BOOST_ROOT}/include"
  -D   Boost_LIBRARY_DIRS:PATH="${OLCF_BOOST_ROOT}/lib"
  -D TPL_ENABLE_LAPACK:BOOL=ON
  -D   LAPACK_LIBRARY_DIRS:PATH="/ccs/home/jjhu/spock/install/lib64"
  -D TPL_ENABLE_MPI:BOOL=ON
  -D   MPI_BASE_DIR:PATH="${CRAY_MPICH_DIR}"
  -D   MPI_USE_COMPILER_WRAPPERS:BOOL=OFF
  -D   MPI_EXEC_NUMPROCS_FLAG:STRING="--ntasks"
  -D   MPI_EXEC:STRING="srun"
  -D TPL_ENABLE_SuperLU:BOOL=OFF
  -D   SuperLU_INCLUDE_DIRS:PATH="/gpfs/alpine/scratch/lberge/csc465/installs/superlu/include"
  -D   SuperLU_LIBRARY_DIRS:PATH="/gpfs/alpine/scratch/lberge/csc465/installs/superlu/lib64"
  -D TPL_ENABLE_Zlib:BOOL=ON
  -D   Zlib_INCLUDE_DIRS:PATH="${OLCF_ZLIB_ROOT}/include"
  -D   Zlib_LIBRARY_DIRS:PATH="${OLCF_ZLIB_ROOT}/lib"
)

cmake "${ARGS[@]}" $TRILINOS_SRC |& tee configure_trilinos_performance.log

NUMTHREADS=127
numactl -C $(seq -s, 0 2 $NUMTHREADS) ninja -j $(seq -s, 0 2 $NUMTHREADS | tr ',' '\n' | wc -l)
