#!/bin/bash
#SBATCH -A csc465
#SBATCH -J configureAndBuildTrilinos
#SBATCH -t 05:00:00
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END

# ExaWind
# #SBATCH -A cfd116

#TRILINOS_SRC="$HOME/crusher/sources/trilinos/Trilinos-muelu"
# FIXME JHU TEMPORARY!!! 2022-07-22
TRILINOS_SRC="$HOME/crusher/sources/trilinos/Trilinos"
# FIXME JHU TEMPORARY!!!
TRILINOS_BUILD="$HOME/crusher/builds/performance-new"
BUILD_TYPE="RELEASE"

#proxies to allow git operations from compute nodes
export all_proxy=socks://proxy.ccs.ornl.gov:3128/
export ftp_proxy=ftp://proxy.ccs.ornl.gov:3128/
export http_proxy=http://proxy.ccs.ornl.gov:3128/
export https_proxy=http://proxy.ccs.ornl.gov:3128/
export no_proxy='localhost,127.0.0.0/8,*.ccs.ornl.gov,*.ncrc.gov',*.olcf.ornl.gov


#source "${TRILINOS_SRC}/packages/tpetra/scripts/PerformanceTesting/crusher/load_modules.sh"
# FIXME JHU TEMPORARY!!!
source "$HOME/crusher/crusher_scripts/load_modules.sh"
module list

if [[ "1" == "0" ]]; then

cd ${TRILINOS_SRC}
#git fetch
#git checkout develop
#git merge
gitCmds=(
  "git fetch"
  "git checkout develop"
  "git merge"
)

for gc in "${gitCmds[@]}"; do
  echo "${gc[@]}"
  realgc=($gc)
  "${realgc[@]}"
  if [[ "$?" == "0" ]]; then
    echo "  \"${realgc[@]}\" succeeded"
  else
    echo "  \"${realgc[@]}\" failed"
  fi
done
fi

cd ${TRILINOS_BUILD}

echo "SOURCE = $TRILINOS_SRC"
echo "BUILD  = $PWD"
echo "BUILD_TYPE = $BUILD_TYPE"

rm -Rf CMakeCache.txt CMakeFiles

# See https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html#compiling for tips on compiling.

#export HIPFLAGS="--amdgpu-target=gfx90a"
#PE_MPICH_GTL_DIR_amd_gfx90a="-L${CRAY_MPICH_ROOTDIR}/gtl/lib"
#PE_MPICH_GTL_LIBS_amd_gfx90a="-lmpi_gtl_hsa"

THE_REAL_NETCDF_C_ROOT="${CRAY_NETCDF_DIR}/crayclang/10.0"

CMAKE_OPTIONS="--debug-trycompile"
#CMAKE_OPTIONS=""

ARGS=(
    -GNinja
    -D Trilinos_TEST_CATEGORIES=PERFORMANCE
    -D Trilinos_ENABLE_Fortran:BOOL=OFF
    -D CMAKE_INSTALL_PREFIX:PATH="${trilinos_install_dir}"
    -D CMAKE_BUILD_TYPE:STRING=$BUILD_TYPE



    #-D MPI_USE_COMPILER_WRAPPERS:BOOL=ON
    # JHU: 2022-02-07 I've been unable to get Trilinos to build with the Cray MPI wrappers.
    # Instead, compile with hipcc directly.
    -D MPI_USE_COMPILER_WRAPPERS:BOOL=OFF
    -D MPI_EXEC_NUMPROCS_FLAG:STRING="--ntasks"
    -D MPI_EXEC:STRING="srun"
    #-D CMAKE_CXX_FLAGS:STRING="-I${ROCM_PATH}/include"
    -D Trilinos_EXTRA_LINK_FLAGS:STRING="-L${ROCM_PATH}/lib -lamdhip64"
    #-D CMAKE_CXX_COMPILER:PATH="$HIP_PATH/bin/hipcc"
    -D CMAKE_CXX_COMPILER:PATH="$HIP_CLANG_PATH/../../bin/hipcc"
    -D CMAKE_CXX_FLAGS:STRING="-I${CRAY_MPICH_DIR}/include"
    #-D CMAKE_C_COMPILER:PATH="$HIP_PATH/bin/hipcc"
    -D CMAKE_C_COMPILER:PATH="$HIP_CLANG_PATH/../../bin/hipcc"
    -D CMAKE_C_FLAGS:STRING="-I${CRAY_MPICH_DIR}/include"
    -D Trilinos_EXTRA_LINK_FLAGS:STRING="-L${CRAY_MPICH_DIR}/lib -lmpi -L${CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa -L/${CRAY_BINUTILS_VERSION}/cce/x86_64/lib -lcraymath -lu -lf"

    -D BUILD_SHARED_LIBS:BOOL=ON


    -D Trilinos_ENABLE_TESTS:BOOL=OFF
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
    -D Trilinos_ENABLE_CXX11:BOOL=ON
    -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
    -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF
    -D Trilinos_ALLOW_NO_PACKAGES:BOOL=OFF
    -D Trilinos_ENABLE_OpenMP:BOOL=OFF

    -D Trilinos_ENABLE_Amesos:BOOL=ON
    -D Trilinos_ENABLE_Amesos2:BOOL=ON
    -D   Amesos2_ENABLE_SuperLU:BOOL=OFF
    -D   Amesos2_ENABLE_KLU2:BOOL=ON
    -D Trilinos_ENABLE_AztecOO:BOOL=ON
    -D Trilinos_ENABLE_Belos:BOOL=ON
    -D Trilinos_ENABLE_Epetra:BOOL=ON
    -D Trilinos_ENABLE_EpetraExt:BOOL=ON
    -D Trilinos_ENABLE_Galeri:BOOL=ON
    -D Trilinos_ENABLE_Ifpack:BOOL=ON
    -D Trilinos_ENABLE_Ifpack2:BOOL=ON
      -D Ifpack2_ENABLE_TESTS:BOOL=ON
    -D Trilinos_ENABLE_Kokkos:BOOL=ON
    -D   Kokkos_ENABLE_OPENMP:BOOL=OFF
    -D   Kokkos_ENABLE_HIP:BOOL=ON
      -D   Kokkos_ARCH_ZEN3:BOOL=ON
      -D   Kokkos_ARCH_VEGA90A:BOOL=ON
    -D   Kokkos_ENABLE_DEPRECATED_CODE:BOOL=OFF
    -D Trilinos_ENABLE_KokkosKernels:BOOL=ON
    -D Trilinos_ENABLE_ML:BOOL=ON
    -D Trilinos_ENABLE_MueLu:BOOL=ON
       -D MueLu_ENABLE_Epetra:BOOL=OFF
       -D MueLu_ENABLE_Kokkos_Refactor:BOOL=ON
       -D MueLu_ENABLE_TESTS:BOOL=ON
       -D MueLu_ENABLE_EXAMPLES:BOOL=ON
    -D Trilinos_ENABLE_Panzer:BOOL=ON
      -D Trilinos_ENABLE_PanzerMiniEM:BOOL=ON
      -D PanzerMiniEM_ENABLE_EXAMPLES:BOOL=ON
      -D PanzerMiniEM_ENABLE_TESTS:BOOL=ON
    -D Trilinos_ENABLE_Percept=OFF
    -D Trilinos_ENABLE_Teko:BOOL=ON
    -D Trilinos_ENABLE_Tpetra:BOOL=ON
    -D   Tpetra_ENABLE_EXAMPLES:BOOL=ON
    -D   Tpetra_ENABLE_TESTS:BOOL=ON
    -D   Tpetra_INST_HIP:BOOL=ON
    -D   Tpetra_INST_SERIAL:BOOL=ON
    -D   Tpetra_INST_OPENMP:BOOL=OFF
    -D   Tpetra_INST_DOUBLE:BOOL=ON
    -D Trilinos_ENABLE_Gtest:BOOL=OFF
    -D Trilinos_ENABLE_Teuchos:BOOL=ON
    -D   Teuchos_KOKKOS_PROFILING:BOOL=ON
    -D Trilinos_ENABLE_Xpetra:BOOL=ON
    -D   Xpetra_ENABLE_Kokkos_Refactor:BOOL=ON
    -D   Xpetra_ENABLE_Epetra:BOOL=OFF
    -D Trilinos_ENABLE_Zoltan:BOOL=ON
    -D Trilinos_ENABLE_Zoltan2:BOOL=ON
    -D TPL_ENABLE_BLAS:BOOL=ON
        -D TPL_BLAS_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_cray.a
    -D TPL_ENABLE_LAPACK:BOOL=ON
        -D TPL_LAPACK_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_cray.a
    -D TPL_ENABLE_Boost:BOOL=OFF
    #-D   Boost_INCLUDE_DIRS:PATH="${OLCF_BOOST_ROOT}/include"
    #-D   Boost_LIBRARY_DIRS:PATH="${OLCF_BOOST_ROOT}/lib"
    #-D   BoostLib_INCLUDE_DIRS:PATH="${OLCF_BOOST_ROOT}/include"
    #-D   BoostLib_LIBRARY_DIRS:PATH="${OLCF_BOOST_ROOT}/lib"
    #-D TPL_ENABLE_HDF5:BOOL=ON
    #-D   HDF5_LIBRARY_DIRS:PATH="${CRAY_HDF5_DIR}/CRAYCLANG/10.0/lib"
    #-D   HDF5_INCLUDE_DIRS:PATH="${CRAY_HDF5_DIR}/CRAYCLANG/10.0/include"
    #-D   HDF5_LIBRARY_DIRS:PATH="${PE_HDF5_PARALLEL_DIR}/CRAYCLANG/10.0/lib"
    #-D   HDF5_INCLUDE_DIRS:PATH="${PE_HDF5_PARALLEL_DIR}/CRAYCLANG/10.0/include"
    #-D   HDF5_NO_SYSTEM_PATHS:BOOL=ON
    -D TPL_ENABLE_MPI:BOOL=ON
    #-D TPL_ENABLE_Pnetcdf:BOOL=ON
    #-D   PNetcdf_INCLUDE_DIRS:PATH="${CRAY_PARALLEL_NETCDF_ROOT}/include"
    #-D   PNetcdf_LIBRARIES:PATH="${CRAY_PARALLEL_NETCDF_ROOT}/lib/libpnetcdf.a"
    #-D TPL_ENABLE_Netcdf:BOOL=ON
    #   -D Netcdf_INCLUDE_DIRS:PATH="${CRAY_NETCDF_HDF5PARALLEL_DIR}/crayclang/10.0/include"
    #   -D Netcdf_LIBRARY_DIRS:PATH="${CRAY_NETCDF_HDF5PARALLEL_DIR}/crayclang/10.0/lib"
    #-D TPL_Netcdf_PARALLEL:BOOL=ON
    #-D TPL_Netcdf_Enables_Netcdf4:BOOL=ON
    #-D TPL_ENABLE_METIS:BOOL=ON
    -D TPL_ENABLE_HDF5=OFF
    -D TPL_ENABLE_Netcdf=OFF
    #-D   METIS_INCLUDE_DIRS:PATH="${OLCF_METIS_ROOT}/include"
    #-D   METIS_LIBRARY_DIRS:PATH="${OLCF_METIS_ROOT}/lib"
    #-D TPL_ENABLE_ParMETIS:BOOL=ON
    #-D   ParMETIS_INCLUDE_DIRS:PATH="${OLCF_PARMETIS_ROOT}/include"
    #-D   ParMETIS_LIBRARY_DIRS:PATH="${OLCF_PARMETIS_ROOT}/lib"
    -D TPL_ENABLE_SuperLU:BOOL=OFF
    #-D   SuperLU_INCLUDE_DIRS:PATH="${OLCF_SUPERLU_ROOT}/include"
    #-D   SuperLU_LIBRARY_DIRS:PATH="${OLCF_SUPERLU_ROOT}/lib"
    -D TPL_ENABLE_Zlib:BOOL=ON
    -D   Zlib_INCLUDE_DIRS:PATH="${OLCF_ZLIB_ROOT}"/include
    -D   Zlib_LIBRARY_DIRS:PATH="${OLCF_ZLIB_ROOT}"/lib
)

CMAKE_START_TIME=`date`
cmake ${CMAKE_OPTIONS} "${ARGS[@]}" $TRILINOS_SRC |& tee configure_trilinos_performance.log
CMAKE_END_TIME=`date`

BUILD_START_TIME=`date`
NUMTHREADS=31
numactl -C $(seq -s, 0 2 $NUMTHREADS) ninja -j $(seq -s, 0 2 $NUMTHREADS | tr ',' '\n' | wc -l) |& tee trilinos_build.log
BUILD_END_TIME=`date`

echo "cmake start: $CMAKE_START_TIME"
echo "cmake end: $CMAKE_END_TIME"
echo ""
echo "build start: $BUILD_START_TIME"
echo "build end: $BUILD_END_TIME"
