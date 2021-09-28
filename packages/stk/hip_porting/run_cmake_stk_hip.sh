#!/bin/bash

rm -rf hip
mkdir hip
cd hip

trilinos_src_dir=${TRILINOS_DIR:-${PWD}/../..}
build_type=${CMAKE_BUILD_TYPE:-release}
fortran_macro=${FORTRAN_MACRO:-FORTRAN_ONE_UNDERSCORE}
cmake_cxx_flags=${CMAKE_CXX_FLAGS}
cmake_exe_linker_flags=${CMAKE_EXE_LINKER_FLAGS}

printf "\nTRILINOS_DIR=${trilinos_src_dir}\n";
printf "CMAKE_BUILD_TYPE=${build_type}\n";
printf "CMAKE_EXE_LINKER_FLAGS=${cmake_exe_linker_flags}\n";
printf "CMAKE_CXX_FLAGS=${cmake_cxx_flags}\n";
printf "FORTRAN_MACRO=${fortran_macro}\n";
printf "\nTo change these vars, set as env vars or pass to this script like 'VAR=value run_cmake_stk'\n\n";

HIP_FLAGS="-I${MPICH_DIR}/include -L${MPICH_DIR}/lib -lmpi ${PE_MPICH_GTL_DIR_amd_gfx908} -lmpi_gtl_hsa"

cmake \
-DCMAKE_BUILD_TYPE=${build_type^^} \
-DCMAKE_C_COMPILER=hipcc \
-DCMAKE_C_FLAGS:STRING="${HIP_FLAGS}" \
-DCMAKE_CXX_COMPILER=hipcc \
-DCMAKE_CXX_STANDARD="14" \
-DCMAKE_CXX_FLAGS:STRING="-DNOT_HAVE_STK_SEACASAPREPRO_LIB -DSTK_NO_BOOST_STACKTRACE -D${fortran_macro} ${cmake_cxx_flags} ${HIP_FLAGS} -Wno-unused-command-line-argument -Werror=dangling-else" \
-DCMAKE_EXE_LINKER_FLAGS="${cmake_exe_linker_flags}" \
-DTPL_ENABLE_MPI=ON \
-DMPI_BASE_DIR:PATH=${MPICH_DIR} \
-DMPI_EXEC=mpiexec \
-DTPL_ENABLE_BLAS=ON \
-DTPL_BLAS_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_cray_mpi.a \
-DTPL_ENABLE_LAPACK=ON \
-DTPL_LAPACK_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_cray_mpi.a \
-DTPL_ENABLE_Boost=ON \
-DTPL_ENABLE_HDF5:BOOL=ON \
-DHDF5_INCLUDE_DIRS=${CRAY_HDF5_PREFIX}/include \
-DHDF5_LIBRARY_DIRS=${CRAY_HDF5_PREFIX}/lib \
-DTPL_ENABLE_Netcdf:BOOL=ON \
-DNetcdf_INCLUDE_DIRS=${CRAY_NETCDF_PREFIX}/include \
-DNetcdf_LIBRARY_DIRS=${CRAY_NETCDF_PREFIX}/lib \
-DTPL_ENABLE_Pnetcdf:BOOL=ON \
-DPnetcdf_INCLUDE_DIRS=${CRAY_PARALLEL_NETCDF_PREFIX}/include \
-DPnetcdf_LIBRARY_DIRS=${CRAY_PARALLEL_NETCDF_PREFIX}/lib \
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTrilinos_ENABLE_TESTS:BOOL=OFF \
-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
-DTrilinos_ALLOW_NO_PACKAGES:BOOL=OFF \
-DTrilinos_ASSERT_MISSING_PACKAGES=OFF \
-DTrilinos_ENABLE_Fortran:BOOL=OFF \
-DTrilinos_ENABLE_Kokkos:BOOL=ON \
-DKokkos_ENABLE_HIP:BOOL=ON \
-DKokkos_ENABLE_ROCM:BOOL=ON \
-DKokkos_ARCH_VEGA908:BOOL=ON \
-DTrilinos_ENABLE_STK:BOOL=ON \
-DTrilinos_ENABLE_STKBalance:BOOL=OFF \
-DSTK_ENABLE_TESTS:BOOL=OFF \
${trilinos_src_dir}/
