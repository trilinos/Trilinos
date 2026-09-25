# cmake_args()
set(Trilinos_ENABLE_ALL_PACKAGES OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~epetra_stack
set(Trilinos_ENABLE_Amesos OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Amesos2 ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +umfpack
set(Amesos2_ENABLE_UMFPACK ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Anasazi ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~epetra_stack
set(Trilinos_ENABLE_AztecOO OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Belos ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~epetra_stack
set(Trilinos_ENABLE_Epetra OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_EpetraExt OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Galeri ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~epetra_stack
set(Trilinos_ENABLE_Ifpack OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Ifpack2 ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~epetra_stack
set(Trilinos_ENABLE_Intrepid OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Intrepid2 ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - Newer than 2025.11.25
set(Intrepid2_ENABLE_Sacado OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Isorropia ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Kokkos ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_KokkosCore ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_KokkosContainers ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_MiniTensor ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~epetra_stack
set(Trilinos_ENABLE_ML OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_MOOCHO ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_MueLu ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~pamgen
set(Trilinos_ENABLE_Pamgen OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +rol
set(Trilinos_ENABLE_ROL ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_RTOp ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Rythmos OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Sacado ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_SEACAS OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Shards ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_ShyLU_NodeTacho ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_ShyLU_NodeHTS ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_STK OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Stokhos OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Teuchos ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Teko ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_ThreadPool ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Thyra ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~epetra_stack
set(Trilinos_ENABLE_ThyraEpetraExtAdapters OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Tpetra ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_TrilinosCouplings OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~epetra_stack
set(Trilinos_ENABLE_Triutils OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Xpetra ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Zoltan ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_Zoltan2 ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +shylu_fastilu
set(Trilinos_ENABLE_ShyLU_NodeFastILU ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +panzerexpreval
set(Trilinos_ENABLE_PanzerExprEval ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Trilinos_ENABLE_Gtest OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Amesos_ENABLE_SuperLU OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Ifpack_ENABLE_SuperLU OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Amesos_ENABLE_SuperLUDist OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +superlu
set(TPL_ENABLE_SuperLU ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_ENABLE_SuperLU5_API ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Amesos2_ENABLE_ShyLU_NodeTacho ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +superlu
set(Amesos2_ENABLE_SuperLU ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Amesos2_ENABLE_SuperLUDist OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(AztecOO_ENABLE_AZLU ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(ML_ENABLE_METIS OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(ML_ENABLE_ParMETIS OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(ML_ENABLE_SuperLU OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - packages
set(Ifpack2_ENABLE_ThyraTpetraAdapters ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +shylu_fastilu
set(Ifpack2_ENABLE_ShyLU_NodeFastILU ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
#TODO: make sure all variants are satisfied
# depends_on(kokkos+wraper+cuda+cuda_relocatable_device_code+cuda_constexpr) from +cuda
set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Kokkos_ENABLE_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +disable_deprecated_code
set(Tpetra_ENABLE_DEPRECATED_CODE OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Kokkos_ENABLE_DEPRECATED_CODE OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Ifpack2_ENABLE_DEPRECATED_CODE OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(MueLu_ENABLE_DEPRECATED_CODE OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Xpetra_ENABLE_DEPRECATED_CODE OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Belos_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Epetra_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Ifpack2_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Panzer_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Phalanx_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(RTop_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(STK_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Teuchos_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Thyra_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Claps_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(GlobiPack_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Tempus_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(SEACASProj_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trios_HIDE_DEPRECATED_CODE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
#TODO: make sure all variants are satisfied
# depends_on(kokkos-kernels+cuda+cublas+cusolver+cusparse+execspace_cuda+memspace_cudaspace) from +cuda
set(KokkosKernels_ENABLE_SUPERNODAL_SPTRSV OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
#TODO: make sure all variants are satisfied
# depends_on(kokkos-kernels memspace_cudauvmspace=False) from +cuda~cuda_default_uvm
# depends_on(kokkos-kernels+superlu) from +superlu"
# depends_on(superlu) from +superlu
# depends_on(umfpack) from +umfpack
# depends_on(kokkos70) from cuda_arch:=70
# depends_on(kokkos-kernels70) from cuda_arch:=70
# cmake_args() - packages
set(ROL_ENABLE_EXAMPLES OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(ROL_ENABLE_TESTS OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Sacado_ENABLE_KokkosCore ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Tacho_ENABLE_INT_INT ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Tpetra_INST_INT_LONG OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Tpetra_INST_INT_LONG_LONG ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~openmp
set(Tpetra_INST_OPENMP OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~rocm
set(Tpetra_INST_HIP OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args() - packages
set(Tpetra_INST_SERIAL ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Zoltan_ENABLE_METIS OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Zoltan_ENABLE_Scotch OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Zoltan2_ENABLE_Experimental ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +scotch
set(Zoltan2_ENABLE_Scotch ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~shared
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_FIND_SHARED_LIBS ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# build_type=Release
set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_DEBUG OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~openmp
set(Trilinos_ENABLE_OpenMP OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +cuda
set(TPL_ENABLE_CUDA ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +cuda
set(TPL_ENABLE_CUSPARSE ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +scotch 
set(TPL_ENABLE_Scotch ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +cuda
set(Tpetra_ASSUME_GPU_AWARE_MPI OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args() - config
set(CMAKE_C_COMPILER $ENV{MPICC} CACHE FILEPATH "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(CMAKE_CXX_COMPILER $ENV{MPICXX} CACHE FILEPATH "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(CMAKE_Fortran_COMPILER $ENV{MPIF90} CACHE FILEPATH "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(CMAKE_CXX_STANDARD 20 CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_COMPLEX TRUE CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_EXPLICIT_INSTANTIATION TRUE CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_LINEAR_SOLVER_FACTORY_REGISTRATION TRUE CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_WARNINGS_AS_ERRORS "" CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args - tpls
set(TPL_ENABLE_MPI ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +cuda~cuda_default_uvm
set(Tpetra_ALLOCATE_IN_SHARED_SPACE OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args()
set(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ASSERT_MISSING_PACKAGES OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(Trilinos_ENABLE_FLOAT ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~complex_float
set(Trilinos_ENABLE_COMPLEX_FLOAT OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ~deprecated_warnings
set(Trilinos_SHOW_DEPRECATED_WARNINGS OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args() - tpls - TPL INCLUDE_DIRS, LIBRARY_DIRS, and NAMES are skipped since the
# container infrastructure should mean either the right ones are found or the build will fail
set(TPL_ENABLE_BLAS ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_ENABLE_Boost ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_ENABLE_LAPACK ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_ENABLE_METIS ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_ENABLE_ParMETIS ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_ENABLE_Pthread OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_ENABLE_yaml OFF CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_ENABLE_Kokkos ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(TPL_ENABLE_KokkosKernels ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +scotch
set(TPL_ENABLE_Scotch ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# +umfpack
set(TPL_ENABLE_UMFPACK ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cmake_args() - tpls
set(TPL_ENABLE_y12m ON CACHE BOOL "Set by ${CMAKE_CURRENT_LIST_FILE}")
# generator
set(CMAKE_GENERATOR Ninja CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cxxflags='-fPIC -fno-semantic-interposition'
set(CMAKE_CXX_FLAGS -fPIC -fno-semantic-interposition CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
# cflags='-fPIC -fno-semantic-interposition'
set(CMAKE_C_FLAGS -fPIC -fno-semantic-interposition CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
# fflags='-fPIC -fno-semantic-interposition'
set(CMAKE_FORTRAN_FLAGS -fPIC -fno-semantic-interposition CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
# ldlibs='-lgfortran -lm -lz -ldl -lpthread'
set(CMAKE_C_STANDARD_LIBRARIES -lgfortran -lm -lz -ldl -lpthread CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(CMAKE_CXX_STANDARD_LIBRARIES -lgfortran -lm -lz -ldl -lpthread CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
set(CMAKE_Fortran_STANDARD_LIBRARIES -lgfortran -lm -lz -ldl -lpthread CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
# release_flags
set(CMAKE_CXX_FLAGS_RELEASE_OVERRIDE -O2 -DNDEBUG CACHE STRING "Set by ${CMAKE_CURRENT_LIST_FILE}")
