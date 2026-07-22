kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(KOKKOSKERNELS_INST_DOUBLE)
  kokkoskernels_add_executable(sparse_pcg       SOURCES KokkosSparse_pcg.cpp)
  kokkoskernels_add_executable(sparse_block_pcg SOURCES KokkosSparse_block_pcg.cpp)
endif()

kokkoskernels_add_executable(sparse_spgemm        SOURCES KokkosSparse_spgemm.cpp)
kokkoskernels_add_executable(sparse_spgemm_jacobi SOURCES KokkosSparse_spgemm_jacobi.cpp)

kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/spmv)

kokkoskernels_add_executable(sparse_spadd              SOURCES KokkosSparse_spadd.cpp)
kokkoskernels_add_executable(sparse_spmv_struct        SOURCES KokkosSparse_spmv_struct.cpp)
kokkoskernels_add_executable(sparse_spmv_struct_tuning SOURCES KokkosSparse_spmv_struct_tuning.cpp)

set(utilities_list)

if(KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
  list(APPEND utilities_list ../PerfTestUtilities.cpp)
endif()

kokkoskernels_add_executable(sparse_spmv
  SOURCES KokkosSparse_spmv.cpp KokkosSparse_spmv_test.cpp spmv/OpenMPSmartStatic_SPMV.cpp ${utilities_list})

kokkoskernels_add_executable(sparse_kk_spmv SOURCES KokkosSparse_kk_spmv.cpp)

if(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
  kokkoskernels_add_executable(sparse_spmv_merge SOURCES KokkosSparse_spmv_merge.cpp)
endif()

kokkoskernels_add_executable(sparse_spmv_bsr       SOURCES KokkosSparse_spmv_bsr.cpp)
kokkoskernels_add_executable(sparse_sptrsv         SOURCES KokkosSparse_sptrsv.cpp)
kokkoskernels_add_executable(sparse_sptrsv_cholmod SOURCES KokkosSparse_sptrsv_cholmod.cpp)

if(KOKKOS_HAS_TRILINOS)
  set(KOKKOSKERNELS_ENABLE_SPTRSV_SUPERLU OFF)
else()
  set(KOKKOSKERNELS_ENABLE_SPTRSV_SUPERLU ON)
endif()
if(KOKKOSKERNELS_ENABLE_SPTRSV_SUPERLU)
  # Disable this perf test with Trilinos builds to workaround
  # -Werror issues error: declaration of xyz with C language linkage
  kokkoskernels_add_executable(sparse_sptrsv_superlu SOURCES KokkosSparse_sptrsv_superlu.cpp)
endif()

kokkoskernels_add_executable(sparse_sptrsv_supernode SOURCES KokkosSparse_sptrsv_supernode.cpp)
kokkoskernels_add_executable(sparse_gs               SOURCES KokkosSparse_gs.cpp)
kokkoskernels_add_executable(sparse_spiluk           SOURCES KokkosSparse_spiluk.cpp)
kokkoskernels_add_executable(sparse_mdf              SOURCES KokkosSparse_mdf.cpp)

# For the sake of build times, don't build this CRS sorting perf test by default.
# It can be enabled if needed by setting -DKokkosKernels_ENABLE_SORT_CRS_PERFTEST=ON.
if(KokkosKernels_ENABLE_SORT_CRS_PERFTEST)
  kokkoskernels_add_executable(sparse_sort_crs SOURCES KokkosSparse_sort_crs.cpp)
endif()
