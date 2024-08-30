KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(KOKKOSKERNELS_INST_DOUBLE)
  KOKKOSKERNELS_ADD_EXECUTABLE(
          sparse_pcg
          SOURCES KokkosSparse_pcg.cpp
  )

  KOKKOSKERNELS_ADD_EXECUTABLE(
          sparse_block_pcg
          SOURCES KokkosSparse_block_pcg.cpp
  )
ENDIF()

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_spgemm
        SOURCES KokkosSparse_spgemm.cpp
)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_spgemm_jacobi
        SOURCES KokkosSparse_spgemm_jacobi.cpp
)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/spmv)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_spadd
        SOURCES KokkosSparse_spadd.cpp
)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_spmv_struct
        SOURCES KokkosSparse_spmv_struct.cpp
)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_spmv_struct_tuning
        SOURCES KokkosSparse_spmv_struct_tuning.cpp
)


set(utilities_list)

IF (KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
    LIST(APPEND utilities_list ../PerfTestUtilities.cpp)
ENDIF ()

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_spmv
        SOURCES KokkosSparse_spmv.cpp KokkosSparse_spmv_test.cpp spmv/OpenMPSmartStatic_SPMV.cpp
        ${utilities_list}
)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_kk_spmv
        SOURCES KokkosSparse_kk_spmv.cpp
)

IF (KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
    KOKKOSKERNELS_ADD_EXECUTABLE(
            sparse_spmv_merge
            SOURCES KokkosSparse_spmv_merge.cpp
    )
ENDIF ()


KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_spmv_bsr
        SOURCES KokkosSparse_spmv_bsr.cpp
)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_sptrsv
        SOURCES KokkosSparse_sptrsv.cpp
)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_sptrsv_cholmod
        SOURCES KokkosSparse_sptrsv_cholmod.cpp
)

IF (KOKKOS_HAS_TRILINOS)
    SET (KOKKOSKERNELS_ENABLE_SPTRSV_SUPERLU OFF)
ELSE ()
    SET (KOKKOSKERNELS_ENABLE_SPTRSV_SUPERLU ON)
ENDIF ()
IF (KOKKOSKERNELS_ENABLE_SPTRSV_SUPERLU)
    # Disable this perf test with Trilinos builds to workaround
    # -Werror issues error: declaration of xyz with C language linkage
    KOKKOSKERNELS_ADD_EXECUTABLE(
            sparse_sptrsv_superlu
            SOURCES KokkosSparse_sptrsv_superlu.cpp
    )
ENDIF ()

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_sptrsv_supernode
        SOURCES KokkosSparse_sptrsv_supernode.cpp
)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_gs
        SOURCES KokkosSparse_gs.cpp
        TESTONLYLIBS kokkoskernelsperf_gtest
)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_spiluk
        SOURCES KokkosSparse_spiluk.cpp
)

KOKKOSKERNELS_ADD_EXECUTABLE(
        sparse_mdf
        SOURCES KokkosSparse_mdf.cpp
)

# For the sake of build times, don't build this CRS sorting perf test by default.
# It can be enabled if needed by setting -DKokkosKernels_ENABLE_SORT_CRS_PERFTEST=ON.
if (KokkosKernels_ENABLE_SORT_CRS_PERFTEST)
  KOKKOSKERNELS_ADD_EXECUTABLE(
    sparse_sort_crs
    SOURCES KokkosSparse_sort_crs.cpp
)
endif ()

if (KokkosKernels_ENABLE_BENCHMARK)
  KOKKOSKERNELS_ADD_BENCHMARK(
    sparse_par_ilut
    SOURCES KokkosSparse_par_ilut.cpp
  )

  # Provide -DGinkgo_DIR to cmake to enable the ginkgo test in sparse_par_ilut. Ginkgo_DIR should
  # point to the dir in the ginkgo install area that contains the GinkgoConfig.cmake file.
  # For me, this was $gingko_install_dir/lib64/cmake/Ginkgo
  if (Ginkgo_DIR)
    find_package(Ginkgo REQUIRED)

    target_compile_definitions(KokkosKernels_sparse_par_ilut PRIVATE "USE_GINKGO")
    target_link_libraries(KokkosKernels_sparse_par_ilut PRIVATE Ginkgo::ginkgo)
  endif()

  KOKKOSKERNELS_ADD_BENCHMARK(
    sparse_spmv_benchmark SOURCES KokkosSparse_spmv_benchmark.cpp
  )

  KOKKOSKERNELS_ADD_BENCHMARK(
    sparse_spmv_bsr_benchmark SOURCES KokkosSparse_spmv_bsr_benchmark.cpp
  )

  # hipcc 5.2 has an underlying clang that has the std::filesystem
  # in an experimental namespace and a different library
  if (Kokkos_CXX_COMPILER_ID STREQUAL HIPCC AND Kokkos_CXX_COMPILER_VERSION VERSION_LESS 5.3)
    target_link_libraries(KokkosKernels_sparse_spmv_bsr_benchmark PRIVATE -lstdc++fs)
  endif()
  # IntelLLVM < 2023.1.0 (and possible higher versions too) have an underlying clang that has the std::filesystem
  # in an experimental namespace and a different library
  if (Kokkos_CXX_COMPILER_ID STREQUAL IntelLLVM AND Kokkos_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 2023.1.0)
    target_link_libraries(KokkosKernels_sparse_spmv_bsr_benchmark PRIVATE -lstdc++fs)
  endif()
endif()
