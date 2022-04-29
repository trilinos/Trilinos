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
        sparse_spmv_blockcrs
        SOURCES KokkosSparse_spmv_blockcrs.cpp
)

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

IF (NOT ${KOKKOS_HAS_TRILINOS})
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
