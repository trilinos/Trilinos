KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_ADD_EXECUTABLE(
  sparse_pcg
  SOURCES KokkosSparse_pcg.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  sparse_block_pcg
  SOURCES KokkosSparse_block_pcg.cpp
  )

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

KOKKOSKERNELS_ADD_EXECUTABLE(
  sparse_spmv
  SOURCES KokkosSparse_spmv.cpp
  )

IF(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
  KOKKOSKERNELS_ADD_EXECUTABLE(
    sparse_spmv_merge
    SOURCES KokkosSparse_spmv_merge.cpp
    )
ENDIF()

KOKKOSKERNELS_ADD_EXECUTABLE(
  sparse_sptrsv
  SOURCES KokkosSparse_sptrsv.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  sparse_sptrsv_cholmod
  SOURCES KokkosSparse_sptrsv_cholmod.cpp 
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  sparse_sptrsv_superlu
  SOURCES KokkosSparse_sptrsv_superlu.cpp 
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  sparse_sptrsv_supernode
  SOURCES KokkosSparse_sptrsv_supernode.cpp 
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  sparse_gs
  SOURCES KokkosSparse_gs.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  sparse_spiluk
  SOURCES KokkosSparse_spiluk.cpp 
  )
