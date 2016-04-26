
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES shyluchol)

# Dense Matrix Objects
# ====================

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseMatrixBase_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseMatrixBase_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseMatrixView_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseMatrixView_Pthreads.cpp
    COMM serial mpi
  )
ENDIF()

IF(Kokkos_ENABLE_Cuda)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseMatrixBase_Cuda
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseMatrixBase_Cuda.cpp
    COMM serial mpi
  )
ENDIF()

# Dense Linear Algebra
# ====================

IF(Kokkos_ENABLE_OpenMP)
IF(TPL_ENABLE_MKL)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseGemmMKL_OpenMP
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseGemmMKL_OpenMP.cpp
    COMM serial mpi
  )
ENDIF()
ENDIF()

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseGemmByBlocks_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseGemmByBlocks_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseHerkByBlocks_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseHerkByBlocks_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseTrsmByBlocks_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseTrsmByBlocks_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseCholByBlocks_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseCholByBlocks_Pthreads.cpp
    COMM serial mpi
  )
ENDIF()

# Crs Matrix Objects
# ==================

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCrsMatrixBase_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleCrsMatrixBase_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCrsMatrixView_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleCrsMatrixView_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleMatrixMarket_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleMatrixMarket_Pthreads.cpp
    COMM serial mpi
  )
ENDIF()

# Graph Tools and Symbolic Factorization
# ======================================

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleGraphTools_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleGraphTools_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleSymbolicFactorization_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleSymbolicFactorization_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleHierMatrix_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleHierMatrix_Pthreads.cpp
    COMM serial mpi
  )
ENDIF()


# Numeric Factorization
# =====================

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholUnblocked_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholUnblocked_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholByBlocks_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholByBlocks_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholSuperNodesByBlocks_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholSuperNodesByBlocks_Pthreads.cpp
    COMM serial mpi
  )
ENDIF()

# IF(Kokkos_ENABLE_QTHREAD) 
#   TRIBITS_ADD_EXECUTABLE(
#     Tacho_ExampleCholSuperNodesByBlocks_Qthreads
#     NOEXEPREFIX
#     SOURCES Tacho_ExampleCholSuperNodesByBlocks_Qthreads.cpp
#     COMM serial mpi
#   )
# ENDIF()


# Solver Interface
# ================

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleSolver_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleSolver_Pthreads.cpp
    COMM serial mpi
  )
ENDIF()


TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUTacho_SimpleSparseTest_File
        SOURCE_FILES test.mtx
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
