
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES shyluchol)

TRIBITS_ADD_EXECUTABLE(
  Tacho_ExampleCholSupernodes
  NOEXEPREFIX
  SOURCES Tacho_ExampleCholSupernodes.cpp
  COMM serial mpi
)
TRIBITS_ADD_EXECUTABLE(
  Tacho_ExampleSolver_double
  NOEXEPREFIX
  SOURCES Tacho_ExampleSolver_double.cpp
  COMM serial mpi
)
TRIBITS_ADD_EXECUTABLE(
  Tacho_ExampleSolver_dcomplex
  NOEXEPREFIX
  SOURCES Tacho_ExampleSolver_dcomplex.cpp
  COMM serial mpi
)
TRIBITS_ADD_EXECUTABLE(
  Tacho_ExampleDenseByBlocks
  NOEXEPREFIX
  SOURCES Tacho_ExampleDenseByBlocks.cpp
  COMM serial mpi
)


IF(Kokkos_ENABLE_OpenMP OR Kokkos_ENABLE_Serial)
  IF(TPL_ENABLE_MKL)
    TRIBITS_ADD_EXECUTABLE(
      Tacho_ExampleCholPardiso
      NOEXEPREFIX
      SOURCES Tacho_ExampleCholPardiso.cpp
      COMM serial mpi
    )
  ENDIF()
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExamplePerfTest
    NOEXEPREFIX
    SOURCES Tacho_ExamplePerfTest.cpp
    COMM serial mpi
  )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUTacho_SimpleSparseTest_File
        SOURCE_FILES test.mtx
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
