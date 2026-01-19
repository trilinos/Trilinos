
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES tacho)

TRIBITS_ADD_EXECUTABLE(
  Tacho_ExampleCombineDataFileToMatrixMarketFile.x
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES Tacho_ExampleCombineDataFileToMatrixMarketFile.cpp
  COMM serial mpi
)

#
# Driver
#
TRIBITS_ADD_EXECUTABLE(
  Tacho_ExampleDriverDouble.x
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES Tacho_ExampleDriver_double.cpp
  CATEGORIES BASIC PERFORMANCE
  COMM serial mpi
)
TRIBITS_ADD_EXECUTABLE(
  Tacho_ExampleDriverDoubleComplex.x
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES Tacho_ExampleDriver_dcomplex.cpp
  COMM serial mpi
)

#
# External Interface
# 
TRIBITS_ADD_EXECUTABLE(
  Tacho_ExampleExternalInterface.x
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES Tacho_ExampleExternalInterface.cpp
  COMM serial mpi
)

#
# NVIDIA cuSolver
# 
IF(Kokkos_ENABLE_CUDA) 
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholCuSolver.x
    NOEXESUFFIX
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholCuSolver.cpp
    COMM serial mpi
  )
ENDIF()

#
# Intel MKL Pardiso 
#
IF(Kokkos_ENABLE_OPENMP OR Kokkos_ENABLE_SERIAL)
  IF(TPL_ENABLE_MKL)
    TRIBITS_ADD_EXECUTABLE(
      Tacho_ExampleCholPardiso.x
      NOEXESUFFIX
      NOEXEPREFIX
      SOURCES Tacho_ExampleCholPardiso.cpp
      COMM serial mpi
    )
  ENDIF()
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tacho_SimpleSparseTest_File
        SOURCE_FILES test2.mtx test.mtx graph.dat weight.dat graph_test.mtx 
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        NOEXEPREFIX
        CATEGORIES BASIC PERFORMANCE
)


#
# Performance Test
#
TRIBITS_ADD_TEST(
  Tacho_ExampleDriverDouble.x
  ARGS "--niters=100 --nfacts=1 --nsolves=5 --default --method=chol --variant=2"
  NAME "Tacho_ExampleDriverDouble_Chol_v2"
  NOEXESUFFIX
  NOEXEPREFIX
  NUM_MPI_PROCS 1
  CATEGORIES BASIC PERFORMANCE
)

TRIBITS_ADD_TEST(
  Tacho_ExampleDriverDouble.x
  ARGS "--niters=100 --nfacts=1 --nsolves=5 --default --method=chol --variant=3"
  NAME "Tacho_ExampleDriverDouble_Chol_v3"
  NOEXESUFFIX
  NOEXEPREFIX
  NUM_MPI_PROCS 1
  CATEGORIES BASIC PERFORMANCE
)

TRIBITS_ADD_TEST(
  Tacho_ExampleDriverDouble.x
  ARGS "--niters=100 --nfacts=1 --nsolves=5 --default --method=lu --variant=2"
  NAME "Tacho_ExampleDriverDouble_LU_v2"
  NOEXESUFFIX
  NOEXEPREFIX
  NUM_MPI_PROCS 1
  CATEGORIES BASIC PERFORMANCE
)

TRIBITS_ADD_TEST(
  Tacho_ExampleDriverDouble.x
  ARGS "--niters=100 --nfacts=1 --nsolves=5 --default --method=lu --variant=3"
  NAME "Tacho_ExampleDriverDouble_LU_v3"
  NOEXESUFFIX
  NOEXEPREFIX
  NUM_MPI_PROCS 1
  CATEGORIES BASIC PERFORMANCE
)
