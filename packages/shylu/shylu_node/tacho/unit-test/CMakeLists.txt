INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tacho_UnitTest_SparseMatrixFile
        SOURCE_FILES test_double.mtx test_dcomplex.mtx
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES tacho)

IF (TACHO_HAVE_KOKKOS_TASK)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tacho_TestUtil
    NOEXEPREFIX
    SOURCES Tacho_TestUtil.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
  )
  IF(Kokkos_ENABLE_SERIAL)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerialDouble
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerialDoubleComplex
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
  ENDIF()
  
  IF(Kokkos_ENABLE_OPENMP)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMPDouble
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMPDoubleComplex
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
  ENDIF()

  IF(Kokkos_ENABLE_CUDA)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCudaDouble
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCudaDoubleComplex
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
  ENDIF()
ENDIF()
