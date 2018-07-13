INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUTacho_UnitTest_SparseMatrixFile
        SOURCE_FILES test_double.mtx test_dcomplex.mtx
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES shyluchol)

IF (TACHO_HAVE_KOKKOS_TASK)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tacho_TestUtil
    NOEXEPREFIX
    SOURCES Tacho_TestUtil.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
  )

  IF(Kokkos_ENABLE_Serial)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerial_double
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerial_dcomplex
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
  ENDIF()
  
  IF(Kokkos_ENABLE_OpenMP)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMP_double
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
  
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMP_dcomplex
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )  
  ENDIF()

  IF(Kokkos_ENABLE_Cuda)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCuda_double
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCuda_dcomplex
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
  ENDIF()
ENDIF()
