
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

ADD_SUBDIRECTORIES(gtest_unit)

#MD: test is moved to gtest folder.
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  CrsMatrixUnitTest
#  SOURCES
#    TestCrsMatrix.cpp
#    ${TEUCHOS_STD_UNIT_TEST_MAIN}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  FAIL_REGULAR_EXPRESSION "  FAILED  "
#  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ArithTraits
  SOURCES test-arith-traits.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  blas1_MV_unit
  SOURCES
    blas1_MV
    KokkosKernels_Blas1_MV_UnitTest_Cuda.cpp
    KokkosKernels_Blas1_MV_UnitTest_OpenMP.cpp
    KokkosKernels_Blas1_MV_UnitTest_Serial.cpp
    KokkosKernels_Blas1_MV_UnitTests.hpp
    KokkosKernels_Blas1_MV_UnitTest_Threads.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

IF (Kokkos_ENABLE_Cuda)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    blas2_MV_Cuda
    SOURCES
      blas2_MV_Cuda
      KokkosKernels_Blas2_MV_GEMV.hpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF ()

IF (Kokkos_ENABLE_OpenMP)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    blas2_MV_OpenMP
    SOURCES
      blas2_MV_OpenMP
      KokkosKernels_Blas2_MV_GEMV.hpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF ()

IF (Kokkos_ENABLE_Serial)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    blas2_MV_Serial
    SOURCES
      blas2_MV_Serial
      KokkosKernels_Blas2_MV_GEMV.hpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF ()

IF (Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    blas2_MV_Threads
    SOURCES
      blas2_MV_Threads
      KokkosKernels_Blas2_MV_GEMV.hpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF ()

ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_MKL)
IF (${PACKAGE_NAME}_ENABLE_MKL)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    blas2_MV_unit_MKL
    SOURCES
      blas2_MV_MKL
      KokkosKernels_Blas2_MV_GEMV_MKL.hpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF ()

ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_quadmath)
IF (${PACKAGE_NAME}_ENABLE_quadmath)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    float128
    SOURCES float128.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF ()

# This test needs to be able to control when Kokkos gets initialized,
# so it can't use the Teuchos unit test "main" .cpp file.
#MD: Test is moved to gtest directory.
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  findRelOffset
#  SOURCES findRelOffset.cpp
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  )

# This test needs to be able to control when Kokkos gets initialized,
# so it can't use the Teuchos unit test "main" .cpp file.
#MD:test is moved to gtest
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  replaceSumInto
#  SOURCES replaceSumInto.cpp
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  )

# This test needs to be able to control when Kokkos gets initialized,
# so it can't use the Teuchos unit test "main" .cpp file.
#MD: test is moved to gtest
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  replaceSumIntoLonger
#  SOURCES replaceSumIntoLonger.cpp
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  )

#This test needs to be able to control when Kokkos gets initialized,
# so it can't use the Teuchos unit test "main" .cpp file.
#MD: test is moved to gtest
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  copyIntegers
#  SOURCES copyIntegers.cpp
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  )

# This test needs to be able to control when Kokkos gets initialized,
# so it can't use the Teuchos unit test "main" .cpp file.
#MD: test is moved to gtest
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  MKL_sparseMatrixHandle
#  SOURCES MKL_sparseMatrixHandle.cpp
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  )
 