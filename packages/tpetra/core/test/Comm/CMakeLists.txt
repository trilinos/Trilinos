
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  iallreduce
  SOURCES
    iallreduce.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  extractMpiCommFromTeuchos
  SOURCES
    extractMpiCommFromTeuchos.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1-4
  COMM mpi
  )

# This test is best done with at least 2 MPI processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  isInterComm
  SOURCES
    isInterComm.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 2-4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  idot
  SOURCES
    idot.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi)

TRIBITS_ADD_TEST(
  idot
  NAME idot
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
)

TRIBITS_ADD_TEST(
  idot
  NAME idot
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
)

IF(Tpetra_INST_CUDA)
  SET(TPETRA_ENABLE_GPU ON)
ELSEIF(Tpetra_INST_HIP)
  SET(TPETRA_ENABLE_GPU ON)
ELSEIF(Tpetra_INST_SYCL)
  SET(TPETRA_ENABLE_GPU ON)
ELSE()
  SET(TPETRA_ENABLE_GPU OFF)
ENDIF()

# If we have GPUs, and we have GPU-aware MPI, also test idot
# with GPU-aware MPI explicitly turned off at runtime.
# This exercises a different code path in idot.
IF (TPETRA_ENABLE_GPU AND TPETRA_ASSUME_GPU_AWARE_MPI)
  TRIBITS_ADD_TEST(
    idot
    NAME idot_no_gpu_aware
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    ENVIRONMENT TPETRA_ASSUME_GPU_AWARE_MPI=OFF
  )

  TRIBITS_ADD_TEST(
    idot
    NAME idot_no_gpu_aware
    COMM mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 4
    ENVIRONMENT TPETRA_ASSUME_GPU_AWARE_MPI=OFF
  )
ENDIF ()

ASSERT_DEFINED (Tpetra_INST_CUDA)
ASSERT_DEFINED (Tpetra_INST_HIP)
ASSERT_DEFINED (Tpetra_INST_SYCL)

# mfh 21 Aug 2017: The "GPU-aware MPI" test still builds and runs if
# TPETRA_ASSUME_GPU_AWARE_MPI is OFF, but it is a trivial run (doing
# nothing and passing) unless the TPETRA_ASSUME_GPU_AWARE_MPI
# environment variable is set (to some "true" value, like "1", "ON",
# "YES", or "TRUE").  This makes the test useful for cross-compilation
# and other situations where the MPI implementation's GPU-awareness
# cannot be known until run time.

IF (TPETRA_ENABLE_GPU AND TPETRA_ASSUME_GPU_AWARE_MPI)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    gpuAwareMpi
    SOURCES
      gpuAwareMpi.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    NUM_MPI_PROCS 2
    COMM mpi
    )
ENDIF ()

IF ( (TPETRA_INST_CUDA OR TPETRA_INST_HIP)  AND TPETRA_ASSUME_GPU_AWARE_MPI)
  # Does not work with SYCL
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    gpuMpi_CanaryTest
    SOURCES
      gpuMpi_CanaryTest.cpp
    NUM_MPI_PROCS 2
    COMM mpi
    )
ENDIF ()
