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

  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   Tacho_TestUtilDeprecatedTaskScheduler 
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
  #   NOEXEPREFIX
  #   SOURCES Tacho_TestUtil.cpp
  #   ARGS PrintItAll
  #   NUM_MPI_PROCS 1
  #   FAIL_REGULAR_EXPRESSION "  FAILED  "
  # )
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   Tacho_TestUtilDeprecatedTaskSchedulerMultiple
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
  #   NOEXEPREFIX
  #   SOURCES Tacho_TestUtil.cpp
  #   ARGS PrintItAll
  #   NUM_MPI_PROCS 1
  #   FAIL_REGULAR_EXPRESSION "  FAILED  "
  # )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tacho_TestUtilTaskScheduler
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_TestUtil.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
  )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tacho_TestUtilTaskSchedulerMultiple
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
    NOEXEPREFIX
    SOURCES Tacho_TestUtil.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
  )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tacho_TestUtilChaseLevTaskScheduler
    TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_TestUtil.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
  )

  IF(Kokkos_ENABLE_Serial)
    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestSerialDoubleDeprecatedTaskScheduler
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestSerial_double.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestSerialDoubleDeprecatedTaskSchedulerMultiple
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestSerial_double.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerialDoubleTaskScheduler
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerialDoubleTaskSchedulerMultiple
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerialDoubleChaseLevTaskScheduler
      TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )

    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestSerialDoubleComplexDeprecatedTaskScheduler
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestSerial_dcomplex.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestSerialDoubleComplexDeprecatedTaskSchedulerMultiple
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestSerial_dcomplex.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerialDoubleComplexTaskScheduler
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerialDoubleComplexTaskSchedulerMultiple
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestSerialDoubleComplexChaseLevTaskScheduler
      TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestSerial_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )

  ENDIF()
  
  IF(Kokkos_ENABLE_OpenMP)
    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestOpenMPDoubleDeprecatedTaskScheduler
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestOpenMP_double.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestOpenMPDoubleDeprecatedTaskSchedulerMultiple
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestOpenMP_double.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMPDoubleTaskScheduler
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMPDoubleTaskSchedulerMultiple
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMPDoubleChaseLevTaskScheduler
      TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )

    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestOpenMPDoubleComplexDeprecatedTaskScheduler
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestOpenMP_dcomplex.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestOpenMPDoubleComplexDeprecatedTaskSchedulerMultiple
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestOpenMP_dcomplex.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMPDoubleComplexTaskScheduler
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMPDoubleComplexTaskSchedulerMultiple
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestOpenMPDoubleComplexChaseLevTaskScheduler
      TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestOpenMP_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
  ENDIF()

  IF(Kokkos_ENABLE_Cuda)
    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestCudaDoubleDeprecatedTaskScheduler
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestCuda_double.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestCudaDoubleDeprecatedTaskSchedulerMultiple
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestCuda_double.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCudaDoubleTaskScheduler
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCudaDoubleTaskSchedulerMultiple
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCudaDoubleChaseLevTaskScheduler
      TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_double.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )

    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestCudaDoubleComplexDeprecatedTaskScheduler
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestCuda_dcomplex.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    # TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #   Tacho_TestCudaDoubleComplexDeprecatedTaskSchedulerMultiple
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
    #   NOEXEPREFIX
    #   SOURCES Tacho_TestCuda_dcomplex.cpp
    #   ARGS PrintItAll
    #   NUM_MPI_PROCS 1
    #   FAIL_REGULAR_EXPRESSION "  FAILED  "
    # )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCudaDoubleComplexTaskScheduler
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCudaDoubleComplexTaskSchedulerMultiple
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tacho_TestCudaDoubleComplexChaseLevTaskScheduler
      TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_TestCuda_dcomplex.cpp
      ARGS PrintItAll
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
    )
  ENDIF()
ENDIF()
