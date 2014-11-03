
IF(Tpetra_ENABLE_Thrust)
  SET(MAXNP NUM_MPI_PROCS 1)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTests
  SOURCES
    MultiVector_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  )


# This test comes from zoltan2/test/temp/multivectorTest.cpp.  It
# mimics the recursive bisection algorithm in Zoltan2, and times data
# migration.
#
# The test itself does not refer to MPI.  Hoewver, it only makes sense
# to run this test in an MPI build.  Otherwise, no data migration is
# necessary.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  rcb
  SOURCES rcb.cpp ${MV_RCB_OBJECTS_CUDA}
  COMM mpi
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  )

SET(TIMING_INSTALLS "")

IF (KokkosClassic_ENABLE_TBB AND Tpetra_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
      GEMMTiming_TBB
      SOURCES GEMMTiming_TBB
      COMM mpi
      )
  APPEND_SET(TIMING_INSTALLS Tpetra_GEMMTiming_TBB)
ENDIF()

IF (KokkosClassic_ENABLE_ThreadPool AND Tpetra_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
      GEMMTiming_TPI
      SOURCES GEMMTiming_TPI
      COMM mpi
      )
  APPEND_SET(TIMING_INSTALLS Tpetra_GEMMTiming_TPI)
ENDIF()

IF (KokkosClassic_ENABLE_Thrust AND KokkosClassic_ENABLE_CUDA_DOUBLE AND Tpetra_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
      GEMMTiming_Thrust
      SOURCES GEMMTiming_Thrust
      COMM mpi
      )
  APPEND_SET(TIMING_INSTALLS Tpetra_GEMMTiming_Thrust)
ENDIF()

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
