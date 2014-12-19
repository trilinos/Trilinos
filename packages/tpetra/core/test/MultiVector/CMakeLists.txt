
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

IF (TpetraClassic_ENABLE_TBB AND ${PACKAGE_NAME}_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
      GEMMTiming_TBB
      SOURCES GEMMTiming_TBB
      COMM mpi
      )
  APPEND_SET(TIMING_INSTALLS ${PACKAGE_NAME}_GEMMTiming_TBB)
ENDIF()

IF (TpetraClassic_ENABLE_ThreadPool AND ${PACKAGE_NAME}_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
      GEMMTiming_TPI
      SOURCES GEMMTiming_TPI
      COMM mpi
      )
  APPEND_SET(TIMING_INSTALLS ${PACKAGE_NAME}_GEMMTiming_TPI)
ENDIF()

IF (TpetraCore_ENABLE_TpetraClassic AND TpetraClassic_ENABLE_Thrust AND TpetraClassic_ENABLE_CUDA_DOUBLE AND ${PACKAGE_NAME}_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
      GEMMTiming_Thrust
      SOURCES GEMMTiming_Thrust
      COMM mpi
      )
  APPEND_SET(TIMING_INSTALLS ${PACKAGE_NAME}_GEMMTiming_Thrust)
ENDIF()

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
