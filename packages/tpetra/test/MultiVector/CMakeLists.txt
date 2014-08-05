INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddExecutable)

IF(Kokkos_ENABLE_Cuda)
  MESSAGE("CUDA is enabled")

  # mfh 01 Jan 2014: This is boilerplate I found in
  # kokkos/core/unit_test/CMakeLists.txt.  It does
  # seem to work but I have no idea how fragile it is.
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  CUDA_COMPILE(MV_UNIT_TEST_OBJECTS_CUDA MultiVector_UnitTests_CU.cu)
  CUDA_COMPILE(MV_RCB_OBJECTS_CUDA rcb_CU.cu)
ENDIF()

IF(Tpetra_ENABLE_Thrust)
  SET(MAXNP NUM_MPI_PROCS 1)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTests
  SOURCES
    MultiVector_UnitTests ${MV_UNIT_TEST_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
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
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  )

SET(TIMING_INSTALLS "")

IF (KokkosClassic_ENABLE_TBB AND Tpetra_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
      GEMMTiming_TBB
      SOURCES GEMMTiming_TBB
      COMM mpi
      DEPLIBS kokkos kokkoslinalg kokkosnodeapi
      )
  APPEND_SET(TIMING_INSTALLS Tpetra_GEMMTiming_TBB)
ENDIF()

IF (KokkosClassic_ENABLE_ThreadPool AND Tpetra_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
      GEMMTiming_TPI
      SOURCES GEMMTiming_TPI
      COMM mpi
      DEPLIBS kokkos kokkoslinalg kokkosnodeapi
      )
  APPEND_SET(TIMING_INSTALLS Tpetra_GEMMTiming_TPI)
ENDIF()

IF (KokkosClassic_ENABLE_Thrust AND KokkosClassic_ENABLE_CUDA_DOUBLE AND Tpetra_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
      GEMMTiming_Thrust
      SOURCES GEMMTiming_Thrust
      COMM mpi
      DEPLIBS kokkos kokkoslinalg kokkosnodeapi
      )
  APPEND_SET(TIMING_INSTALLS Tpetra_GEMMTiming_Thrust)
ENDIF()

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
