INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageAddExecutable)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTests
  SOURCES
    MultiVector_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  )

SET(TIMING_INSTALLS "")

IF (Kokkos_ENABLE_TBB AND Tpetra_ENABLE_MPI)
  PACKAGE_ADD_EXECUTABLE(
      GEMMTiming_TBB
      SOURCES GEMMTiming_TBB
      COMM mpi
      DEPLIBS kokkos kokkoslinalg kokkosnodeapi
      )
  APPEND_SET(TIMING_INSTALLS Tpetra_GEMMTiming_TBB)
ENDIF()

IF (Kokkos_ENABLE_ThreadPool AND Tpetra_ENABLE_MPI)
  PACKAGE_ADD_EXECUTABLE(
      GEMMTiming_TPI
      SOURCES GEMMTiming_TPI
      COMM mpi
      DEPLIBS kokkos kokkoslinalg kokkosnodeapi
      )
  APPEND_SET(TIMING_INSTALLS Tpetra_GEMMTiming_TPI)
ENDIF()

IF (Kokkos_ENABLE_Thrust AND Kokkos_ENABLE_CUDA_DOUBLE AND Tpetra_ENABLE_MPI)
  PACKAGE_ADD_EXECUTABLE(
      GEMMTiming_Thrust
      SOURCES GEMMTiming_Thrust
      COMM mpi
      DEPLIBS kokkos kokkoslinalg kokkosnodeapi
      )
  APPEND_SET(TIMING_INSTALLS Tpetra_GEMMTiming_Thrust)
ENDIF()

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
