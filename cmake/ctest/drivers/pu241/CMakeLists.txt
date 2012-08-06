# pu241 driver

TRILINOS_DRIVER_SETUP()

# ToDo: I think we can remove these two variables
set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

if(${DAYOFWEEK} MATCHES "Saturday")
  set(RUN_REGULAR_TESTS ON)
elseif(${DAYOFWEEK} MATCHES "Sunday")
  set(RUN_REGULAR_TESTS ON)
else()
  set(RUN_REGULAR_TESTS ON)
endif()

SET(TIMEOUT_MINUTES 300) # 5 hours

IF (RUN_REGULAR_TESTS)

  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_DEBUG_GCC
    ctest_linux_nightly_mpi_debug_gcc_4.6.1_trilinos.cmake
    TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_RELEASE_GCC
    ctest_linux_nightly_mpi_release_gcc_4.6.1_trilinos.cmake
    TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_RELEASE_GCC
    ctest_linux_nightly_serial_release_gcc_4.6.1_trilinos.cmake
    TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
    )

ENDIF (RUN_REGULAR_TESTS)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
