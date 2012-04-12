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

SET_DEFAULT_AND_FROM_ENV(RUN_CI_SERVER FALSE)
SET_DEFAULT_AND_FROM_ENV(INTEL12_BUILD FALSE)

SET(TIMEOUT_MINUTES 300) # 5 hours

IF (RUN_REGULAR_TESTS)

  IF (INTEL12_BUILD)
    #
    # Run the builds of just Trilinos packages with Intel 12.0.4
    # and send results to main Trilinos CDash server to help
    # Trilinos developers address issues.
    #

    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_DEBUG_ICPC_TRILINOS
      ctest_linux_nightly_serial_debug_icpc_12.0.4_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )

    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_RELEASE_ICPC_TRILINOS
      ctest_linux_nightly_serial_release_icpc_12.0.4_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )

  ELSE (INTEL12_BUILD)

    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_DEBUG_GCC_TRILINOS
      ctest_linux_nightly_mpi_debug_gcc_4.6.1_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_RELEASE_GCC_TRILINOS
      ctest_linux_nightly_mpi_release_gcc_4.6.1_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_RELEASE_GCC_TRILINOS
      ctest_linux_nightly_serial_release_gcc_4.6.1_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_DEBUG_GCC_TRIKOTA
      ctest_linux_nightly_mpi_debug_gcc_4.6.1_trikota_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )

  ENDIF (INTEL12_BUILD)
  

ENDIF (RUN_REGULAR_TESTS)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
