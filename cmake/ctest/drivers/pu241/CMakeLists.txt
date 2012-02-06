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
  
  #
  # Run the builds of just Trilinos packages with GCC 4.5.1
  # and send results to main Trilinos CDash server to help
  # Trilinos developers address issues.
  #

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
  ELSE(INTEL12_BUILD)

    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_DEBUG_GCC_TRILINOS
      ctest_linux_nightly_mpi_debug_gcc_4.5.1_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_RELEASE_GCC_TRILINOS
      ctest_linux_nightly_mpi_release_gcc_4.5.1_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_RELEASE_GCC_TRILINOS
      ctest_linux_nightly_serial_release_gcc_4.5.1_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_DEBUG_GCC_TRIKOTA
      ctest_linux_nightly_mpi_debug_gcc_4.5.1_trikota_trilinos.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )
  ENDIF(INTEL12_BUILD)
  

ENDIF(RUN_REGULAR_TESTS)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()


# NOTE: Above, we run the Trilinos builds first so that if there are Trilinos
# errors, they will get submitted to the main Trilinos CDash site and send out
# emails before the CASL builds are done that are sent to the CASL CDash site.
# This will allow Trilinos developers to see the problems on a CDash site they
# can get to before they get emails about a CDash site they can't get to.  We
# are also running the CASL add-on packages last so that when we get to the
# expensive tests in VRIPSS, these can be left to run for some time and since
# they don't take a lot of processes, can overlap with CI testing that starts
# up in the morning.
