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

IF (RUN_CI_SERVER)
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_DEBUG_GCC_CI_CASLDEV
      ctest_linux_continuous_mpi_debug_gcc_4.5.1_casl_dev.cmake
      TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
      )

ELSE()

  IF (RUN_REGULAR_TESTS)
    
    #
    # Run the builds of just Trilinos packages with GCC 4.5.1
    # and send results to main Trilinos CDash server to help
    # Trilinos developers address issues.
    #

    SET(DO_TRILINOS_BUILDS FALSE)

    IF (DO_TRILINOS_BUILDS)
    
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

    ENDIF() # DO_TRILINOS_BUILDS
    
    #
    # Run the builds of CASL-specific add-on packages with GCC 4.5.1 and
    # send results to the CASL CDash server on casl-dev.  Upstream
    # Trilinos packages will still be processed so that errors in
    # upstream Trilinos package library build do not show up in
    # CASL-specific add-on packages.
    #

    SET(DO_CASL_BUILDS TRUE)

    IF (DO_CASL_BUILDS)
    
      IF(INTEL12_BUILD)
        #
        # Run builds of just CASL add-on packages and LIME components 
        # with the Intel 12.0.4 compiler and submit to the CASL CDash 
        # server
        #
    
        TRILINOS_DRIVER_ADD_DASHBOARD(
          SERIAL_DEBUG_INTEL12_CASLDEV
          ctest_linux_nightly_serial_debug_icpc_12.0.4_casl_dev.cmake
          TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
          )
    
        TRILINOS_DRIVER_ADD_DASHBOARD(
          SERIAL_RELEASE_INTEL12_CASLDEV
          ctest_linux_nightly_serial_release_icpc_12.0.4_casl_dev.cmake
          TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
          )

        TRILINOS_DRIVER_ADD_DASHBOARD(
          SERIAL_RELEASE_HYBRID12_CASLDEV
          ctest_linux_nightly_serial_release_gcc_4.5.1_ifort12_casl_dev.cmake
          TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
          )

        TRILINOS_DRIVER_ADD_DASHBOARD(
          MPI_DEBUG_HYBRID12_CASLDEV
          ctest_linux_nightly_mpi_debug_gcc_4.5.1_ifort12_casl_dev.cmake
          TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
          )

      ELSE(INTEL12_BUILD)

        TRILINOS_DRIVER_ADD_DASHBOARD(
          MPI_DEBUG_HYBRID11_CASLDEV
          ctest_linux_nightly_mpi_debug_gcc_4.5.1_casl_dev.cmake
          TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
          )
    
        TRILINOS_DRIVER_ADD_DASHBOARD(
          MPI_RELEASE_HYBRID11_CASLDEV
          ctest_linux_nightly_mpi_release_gcc_4.5.1_casl_dev.cmake
          TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
          )
    
        TRILINOS_DRIVER_ADD_DASHBOARD(
          SERIAL_RELEASE_HYBRID11_CASLDEV
          ctest_linux_nightly_serial_release_gcc_4.5.1_casl_dev.cmake
          TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
          )
    
      ENDIF(INTEL12_BUILD)

    ENDIF() # DO_CASL_BUILDS

  ENDIF()
  
ENDIF()

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
