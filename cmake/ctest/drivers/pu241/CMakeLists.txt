# pu241 driver

TRILINOS_DRIVER_SETUP()

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

IF (RUN_CI_SERVER)
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_DEBUG_GCC_CI_CASLDEV
      ctest_linux_continuous_mpi_debug_gcc_4.5.1_casl_dev.cmake
      TIMEOUT_MINUTES 240
      )

ELSE()

  IF (RUN_REGULAR_TESTS)
  
    #
    # Run the builds of just Trilinos packages with GCC 4.5.1
    # and send results to main Trilinos CDash server to help
    # Trilinos developers address issues.
    #
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_DEBUG_GCC_TRILINOS
      ctest_linux_nightly_mpi_debug_gcc_4.5.1_trilinos.cmake
      TIMEOUT_MINUTES 240
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_RELEASE_GCC_TRILINOS
      ctest_linux_nightly_serial_release_gcc_4.5.1_trilinos.cmake
      TIMEOUT_MINUTES 240
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_DEBUG_GCC_TRIKOTA
      ctest_linux_nightly_mpi_debug_gcc_4.5.1_trikota_trilinos.cmake
      TIMEOUT_MINUTES 240
      )
  
    #
    # Run the builds of CASL-specific add-on packages with GCC 4.5.1 and
    # send results to the CASL CDash server on casl-dev.  Upstream
    # Trilinos packages will still be processed so that errors in
    # upstream Trilinos package library build do not show up in
    # CASL-specific add-on packages.
    #
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      MPI_DEBUG_GCC_CASLDEV
      ctest_linux_nightly_mpi_debug_gcc_4.5.1_casl_dev.cmake
      TIMEOUT_MINUTES 240
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_RELEASE_GCC_CASLDEV
      ctest_linux_nightly_serial_release_gcc_4.5.1_casl_dev.cmake
      TIMEOUT_MINUTES 240
      )
  
    #
    # Run builds of just CASL add-on packages and LIME components 
    # with the Intel 12.0.4 compiler and submit to the CASL CDash 
    # server
    #
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_DEBUG_ICPC_CASLDEV
      ctest_linux_nightly_serial_debug_icpc_12.0.4_casl_dev.cmake
      TIMEOUT_MINUTES 240
      )
  
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_RELEASE_ICPC_CASLDEV
      ctest_linux_nightly_serial_release_icpc_12.0.4_casl_dev.cmake
      TIMEOUT_MINUTES 240
      )
  
    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_RELEASE_ICPC_CASLDEV_NATIVEBUILDS
      ctest_linux_nightly_serial_release_icpc_12.0.4_casl_dev_native.cmake
      TIMEOUT_MINUTES 240
      )

  ENDIF()
  
ENDIF()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
