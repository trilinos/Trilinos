# u233 driver

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

if (RUN_REGULAR_TESTS)

  #
  # Run builds of just Trilinos packages with the Intel 11.???
  # compiler and submit to the main Trilinos CDash sever so that
  # Trilinos developers can help maintain their packages on this
  # platform.
  #

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_ICPC
    ctest_linux_nightly_serial_debug_icpc_u233_trilinos.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_RELEASE_ICPC
    ctest_linux_nightly_serial_release_icpc_u233_trilinos.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
    )

endif()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
