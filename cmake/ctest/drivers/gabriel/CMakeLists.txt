TRILINOS_DRIVER_SETUP()

# Set a timeout of 7 hours for an individual build run in parallel with other
# builds
SET(TIMEOUT_MINUTES 420)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_DEBUG_PS
  ctest_linux_nightly_mpi_debug_gabriel.cmake
  CTEST_INSTALLER_TYPE dev
  TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_RELEASE_PS
  ctest_linux_nightly_serial_release_gabriel.cmake
  CTEST_INSTALLER_TYPE dev
  TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_NOTRACE
  ctest_linux_nightly_serial_debug_notracing_gabriel.cmake
  CTEST_INSTALLER_TYPE dev
  TIMEOUT_MINUTES 150
  )

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
