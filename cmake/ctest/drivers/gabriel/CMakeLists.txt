TRILINOS_DRIVER_SETUP()

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_DEBUG_PS
  ctest_linux_nightly_mpi_debug_gabriel.cmake
  TIMEOUT_MINUTES 300
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_RELEASE_PS
  ctest_linux_nightly_serial_release_gabriel.cmake
  TIMEOUT_MINUTES 300
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_NOTRACE
  ctest_linux_nightly_serial_debug_notracing_gabriel.cmake
  TIMEOUT_MINUTES 300
  )

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
