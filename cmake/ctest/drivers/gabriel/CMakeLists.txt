TRILINOS_DRIVER_SETUP()

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_DEBUG_PS
  ctest_linux_nightly_mpi_debug_gabriel.cmake
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_RELEASE_PS
  ctest_linux_nightly_serial_release_gabriel.cmake
  )

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  SERIAL_DEBUG_BOOST_TRACE
#  ctest_linux_nightly_serial_debug_boost_tracing_gabriel.cmake
#  )

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
