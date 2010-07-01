TRILINOS_DRIVER_SETUP()

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_OPT_DEV
  ctest_linux_nightly_mpi_development_opt_trilinos-test2.cmake
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_DEBUG_DEV
  ctest_linux_nightly_mpi_development_debug_trilinos-test2.cmake
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_OPT_DEV_SHARED
  ctest_linux_nightly_mpi_development_opt_shared_trilinos-test2.cmake
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_DEV
  ctest_linux_nightly_serial_development_debug_trilinos-test2.cmake
)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
