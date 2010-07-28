TRILINOS_DRIVER_SETUP()

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  MPI_OPT_DEV
#  ctest_linux_nightly_mpi_development_opt_trilinos-test.cmake
#  CTEST_INSTALLER_TYPE dev
#  TIMEOUT_MINUTES 300
#)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_OPT_DEV
  ctest_linux_nightly_serial_development_opt_gcc_4.5.0_trilinos-test.cmake
  CTEST_INSTALLER_TYPE dev
  TIMEOUT_MINUTES 300
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_DEBUG_DEV
  ctest_linux_nightly_mpi_development_debug_trilinos-test.cmake
  CTEST_INSTALLER_TYPE dev
  TIMEOUT_MINUTES 300
)

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  MPI_OPT_DEV_SHARED
#  ctest_linux_nightly_mpi_development_opt_shared_trilinos-test.cmake
#  CTEST_INSTALLER_TYPE dev
#  TIMEOUT_MINUTES 300
#)

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  SERIAL_DEBUG_DEV
#  ctest_linux_nightly_serial_development_debug_trilinos-test.cmake
#  CTEST_INSTALLER_TYPE dev
#  TIMEOUT_MINUTES 300
#)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
