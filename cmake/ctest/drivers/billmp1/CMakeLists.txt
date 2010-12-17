TRILINOS_DRIVER_SETUP()

# Set a timeout of 6 hours for an individual build run in parallel with other
# builds
SET(TIMEOUT_MINUTES 120)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_DEBUG_PS
  ctest_osx_nightly_mpi_debug_billmp1.cmake
  CTEST_INSTALLER_TYPE dev
  TIMEOUT_MINUTES ${TIMEOUT_MINUTES}
  )

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
