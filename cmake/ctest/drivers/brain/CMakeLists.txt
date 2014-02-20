TRILINOS_DRIVER_SETUP()

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  MPI_OPT_DEV2 
#  ctest_linux_nightly_mpi_development_opt_gcc_4.6.3_brain.cmake
#  CTEST_INSTALLER_TYPE release
#  RUN_SERIAL
#  TIMEOUT_MINUTES 300
#)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG
  ctest_linux_nightly_serial_debug_gcc-4.6.3_cxx11_brain.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 360
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_MEMCHECK
  ctest_linux_nightly_serial_debug_gcc-4.1.2_memcheck_brain.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 360
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_NOEI
  ctest_linux_nightly_serial_development_debug_gcc_4.6.3_noei_brain.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 360
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  OPENMP_OPT
  ctest_linux_nightly_openmp_development_opt_brain.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 180
)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()

