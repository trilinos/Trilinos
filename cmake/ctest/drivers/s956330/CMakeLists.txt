INCLUDE("../../../../Version.cmake")

TRILINOS_DRIVER_SETUP()

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Get the day of the week
#execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

IF(${Trilinos_REPOSITORY_BRANCH} STREQUAL "master")
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_DEBUG_DEV
    ctest_linux_nightly_mpi_development_debug_gcc_4.8.2_s956330.cmake
    CTEST_INSTALLER_TYPE release
    TIMEOUT_MINUTES 480
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT_DEV
    ctest_linux_nightly_serial_development_opt_gcc_4.8.2_s956330.cmake
    CTEST_INSTALLER_TYPE release
    TIMEOUT_MINUTES 480
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_OPT_DEV_SHARED
    ctest_linux_nightly_mpi_development_opt_shared_gcc_4.8.2_s956330.cmake
    CTEST_INSTALLER_TYPE release
    TIMEOUT_MINUTES 480
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT_DEV_CLANG
    ctest_linux_nightly_serial_development_opt_clang_3.4_s956330.cmake
    CTEST_INSTALLER_TYPE release
    TIMEOUT_MINUTES 480
  )

ELSE()
#release tests would go here
ENDIF()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
