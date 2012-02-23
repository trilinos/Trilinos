INCLUDE("../../../../Version.cmake")

TRILINOS_DRIVER_SETUP()

IF(${Trilinos_REPOSITORY_BRANCH} STREQUAL "master")
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_OPT
    ctest_mac_nightly_mpi_opt_rouson-mp.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG
    ctest_mac_nightly_serial_debug_rouson-mp.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 330
  )
ELSE()
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_OPT_RELEASE
    ctest_mac_nightly_mpi_opt_release_rouson-mp.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_RELEASE
    ctest_mac_nightly_serial_debug_release_rouson-mp.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 330
  )
ENDIF()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
