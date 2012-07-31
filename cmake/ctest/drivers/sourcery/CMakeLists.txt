INCLUDE("../../../../Version.cmake")

TRILINOS_DRIVER_SETUP()

IF(${Trilinos_REPOSITORY_BRANCH} STREQUAL "master")
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_OPT_GCC
    ctest_mac_nightly_mpi_opt_gcc_sourcery.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_OPT_NAG
    ctest_mac_nightly_mpi_opt_nag_sourcery.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_GCC
    ctest_mac_nightly_serial_debug_gcc_sourcery.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 330
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_NAG
    ctest_mac_nightly_serial_debug_nag_sourcery.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 330
  )

ELSE()
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_OPT_RELEASE_NAG
    ctest_mac_nightly_mpi_opt_release_nag_sourcery.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_RELEASE_NAG
    ctest_mac_nightly_serial_debug_release_nag_sourcery.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 330
  )
ENDIF()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
