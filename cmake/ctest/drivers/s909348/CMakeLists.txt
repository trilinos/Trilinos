INCLUDE("../../../TrilinosVersion.cmake")

TRILINOS_DRIVER_SETUP()

IF(${Trilinos_REPOSITORY_BRANCH} STREQUAL "master")
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_OPT
    ctest_mac_nightly_mpi_opt_s909348.cmake
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG
    ctest_mac_nightly_serial_debug_s909348.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
ELSE()
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_RELEASE
    ctest_mac_nightly_mpi_release_s909348.cmake
    TIMEOUT_MINUTES 300
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_RELEASE
    ctest_mac_nightly_serial_release_s909348.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
  
ENDIF()


TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
