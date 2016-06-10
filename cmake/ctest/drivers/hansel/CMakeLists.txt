INCLUDE("../../../../Version.cmake")

TRILINOS_DRIVER_SETUP()

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

IF(${Trilinos_REPOSITORY_BRANCH} STREQUAL "develop")
  MESSAGE("No testing defined for develop branch yet!")
ELSE()
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_RELEASE
    ctest_linux_nightly_mpi_release_hansel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_RELEASE_SHARED
    ctest_linux_nightly_mpi_release_shared_hansel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_RELEASE
    ctest_linux_nightly_serial_release_hansel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
ENDIF()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
