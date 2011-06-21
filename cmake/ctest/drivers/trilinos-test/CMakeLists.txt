INCLUDE("../../../TrilinosVersion.cmake")

TRILINOS_DRIVER_SETUP()

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

IF(${Trilinos_REPOSITORY_BRANCH} STREQUAL "master")
  #TRILINOS_DRIVER_ADD_DASHBOARD(
  #  MPI_OPT_DEV
  #  ctest_linux_nightly_mpi_development_opt_trilinos-test.cmake
  #  CTEST_INSTALLER_TYPE dev
  #  TIMEOUT_MINUTES 300
  #)

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT_DEV
    ctest_linux_nightly_serial_development_opt_gcc_4.5.1_trilinos-test.cmake
    CTEST_INSTALLER_TYPE dev
    RUN_SERIAL
    TIMEOUT_MINUTES 210
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_DEBUG_DEV
    ctest_linux_nightly_mpi_development_debug_trilinos-test.cmake
    CTEST_INSTALLER_TYPE dev
    RUN_SERIAL
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

  #TRILINOS_DRIVER_ADD_DASHBOARD(
  #  SERIAL_OPT_DEV2
  #  ctest_linux_nightly_serial_development_opt_pgi_10.9_trilinos-test.cmake 
  #  CTEST_INSTALLER_TYPE dev
  #  RUN_SERIAL
  #  TIMEOUT_MINUTES 450
  #)

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT_DEV1
    ctest_linux_nightly_serial_development_opt_pgi_11.1_stable_trilinos-test.cmake
    CTEST_INSTALLER_TYPE dev
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_DEV2
    ctest_linux_nightly_serial_development_debug_gcc_4.6.0_trilinos-test.cmake
    CTEST_INSTALLER_TYPE dev
    RUN_SERIAL
    TIMEOUT_MINUTES 210
  )

#Only run this build once a week
if(${DAYOFWEEK} MATCHES "Saturday")
  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT_DEV3
    ctest_linux_nightly_serial_development_opt_pgi_11.1_trilinos-test.cmake
    CTEST_INSTALLER_TYPE dev
    RUN_SERIAL
    TIMEOUT_MINUTES 600
  )
endif()

ELSE()
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_RELEASE
    ctest_linux_nightly_mpi_release_trilinos-test.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_RELEASE_SHARED
    ctest_linux_nightly_mpi_release_shared_trilinos-test.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_RELEASE
    ctest_linux_nightly_serial_release_trilinos-test.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
ENDIF()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
