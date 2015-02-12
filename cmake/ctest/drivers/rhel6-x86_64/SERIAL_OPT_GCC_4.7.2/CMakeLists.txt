INCLUDE("../../../../Version.cmake")

TRILINOS_DRIVER_SETUP()

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

IF(${Trilinos_REPOSITORY_BRANCH} STREQUAL "master")
  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT_DEV
    ctest_linux_nightly_serial_development_opt_gcc_4.7.2_rhel6-x86_64.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )

#if(${DAYOFWEEK} MATCHES "Saturday")
#endif()

ELSE()
#release tests would go here
ENDIF()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()

