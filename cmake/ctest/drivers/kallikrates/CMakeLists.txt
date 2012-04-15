INCLUDE("../../../../Version.cmake")

TRILINOS_DRIVER_SETUP()

IF(${Trilinos_REPOSITORY_BRANCH} STREQUAL "master")
  #Nightly tests
  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT
    ctest_windows_nightly_serial_development.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
  
   TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT_SHARED
    ctest_windows_nightly_serial_development_shared.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 60
  )
 
ELSE()
  #Release tests
  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT_RELEASE
    ctest_windows_nightly_serial_release.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_OPT_RELEASE_SHARED
    ctest_windows_nightly_serial_release_shared.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 60
  )
ENDIF()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
