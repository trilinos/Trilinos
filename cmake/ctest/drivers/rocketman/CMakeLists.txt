TRILINOS_DRIVER_SETUP()

if(DEFINED ENV{DO_TPETRA_TESTING})

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_TPETRA_DEPRECATED_CODE_OFF
    ctest_linux_nightly_mpi_release_tpetra_rocketman.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

else()

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_DEFAULT
    ctest_linux_nightly_mpi_release_muelu_rocketman.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )


  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_AVATAR
    ctest_linux_experimental_mpi_release_avatar_rocketman.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

endif()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
