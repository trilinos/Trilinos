TRILINOS_DRIVER_SETUP()

if(DEFINED ENV{DO_TPETRA_TESTING})

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_TPETRA_DEPRECATED_CODE_OFF
    ctest_linux_nightly_mpi_release_tpetra_rocketman.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_TPETRA_DEPRECATED_CODE_OFF_DOWNSTREAM_ENABLED_NO_EPETRA
    mpi_release_tpetra_deprecated_code_off_downstream_enabled_no_epetra.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_TPETRA_DEPRECATED_CODE_OFF_DOWNSTREAM_ENABLED
    mpi_release_tpetra_deprecated_code_off_downstream_enabled.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_TPETRA_DEPRECATED_CODE_OFF_DOWNSTREAM_ENABLED_GO_INT
    mpi_release_tpetra_deprecated_code_off_downstream_enabled_GO_int.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_TPETRA_PERFORMANCE
    ctest_linux_experimental_mpi_release_tpetra_performance_rocketman.cmake
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
