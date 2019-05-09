TRILINOS_DRIVER_SETUP()

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

if($ENV{CTEST_CONFIGURATION} MATCHES "default")

  # TRILINOS_DRIVER_ADD_DASHBORAD(
  #   OPENMPI_1.10.0_RELEASE_Tpetra
  #   ctest_linux_nightly_mpi_release_tpetra_trappist.gcc.cmake
  #   CTEST_INSTALLER_TYPE release
  #   RUN_SERIAL
  #   TIMEOUT_MINUTES 330
  #   )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.10.0_RELEASE_MueLu
    ctest_linux_nightly_mpi_release_muelu_trappist.gcc.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  # ETI tests
  TRILINOS_DRIVER_ADD_DASHBOARD(
    gcc_5.3.0_OPENMPI_1.10.1_RELEASE_NO-INT_NO-SERIAL_EXPERIMENTAL
    ctest_linux_nightly_mpi_release_muelu_no_int_no_serial_openmp_trappist.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  # REFACTOR tests
  TRILINOS_DRIVER_ADD_DASHBOARD(
    Linux_gcc_5.3.0_OPENMPI_1.10.1_KOKKOS_REFACTOR_OPENMP_EXPERIMENTAL
    ctest_linux_nightly_mpi_experimental_openmp_refactor_muelu_trappist.gcc.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

elseif($ENV{CTEST_CONFIGURATION} MATCHES "clang")

  TRILINOS_DRIVER_ADD_DASHBOARD(
    CLANG_OPENMPI_1.10.0_RELEASE
    ctest_linux_nightly_mpi_release_muelu_trappist.clang.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )
endif()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
