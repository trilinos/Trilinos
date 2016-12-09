TRILINOS_DRIVER_SETUP()

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

# serial valgrind tests run on Sunday night
# parallel valgrind tests run every night except Sunday
# coverage tests run on Saturday night

if(${DAYOFWEEK} MATCHES "Sunday")
  SET(RUN_HEAVY_MEMORY_TESTS ON)
  SET(RUN_COVERAGE_TESTS OFF)
  SET(RUN_CUDA_OLD_TESTS ON)
elseif(${DAYOFWEEK} MATCHES "Saturday")
  SET(RUN_HEAVY_MEMORY_TESTS OFF)
  SET(RUN_COVERAGE_TESTS ON)
  SET(RUN_CUDA_OLD_TESTS OFF)
else()
  SET(RUN_HEAVY_MEMORY_TESTS OFF)
  SET(RUN_COVERAGE_TESTS OFF)
  SET(RUN_CUDA_OLD_TESTS OFF)
endif()

SET(RUN_CUDA_OLD_TESTS ON)

if($ENV{CTEST_CONFIGURATION} MATCHES "default")

  # Regular tests
  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_MueLu_Tpetra
    ctest_linux_nightly_serial_debug_muelu_tpetra_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_MueLu_Epetra
    ctest_linux_nightly_serial_debug_muelu_epetra_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.10.0_RELEASE_MueLu
    ctest_linux_nightly_mpi_release_muelu_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  # ETI tests
  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.10.0_RELEASE_DEV_MueLu_Tpetra_NO_INT_NO_SERIAL
    ctest_linux_nightly_mpi_release_muelu_tpetra_no_int_no_serial_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.10.0_RELEASE_DEV_MueLu_NO_SERIAL
    ctest_linux_nightly_mpi_release_muelu_no_serial_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.10.0_RELEASE_DEV_MueLu_NO_INT
    ctest_linux_nightly_mpi_release_muelu_no_int_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.10.0_RELEASE_DEV_MueLu_NO_EPETRA_NO_SERIAL
    ctest_linux_nightly_mpi_release_muelu_no_epetra_no_serial_openmp_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.10.0_RELEASE_DEV_MueLu_NO_INT_NO_SERIAL
    ctest_linux_nightly_mpi_release_muelu_no_int_no_serial_openmp_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.10.0_RELEASE_DEV_MueLu_NO_INT
    ctest_linux_nightly_mpi_release_muelu_no_int_openmp_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  # Refactor (kokkos) tests
  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.10.0_DEBUG_MueLu_Kokkos_Refactor_OpenMP
    ctest_linux_nightly_mpi_debug_muelu_kokkos_refactor_openmp_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 420
    )

  # Misc tests
  if(RUN_HEAVY_MEMORY_TESTS)

    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_DEBUG_MueLu_Valgrind
      ctest_linux_nightly_serial_debug_valgrind_muelu_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

    TRILINOS_DRIVER_ADD_DASHBOARD(
      OPENMPI_1.10.0_DEBUG_MueLu
      ctest_linux_nightly_mpi_debug_muelu_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

  else()

    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_DEBUG_MueLu
      ctest_linux_nightly_serial_debug_muelu_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

    TRILINOS_DRIVER_ADD_DASHBOARD(
      OPENMPI_1.10.0_DEBUG_MueLu_Valgrind
      ctest_linux_nightly_mpi_debug_valgrind_muelu_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

  endif()

  if(RUN_COVERAGE_TESTS)

    TRILINOS_DRIVER_ADD_DASHBOARD(
      OPENMPI_1.10.0_DEBUG_MueLu_Coverage
      ctest_linux_nightly_mpi_debug_muelu_coverage_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

  endif()

elseif($ENV{CTEST_CONFIGURATION} MATCHES "broken")

  # ETI tests
  TRILINOS_DRIVER_ADD_DASHBOARD(
    BROKEN_Tpetra_NO_INT_NO_SERIAL
    ctest_linux_nightly_mpi_release_muelu_tpetra_no_int_no_serial_broken_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    BROKEN_NO_SERIAL
    ctest_linux_nightly_mpi_release_muelu_no_serial_broken_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    BROKEN_NO_INT
    ctest_linux_nightly_mpi_release_muelu_no_int_broken_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    BROKEN_NO_EPETRA_NO_SERIAL
    ctest_linux_nightly_mpi_release_muelu_no_epetra_no_serial_openmp_broken_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    BROKEN_NO_INT_NO_SERIAL
    ctest_linux_nightly_mpi_release_muelu_no_int_no_serial_openmp_broken_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    BROKEN_NO_INT
    ctest_linux_nightly_mpi_release_muelu_no_int_openmp_broken_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

elseif($ENV{CTEST_CONFIGURATION} MATCHES "nvcc_wrapper")

  if(RUN_CUDA_OLD_TESTS)

    #TRILINOS_DRIVER_ADD_DASHBOARD(
    #  OPENMPI_1.10.0_DEBUG_MueLu_Kokkos_Refactor_CUDA
    #  ctest_linux_nightly_mpi_debug_muelu_kokkos_refactor_cuda_geminga.cmake
    #  CTEST_INSTALLER_TYPE release
    #  RUN_SERIAL
    #  TIMEOUT_MINUTES 420
    #  )

    TRILINOS_DRIVER_ADD_DASHBOARD(
      OPENMPI_1.10.0_RELEASE_MueLu_Kokkos_Refactor_CUDA
      ctest_linux_nightly_mpi_release_muelu_kokkos_refactor_cuda_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 420
      )

  endif()


endif()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
