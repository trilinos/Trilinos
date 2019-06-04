TRILINOS_DRIVER_SETUP()

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

# serial valgrind tests run on Sunday night
# parallel valgrind tests run every night except Sunday
# coverage tests run on Saturday night

if(${DAYOFWEEK} MATCHES "Saturday")
  SET(RUN_HEAVY_MEMORY_TESTS ON)
  SET(RUN_COVERAGE_TESTS OFF)
  SET(RUN_CUDA_OLD_TESTS ON)
elseif(${DAYOFWEEK} MATCHES "Sunday")
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
    SERIAL_DEBUG_NO-EPETRA
    ctest_linux_nightly_serial_debug_muelu_tpetra_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_NO-TPETRA
    ctest_linux_nightly_serial_debug_muelu_epetra_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  #JHU: 2018-Apr-11 migrated to enigma
  #TRILINOS_DRIVER_ADD_DASHBOARD(
  #  OPENMPI-1.10.1_RELEASE_DEFAULT
  #  ctest_linux_nightly_mpi_release_muelu_geminga.cmake
  #  CTEST_INSTALLER_TYPE release
  #  RUN_SERIAL
  #  TIMEOUT_MINUTES 330
  #  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_COMPLEX_DEFAULT
    ctest_linux_nightly_mpi_release_complex_muelu_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  # ETI tests
  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_NO-EPETRA_NO-INT_NO-SERIAL_OPENMP_EXPERIMENTAL
    ctest_linux_nightly_mpi_release_muelu_tpetra_no_int_no_serial_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  #JHU 2018-Apr-11 migrated to enigma
  #TRILINOS_DRIVER_ADD_DASHBOARD(
  #  OPENMPI-1.10.1_RELEASE_NO-SERIAL_OPENMP_EXPERIMENTAL
  #  ctest_linux_nightly_mpi_release_muelu_no_serial_geminga.cmake
  #  CTEST_INSTALLER_TYPE release
  #  RUN_SERIAL
  #  TIMEOUT_MINUTES 330
  #  )

  #JHU 2018-Apr-11 migrated to enigma
  #TRILINOS_DRIVER_ADD_DASHBOARD(
  #  OPENMPI-1.10.1_RELEASE_NO-EPETRA_NO-INT_EXPERIMENTAL
  #  ctest_linux_nightly_mpi_release_muelu_no_int_geminga.cmake
  #  CTEST_INSTALLER_TYPE release
  #  RUN_SERIAL
  #  TIMEOUT_MINUTES 330
  #  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_NO-EPETRA_NO-SERIAL_OPENMP_EXPERIMENTAL
    ctest_linux_nightly_mpi_release_muelu_no_epetra_no_serial_openmp_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )

  #JHU 2018-Apr-12 migrated to enigma
  #TRILINOS_DRIVER_ADD_DASHBOARD(
  #  OPENMPI-1.10.1_RELEASE_NO-INT_NO-SERIAL_OPENMP_EXPERIMENTAL
  #  ctest_linux_nightly_mpi_release_muelu_no_int_no_serial_openmp_geminga.cmake
  #  CTEST_INSTALLER_TYPE release
  #  RUN_SERIAL
  #  TIMEOUT_MINUTES 330
  #  )

  #JHU 2018-Apr-12 migrated to enigma
  #TRILINOS_DRIVER_ADD_DASHBOARD(
  #  OPENMPI-1.10.1_RELEASE_NO-INT_OPENMP_EXPERIMENTAL
  #  ctest_linux_nightly_mpi_release_muelu_no_int_openmp_geminga.cmake
  #  CTEST_INSTALLER_TYPE release
  #  RUN_SERIAL
  #  TIMEOUT_MINUTES 330
  #  )

  # Tpetra Experimental Testing
  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI-1.10.1_RELEASE_NO-EPETRA_NO-INT_EXPERIMENTAL
    ctest_linux_nightly_mpi_release_tpetra_no_int_experimental_geminga.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
    )
  
  # Misc tests
  if(RUN_HEAVY_MEMORY_TESTS)

    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_DEBUG_VALGRIND
      ctest_linux_nightly_serial_debug_valgrind_muelu_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

    TRILINOS_DRIVER_ADD_DASHBOARD(
      OPENMPI-1.10.1_DEBUG_DEFAULT
      ctest_linux_nightly_mpi_debug_muelu_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

  else()

    TRILINOS_DRIVER_ADD_DASHBOARD(
      SERIAL_DEBUG_DEFAULT
      ctest_linux_nightly_serial_debug_muelu_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

    TRILINOS_DRIVER_ADD_DASHBOARD(
      OPENMPI-1.10.1_DEBUG_VALGRIND
      ctest_linux_nightly_mpi_debug_valgrind_muelu_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

  endif()

  if(RUN_COVERAGE_TESTS)

    TRILINOS_DRIVER_ADD_DASHBOARD(
      OPENMPI-1.10.1_DEBUG_COVERAGE_EXPERIMENTAL
      ctest_linux_nightly_mpi_debug_muelu_coverage_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 330
      )

  endif()

elseif($ENV{CTEST_CONFIGURATION} MATCHES "nvcc_wrapper")

  if(RUN_CUDA_OLD_TESTS)

    #TRILINOS_DRIVER_ADD_DASHBOARD(
    #  OPENMPI-1.8.7_DEBUG_KOKKOS-REFACTOR_CUDA_8.0.44
    #  ctest_linux_nightly_mpi_debug_muelu_kokkos_refactor_cuda_geminga.cmake
    #  CTEST_INSTALLER_TYPE release
    #  RUN_SERIAL
    #  TIMEOUT_MINUTES 420
    #  )

    TRILINOS_DRIVER_ADD_DASHBOARD(
      OPENMPI-1.8.7_RELEASE_KOKKOS-REFACTOR_CUDA-8.0.44
      ctest_linux_nightly_mpi_release_muelu_kokkos_refactor_cuda_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 420
      )

    TRILINOS_DRIVER_ADD_DASHBOARD(
      RELEASE_AMGX_CUDA-8.0.44
      ctest_linux_nightly_mpi_release_muelu_amgx_cuda_geminga.cmake
      CTEST_INSTALLER_TYPE release
      RUN_SERIAL
      TIMEOUT_MINUTES 420
      )

  endif()


endif()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
