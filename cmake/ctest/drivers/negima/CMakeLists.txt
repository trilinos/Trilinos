TRILINOS_DRIVER_SETUP()

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

# serial valgrind tests run on Sunday night
# parallel valgrind tests run every night except Sunday
# coverage tests run on Saturday night

if(${DAYOFWEEK} MATCHES "Sunday")
  SET(RUN_HEAVY_MEMORY_TESTS ON)
  SET(RUN_COVERAGE_TESTS OFF)
elseif(${DAYOFWEEK} MATCHES "Saturday")
  SET(RUN_HEAVY_MEMORY_TESTS OFF)
  SET(RUN_COVERAGE_TESTS ON)
else()
  SET(RUN_HEAVY_MEMORY_TESTS OFF)
  SET(RUN_COVERAGE_TESTS OFF)
endif()

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_MueLu_Tpetra
  ctest_linux_nightly_serial_debug_muelu_tpetra_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_MueLu_Epetra
  ctest_linux_nightly_serial_debug_muelu_epetra_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_MueLu_ExtraTypes
  ctest_linux_nightly_serial_debug_muelu_extratypes_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_MueLu_ExtraTypes_EI
  ctest_linux_nightly_mpi_debug_muelu_extratypes_ei_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  SERIAL_DEBUG_MueLu_ExtraTypes_EI_CLANG
#  ctest_linux_nightly_mpi_debug_muelu_extratypes_ei_negima_clang.cmake
#  CTEST_INSTALLER_TYPE release
#  RUN_SERIAL
#  TIMEOUT_MINUTES 330
#)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_RELEASE_MueLu_Experimental
  ctest_linux_nightly_serial_release_muelu_experimental_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  OPENMPI_1.8.3_RELEASE_MueLu
  ctest_linux_nightly_mpi_release_muelu_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  OPENMPI_1.8.3_DEBUG_MueLu_KLU2
  ctest_linux_nightly_mpi_debug_muelu_klu2_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

if(RUN_HEAVY_MEMORY_TESTS)

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_MueLu_Valgrind
    ctest_linux_nightly_serial_debug_valgrind_muelu_negima.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
  )

else()

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_MueLu
    ctest_linux_nightly_serial_debug_muelu_negima.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.8.3_DEBUG_MueLu_Valgrind
    ctest_linux_nightly_mpi_debug_valgrind_muelu_negima.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
  )

endif()

if(RUN_COVERAGE_TESTS)

  TRILINOS_DRIVER_ADD_DASHBOARD(
    OPENMPI_1.8.3_DEBUG_MueLu_Coverage
    ctest_linux_nightly_mpi_debug_muelu_coverage_negima.cmake
    CTEST_INSTALLER_TYPE release
    RUN_SERIAL
    TIMEOUT_MINUTES 330
  )

endif()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
