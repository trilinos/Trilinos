TRILINOS_DRIVER_SETUP()

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_MueLu
  ctest_linux_nightly_serial_debug_muelu_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  OPENMPI_1.8.3_DEBUG_MueLu
  ctest_linux_nightly_mpi_debug_muelu_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

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
  OPENMPI_1.8.3_RELEASE_MueLu
  ctest_linux_nightly_mpi_debug_muelu_klu2_negima.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
