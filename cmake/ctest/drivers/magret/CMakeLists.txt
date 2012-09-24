TRILINOS_DRIVER_SETUP()

# TRILINOS_DRIVER_ADD_DASHBOARD(
#   SERIAL_DEBUG_MueLu
#   ctest_linux_nightly_serial_debug_muelu_magret.cmake
#   CTEST_INSTALLER_TYPE release
#   RUN_SERIAL
#   TIMEOUT_MINUTES 330
# )

# TRILINOS_DRIVER_ADD_DASHBOARD(
#   SERIAL_DEBUG_MueLu_Tpetra
#   ctest_linux_nightly_serial_debug_muelu_tpetra_magret.cmake
#   CTEST_INSTALLER_TYPE release
#   RUN_SERIAL
#   TIMEOUT_MINUTES 330
# )

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_MueLu_Epetra
  ctest_linux_nightly_serial_debug_muelu_epetra_magret.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

# TRILINOS_DRIVER_ADD_DASHBOARD(
#   SERIAL_DEBUG_MueLu_ExtraTypes
#   ctest_linux_nightly_serial_debug_muelu_extratypes_magret.cmake
#   CTEST_INSTALLER_TYPE release
#   RUN_SERIAL
#   TIMEOUT_MINUTES 330
# )

# TRILINOS_DRIVER_ADD_DASHBOARD(
#   SERIAL_DEBUG_MueLu_ExtraTypes_EI
#   ctest_linux_nightly_serial_debug_muelu_extratypes_ei_magret.cmake
#   CTEST_INSTALLER_TYPE release
#   RUN_SERIAL
#   TIMEOUT_MINUTES 330
# )

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
