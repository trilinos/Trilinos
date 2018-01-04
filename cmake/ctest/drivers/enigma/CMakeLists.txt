TRILINOS_DRIVER_SETUP()

TRILINOS_DRIVER_ADD_DASHBOARD(
  OPENMPI_1.10.1_RELEASE_DEFAULT
  ctest_linux_nightly_mpi_release_muelu_enigma.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  OPENMPI_1.10.1_DEBUG_NO-SUPERLU_BASKER
  ctest_linux_nightly_mpi_debug_muelu_basker_enigma.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  OPENMPI_1.10.1_DEBUG_NO-SUPERLU_KLU2
  ctest_linux_nightly_mpi_debug_muelu_klu2_enigma.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_LONG-LONG
  ctest_linux_nightly_mpi_debug_muelu_extratypes_ei_enigma.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_LONG-LONG_COMPLEX-DOUBLE
  ctest_linux_nightly_serial_debug_muelu_extratypes_enigma.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_RELEASE_EXPERIMENTAL
  ctest_linux_nightly_serial_release_muelu_experimental_enigma.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
  )


#TRILINOS_DRIVER_ADD_DASHBOARD(
#  SERIAL_RELEASE_MueLu_Experimental
#  ctest_linux_nightly_serial_release_muelu_experimental_typhon.cmake
#  CTEST_INSTALLER_TYPE release
#  RUN_SERIAL
#  TIMEOUT_MINUTES 330
#)

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  MPICH2_1.4.1p1_RELEASE_MueLu
#  ctest_linux_nightly_mpi_release_muelu_typhon.cmake
#  CTEST_INSTALLER_TYPE release
#  RUN_SERIAL
#  TIMEOUT_MINUTES 330
#)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
