TRILINOS_DRIVER_SETUP()

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_OPT_DEV
  ctest_linux_nightly_mpi_development_opt_trilinos-test2.cmake
  CTEST_INSTALLER_TYPE dev
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_OPT_DEV_NOEI
  ctest_linux_nightly_mpi_development_opt_noei_trilinos-test2.cmake
  CTEST_INSTALLER_TYPE dev
)

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  MPI_DEBUG_DEV
#  ctest_linux_nightly_mpi_development_debug_trilinos-test2.cmake
#  CTEST_INSTALLER_TYPE dev
#)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_OPT_DEV_SHARED
  ctest_linux_nightly_mpi_development_opt_shared_trilinos-test2.cmake
  CTEST_INSTALLER_TYPE dev
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_DEV
  ctest_linux_nightly_serial_development_debug_trilinos-test2.cmake
  CTEST_INSTALLER_TYPE dev
)

SET(INSTALLATION_DIR "${CMAKE_BINARY_DIR}/../MPI_OPT_DEV_INSTALL/installed")
SET(INSTALLATION_WORKING_DIR "${CMAKE_BINARY_DIR}/../MPI_OPT_DEV_INSTALL")
message("MPI_BASE_DIR = '${MPI_BASE_DIR}'")
find_program(PYTHON_EXE python)
ADD_TEST(Install_Trilinos ${PYTHON_EXE}
  "${CMAKE_CURRENT_SOURCE_DIR}/nightly_create_tarball.py"
  "--enable-mpi"
  "--mpi-base-path=/home/trilinos"
  "--no-dashboard"
  "--working-dir=${INSTALLATION_WORKING_DIR}"
  "--install-dir=${INSTALLATION_DIR}"
  "--boost-dir=/home/trilinos/include"
  "--verbose-errors"
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_OPT_DEV_INSTALL
  ctest_linux_nightly_mpi_development_opt_install_trilinos-test2.cmake
  CTEST_INSTALLER_TYPE dev
  DEPENDS Install_Trilinos
  REQUIRED_FILES ${INSTALLATION_DIR}/lib/cmake/Trilinos/TrilinosConfig.cmake
)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
