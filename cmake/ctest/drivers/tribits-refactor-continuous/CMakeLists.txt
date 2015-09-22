TRILINOS_DRIVER_SETUP()

find_program(PYTHON_EXE python)

TRILINOS_DRIVER_ADD_DASHBOARD(
  TRIBITS_MPI 
  ctest_tribits_refactor.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 240
)

SET(INSTALLATION_DIR "${CMAKE_BINARY_DIR}/../installation/installed")
SET(INSTALLATION_WORKING_DIR "${CMAKE_BINARY_DIR}/../installation")
ADD_TEST(Install_Trilinos ${PYTHON_EXE}
  "${CMAKE_CURRENT_SOURCE_DIR}/build_tarball.py"
  "--no-dashboard"
  "--working-dir=${INSTALLATION_WORKING_DIR}"
  "--install-dir=${INSTALLATION_DIR}"
  "--disable-boost"
  "--disable-netcdf"
  "--repository=software.sandia.gov:/home/rabartl/PROJECTS/backups/Trilinos"
  "--branch=tribits_refactor"
  )

TRILINOS_DRIVER_ADD_DASHBOARD(
  TRIBITS_INSTALL
  ctest_tribits_installation_testing.cmake
  CTEST_INSTALLER_TYPE release
  DEPENDS Install_Trilinos
  REQUIRED_FILES ${INSTALLATION_DIR}/lib/cmake/Trilinos/TrilinosConfig.cmake
)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()

