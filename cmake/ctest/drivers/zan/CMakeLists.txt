INCLUDE("../../../../Version.cmake")

TRILINOS_DRIVER_SETUP()

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

IF(${Trilinos_REPOSITORY_BRANCH} STREQUAL "master")


TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_OPT_DEV
  ctest_linux_nightly_mpi_development_opt_zan.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 300
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_DEBUG_DEV
  ctest_linux_nightly_serial_development_debug_zan.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 360
)

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  SERIAL_DEBUG_ICPC
#  ctest_linux_nightly_serial_debug_icpc_zan.cmake
#  CTEST_INSTALLER_TYPE release
#  RUN_SERIAL
#  TIMEOUT_MINUTES 540
#)

#Only run these builds once a week
if(${DAYOFWEEK} MATCHES "Saturday")
TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_DEBUG_ICPC
  ctest_linux_nightly_mpi_debug_icpc_zan.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 300
)
endif()

if(${DAYOFWEEK} MATCHES "Sunday")
TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_OPT_ICPC 
  ctest_linux_nightly_serial_opt_icpc_zan.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 240
)
endif()

ELSE()

#TRILINOS_DRIVER_ADD_DASHBOARD(
#  SERIAL_OPT_RELEASE_ICPC 
#  ctest_linux_nightly_serial_opt_release_icpc_zan.cmake
#  CTEST_INSTALLER_TYPE release
#  RUN_SERIAL
#  TIMEOUT_MINUTES 240
#)

ENDIF()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
