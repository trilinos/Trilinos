# godel driver

TRILINOS_DRIVER_SETUP()

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

if(${DAYOFWEEK} MATCHES "Saturday")
  set(RUN_REGULAR_TESTS ON)
  set(RUN_COVERAGE_TESTS ON)
  set(RUN_MEMCHECK_TESTS OFF)
elseif(${DAYOFWEEK} MATCHES "Sunday")
  set(RUN_REGULAR_TESTS ON)
  set(RUN_COVERAGE_TESTS OFF)
  set(RUN_MEMCHECK_TESTS OFF)
else()
  set(RUN_REGULAR_TESTS ON)
  set(RUN_COVERAGE_TESTS OFF)
  set(RUN_MEMCHECK_TESTS OFF)
endif()

if(RUN_REGULAR_TESTS)
  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_PERFORMANCE
    ctest_linux_nightly_serial_performance_godel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_OPT_ZOLTAN_C
    ctest_linux_nightly_mpi_opt_zoltan_c_godel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_OPT_ZOLTAN_C_64
    ctest_linux_nightly_mpi_opt_zoltan_c64_godel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_ICPC
    ctest_linux_nightly_serial_debug_icpc_godel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_RELEASE_ICPC
    ctest_linux_nightly_serial_release_icpc_godel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )
  
endif()

if(RUN_COVERAGE_TESTS)
  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_COVERAGE
    ctest_linux_nightly_serial_debug_coverage_godel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_DEBUG_COVERAGE
    ctest_linux_nightly_mpi_debug_coverage_godel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 300
  )
  
endif()

if(RUN_MEMCHECK_TESTS)

  TRILINOS_DRIVER_ADD_DASHBOARD(
    SERIAL_DEBUG_MEMCHECK
    ctest_linux_nightly_serial_debug_memcheck_godel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )

  TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_DEBUG_MEMCHECK
    ctest_linux_nightly_mpi_debug_memcheck_godel.cmake
    RUN_SERIAL
    TIMEOUT_MINUTES 240
  )
 
#not sure what to do with these maybe just a command? I believe that the new
#testing system is suposed to clean up spawned processes completely.
#  
#  add_test(kill_memcheck_processes killall -s 9 memcheck)
#  set_property(TEST kill_memcheck_processes PROPERTY DEPENDS mpi_debug_memcheck serial_debug_memcheck)
endif()

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
