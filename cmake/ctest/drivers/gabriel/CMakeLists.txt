# gabriel driver

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")

set(workingdir "${TrilinosDriver_BINARY_DIR}")
set(DRIVER_SCRIPT_DIR="${workingdir}/Trilinos/cmake/ctest/drivers/gabriel")
set(BRANCH trilinos-release-10-0-branch)
set(TRILINOS_REPOSITORY_LOCATION software.sandia.gov:/space/git/nightly/Trilinos)

# checkout the basics from Trilinos needed to run the dashboard including
# this script.
if(NOT EXISTS "${workingdir}/Trilinos")
  message("Cloning the repository because none exists yet.")
  execute_process(COMMAND git clone ${TRILINOS_REPOSITORY_LOCATION}
    WORKING_DIRECTORY "${workingdir}")
  execute_process(COMMAND git checkout --track origin/${BRANCH}
    WORKING_DIRECTORY "${workingdir}/Trilinos")
endif()

# Add targets to update the local repository
add_custom_target(gitpull ALL
  COMMAND eg pull
  WORKING_DIRECTORY "${workingdir}/Trilinos"
  COMMENT "Updating the local git repo")

# Now we add our tests
add_test(mpi_debug
  "${CTEST_EXE}"
  -S
  "${DRIVER_SCRIPT_DIR}/ctest_linux_nightly_mpi_debug_gabriel.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_linux_nightly_mpi_debug_gabriel.out"
)

add_test(kill_orted killall -s 9 orted)
set_tests_properties(kill_orted PROPERTIES DEPENDS mpi_debug)

add_test(serial_release
  "${CTEST_EXE}
  -S
  "${DRIVER_SCRIPT_DIR}/ctest_linux_nightly_serial_release_gabriel.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_linux_nightly_serial_release_gabriel.out"
)

add_test(serial_debug_boost_tracing
  "${CTEST_EXE}"
  -S
  "${DRIVER_SCRIPT_DIR}/ctest_linux_nightly_serial_debug_boost_tracing_gabriel.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_linux_nightly_serial_debug_boost_tracing_gabriel.out"
)

add_test(notify_finish /home/rabartl/mailmsg.py "Finished nightly Trilinos CMake tests gabriel: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos")
set_tests_properties(notify_finish DEPENDS mpi_debug serial_release serial_debug_boost_tracing)
