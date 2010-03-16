# gabriel driver

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Now we add our tests
add_test(mpi_debug
  "${CTEST_EXE}"
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_mpi_debug_gabriel.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_linux_nightly_mpi_debug_gabriel.out"
)

add_test(kill_orted killall -s 9 orted)
set_property(TEST kill_orted PROPERTY DEPENDS mpi_debug)

add_test(serial_release
  "${CTEST_EXE}"
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_serial_release_gabriel.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_linux_nightly_serial_release_gabriel.out"
)

add_test(serial_debug_boost_tracing
  "${CTEST_EXE}"
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_serial_debug_boost_tracing_gabriel.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_linux_nightly_serial_debug_boost_tracing_gabriel.out"
)

add_test(notify_finish /home/rabartl/mailmsg.py "Finished nightly Trilinos CMake tests gabriel: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos")
set_property(TEST notify_finish PROPERTY DEPENDS mpi_debug serial_release serial_debug_boost_tracing)
