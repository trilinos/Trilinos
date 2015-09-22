# p90n03 driver

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

add_test(serial_debug
  "${CTEST_EXE}"
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_aix_nightly_serial_debug.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_aix_nightly_serial_debug_p90n03.out"
)
