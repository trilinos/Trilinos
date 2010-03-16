# dash2win64 driver

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

add_test(serial_release
  "${CTEST_EXE}"
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_windows_nightly_serial_release.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_msvc_nightly_serial_optimized_dash2win64.out"
)
