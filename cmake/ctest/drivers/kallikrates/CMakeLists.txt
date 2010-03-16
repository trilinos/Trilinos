# kallikrates driver

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Now run ctest on each of the ctest build scripts for this machine
add_test(serial_release
  "${CTEST_EXE}"
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_windows_nightly_serial_release.cmake"
  -VV
  --output-log
  "${workingdir}/ctest_msvc_nightly_serial_release_kallikrates.out"
)

# disabling temporarily since the mpi build is unstable on windows right now.
#add_test(mpi_release
#  "${CTEST_EXE}"
#  -S
#  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_windows_nightly_mpi_release.cmake"
#  -VV
#  --output-log
#  "${workingdir}/ctest_msvc_nightly_mpi_release_kallikrates.out"
#)
