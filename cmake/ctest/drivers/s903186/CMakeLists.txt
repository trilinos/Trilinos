# s903186 driver

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

add_test(mpi_optimized
  ${CTEST_EXE}
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_mac_nightly_mpi_opt_s903186.cmake"
  -VV
  --output-log
  "${workingdir}/MPI_OPT_DEV.log"
  )

add_test(serial_debug
  ${CTEST_EXE}
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_mac_nightly_serial_debug_s903186.cmake"
  -VV
  --output-log
  "${workingdir}/SERIAL_DEBUG_DEV.log"
  )
