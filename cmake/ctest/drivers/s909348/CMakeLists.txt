# s909348 driver

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

add_test(mpi_release
  ${CTEST_EXE}
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_mac_nightly_mpi_release_s909348.cmake"
  -VV
  --output_log
  "${workingdir}/MPI_RELEASE_10.0.log"
  )

add_test(serial_release
  ${CTEST_EXE}
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_mac_nightly_serial_release_s909348.cmake"
  -VV
  --output-log
  "${workingdir}/SERIAL_RELEASE_10.0.log"
  )

add_test(mpi_optimized
  ${CTEST_EXE}
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_mac_nightly_mpi_opt_s909348.cmake"
  -VV
  --output_log
  "${workingdir}/MPI_OPT_DEV.log"
  )

add_test(serial_debug
  ${CTEST_EXE}
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_mac_nightly_serial_debug_s909348.cmake"
  -VV
  --output-log
  "${workingdir}/SERIAL_DEBUG_DEV.log"
  )
