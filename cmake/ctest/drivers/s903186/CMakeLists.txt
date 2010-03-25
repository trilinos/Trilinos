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

# Ross Bartlett: We need to define a helper function called something like
# SET_CTEST_DRIVER_TEST(...) that takes care of setting the correct test
# properties logic.  Also, we need to inject the update of the Trilinos
# sources and the inner CMake/CTest.
#
# There is just too much duplication in the raw add_test(...) examples above.
#
# The updated CMakeLists.txt file should look something like:
#
#
#  INCLUDE(TrilinosStandardDriverSupport)
# 
#  # Set up default Trilinos and CMake update tests etc.
#  TRILINOS_STANDARD_DRIVER_DEFAULT_SETUP()
#
#  SET_CTEST_DRIVER_TEST(MPI_OPT
#      ctest_mac_nightly_mpi_opt_s903186.cmake
#      TIMEOUT_MIN 60
#      )
#
#  SET_CTEST_DRIVER_TEST(SERIAL_DEBUG
#      ctest_mac_nightly_serial_debug_s903186.cmake
#      TIMEOUT_MIN 90
#      )
#
#
