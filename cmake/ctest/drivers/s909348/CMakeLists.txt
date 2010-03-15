# s909348 driver

#get the date for use in log files
execute_process(COMMAND date "+%m-%d-%Y" OUTPUT_VARIABLE DATE)

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(EG_EXE /Users/bmpersc/bin/eg)
set(BASEDIR /Users/bmpersc/nightly/Trilinos.base/development)
set(workingdir "${TrilinosDriver_BINARY_DIR}")
set(DRIVER_SCRIPT_DIR "${TrilinosDriver_SOURCE_DIR}/s909348")
set(TRILINOS_REPOSITORY_LOCATION "software.sandia.gov:/space/git/Trilinos")

if(TEST_RELEASE)
  add_test(mpi_release
    ${CTEST_EXE}
    -S
    "${DRIVER_SCRIPT_DIR}/ctest_mac_nightly_mpi_release_s909348.cmake"
    -VV
    --output_log
    "${workingdir}/MPI_RELEASE_10.0_${DATE}.log"
    )

  add_test(serial_release
    ${CTEST_EXE}
    -S
    "${DRIVER_SCRIPT_DIR}/ctest_mac_nightly_serial_release_s909348.cmake"
    -VV
    --output-log
    "${workingdir}/SERIAL_RELEASE_10.0_${DATE}.log"
    )
else()
  add_test(mpi_optimized
    ${CTEST_EXE}
    -S
    "${DRIVER_SCRIPT_DIR}/ctest_mac_nightly_mpi_opt_s909348.cmake"
    -VV
    --output_log
    "${workingdir}/MPI_OPT_DEV_${DATE}.log"
    )

  add_test(serial_debug
    ${CTEST_EXE}
    -S
    "${DRIVER_SCRIPT_DIR}/ctest_mac_nightly_serial_debug_s909348.cmake"
    -VV
    --output-log
    "${workingdir}/SERIAL_DEBUG_DEV_${DATE}.log"
    )
endif()
