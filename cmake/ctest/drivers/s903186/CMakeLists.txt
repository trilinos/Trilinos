# s903186 driver

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")

#get the date for use in log files
execute_process(COMMAND date "+%m-%d-%Y" OUTPUT_VARIABLE DATE)

set(EG_EXE /Users/jmwille/bin/eg)
set(BASEDIR /Users/jmwille/TrilinosTestHarness)
set(DRIVER_SCRIPT_DIR "${TrilinosDriver_SOURCE_DIR}/s903186")
set(TRILINOS_REPOSITORY_LOCATION "software.sandia.gov:/space/git/Trilinos")
set(workingdir "${TrilinosDriver_BINARY_DIR}")
#export PYTHONPATH=/Users/jmwille/install/lib/python2.5/site-packages

add_test(mpi_optimized
  ${CTEST_EXE}
  -S
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest_mac_nightly_mpi_opt_s903186.cmake"
  -VV
  --output-log
  "${workingdir}/MPI_OPT_DEV_${DATE}.log"
  )

add_test(serial_debug
 ${CTEST_EXE}
 -S
 "${CMAKE_CURRENT_SOURCE_DIR}/ctest_mac_nightly_serial_debug_s903186.cmake"
 -VV
 --output-log
 "${workingdir}/SERIAL_DEBUG_DEV_${DATE}.log"
  )
