# godel driver

set(CTEST_EXE "${CMAKE_CTEST_COMMAND}")
set(workingdir "${TrilinosDriver_BINARY_DIR}")

# Get the day of the week
execute_process(COMMAND date +%A OUTPUT_VARIABLE DAYOFWEEK)

if(${DAYOFWEEK} MATCHES "Saturday")
  set(RUN_REGULAR_TESTS ON)
  set(RUN_COVERAGE_TESTS ON)
  set(RUN_MEMCHECK_TESTS OFF)
elseif(${DAYOFWEEK} MATCHES "Sunday")
  set(RUN_REGULAR_TESTS ON)
  set(RUN_COVERAGE_TESTS OFF)
  set(RUN_MEMCHECK_TESTS OFF)
else()
  set(RUN_REGULAR_TESTS ON)
  set(RUN_COVERAGE_TESTS OFF)
  set(RUN_MEMCHECK_TESTS OFF)
endif()

#add_test(email_results /home/rabartl/mailmsg.py "Finished nightly Trilinos tests godel: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos")

if(RUN_REGULAR_TESTS)
  add_test(serial_performance
    ${CTEST_EXE}
    -S
    "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_serial_performance_godel.cmake"
    -VV
    --output-log
    "${workingdir}/ctest_linux_nightly_serial_performance_godel.out"
  )

#echo
#echo "Doing mpi optimized build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_optimized_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_mpi_optimized_godel.out
#
#echo
#echo "Doing serial debug build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_debug_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_serial_debug_godel.out
#
#echo
#echo "Doing mpi optimized shared library build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_optimized_shared_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_mpi_optimized_shared_godel.out

#echo
#echo "Doing serial optimized implicit instantiation build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_opt_impl_instant_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_serial_opt_impl_instant_godel.out

#echo
#echo "Doing mpi debug build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_debug_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_mpi_debug_godel.out

  add_test(mpi_opt_zoltan_c
    ${CTEST_EXE}
    -S
    "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_mpi_opt_zoltan_c_godel.cmake"
    -VV
    --output-log
    "${workingdir}/ctest_linux_nightly_mpi_opt_zoltan_c_godel.out"
  )

  add_test(serial_debug_icpc
    ${CTEST_EXE}
    -S
    "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_serial_debug_icpc_godel.cmake"
    -VV
    --output-log
    "${workingdir}/ctest_linux_nightly_serial_debug_icpc_godel.out"
  )

  add_test(serial_release
    ${CTEST_EXE}
    -S
    "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_serial_release_icpc_godel.cmake"
    -VV
    --output-log
    "${workingdir}/ctest_linux_nightly_serial_release_icpc_godel.out"
  )
  
  #set_property(TEST email_results APPEND PROPERTY DEPENDS
    #serial_performance mpi_opt_zoltan_c serial_debug_icpc serial_release)
endif(RUN_REGULAR_TESTS)

if(RUN_COVERAGE_TESTS)
  add_test(serial_debug_coverage
    ${CTEST_EXE}
    -S
    "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_serial_debug_coverage_godel.cmake"
    -VV
    --output-log
    "${workingdir}/ctest_linux_nightly_serial_debug_coverage_godel.out"
  )

  add_test(mpi_debug_coverage
    ${CTEST_EXE}
    -S
    "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_mpi_debug_coverage_godel.cmake"
    -VV
    --output-log
    "${workingdir}/ctest_linux_nightly_mpi_debug_coverage_godel.out"
  )
  
  #set_property(TEST email_results APPEND PROPERTY DEPENDS serial_debug_coverage mpi_debug_coverage)
endif(RUN_COVERAGE_TESTS)

if(RUN_MEMCHECK_TESTS)

  add_test(serial_debug_memcheck
    ${CTEST_EXE}
    -S
    "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_serial_debug_memcheck_godel.cmake"
    -VV
    --output-log
    "${workingdir}/ctest_linux_nightly_serial_debug_memcheck_godel.out"
  )

  add_test(mpi_debug_memcheck
    ${CTEST_EXE}
    -S
    "${CMAKE_CURRENT_SOURCE_DIR}/ctest_linux_nightly_mpi_debug_memcheck_godel.cmake"
    -VV
    --output-log
    "${workingdir}/ctest_linux_nightly_mpi_debug_memcheck_godel.out"
  )
  
  add_test(kill_memcheck_processes killall -s 9 memcheck)
  set_property(TEST kill_memcheck_processes PROPERTY DEPENDS mpi_debug_memcheck serial_debug_memcheck)
  #set_property(TEST email_results APPEND PROPERTY DEPENDS mpi_debug_memcheck serial_debug_memcheck)
endif(RUN_MEMCHECK_TESTS)

#forcing all files/directories to be group accessible
#chgrp -R trilinos-dev *
#chmod -R g+w *
