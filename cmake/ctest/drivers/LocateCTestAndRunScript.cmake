# This script is a cmake -P script that should be called with the following
# variables defined:
#
#    -D binary_dir=${CMAKE_CURRENT_BINARY_DIR}
#    -D source_dir=${CMAKE_CURRENT_SOURCE_DIR}
#    -D ctest_type=${ctest_type}
#    -D scriptname=${scriptname}
#    -D TD_BASE_DIR=${TD_BASE_DIR}
#    -D testname=${testname}
#
# It looks recursively under ${TD_BASE_DIR}/tools/cmake-${ctest_type} for a
# ctest executable and uses it to drive a ctest -S script to run a dashboard
# for Trilinos. It has to be run this way indirectly through a cmake -P script
# because the desired ctest executable location is unknown at CMake configure
# time of the TrilinosDriver project.

if(NOT CTEST_EXE)
  if(WIN32)
    set(ctest_filename "ctest.exe")
  else()
    set(ctest_filename "ctest")
  endif()

  message("globbing for '${ctest_filename}'...")

  file(GLOB_RECURSE CTEST_EXE
    "${TD_BASE_DIR}/tools/cmake-${ctest_type}/${ctest_filename}")
endif()

if(NOT CTEST_EXE)
  message(FATAL_ERROR "error: '${ctest_type}' ctest could not be found...")
endif()

if(NOT EXISTS "${CTEST_EXE}")
  message(FATAL_ERROR "error: CTEST_EXE='${CTEST_EXE}' does not exist...")
endif()

message("CTEST_EXE='${CTEST_EXE}'")

message("=========== variables ===========")
message("binary_dir='${binary_dir}'")
message("source_dir='${source_dir}'")
message("ctest_type='${ctest_type}'")
message("scriptname='${scriptname}'")
message("TD_BASE_DIR='${TD_BASE_DIR}'")
message("testname='${testname}'")
message("=================================")

message("========== environment ==========")
execute_process(COMMAND ${CMAKE_COMMAND} -E environment)
message("=================================")

message("============ script =============")
message("executing ctest -S '${scriptname}' for test '${testname}'...")

execute_process(COMMAND ${CTEST_EXE}
  -S
  "${source_dir}/${scriptname}"
  -V
  --output-log
  "${binary_dir}/${testname}.log"
  RESULT_VARIABLE rv
  )

if(NOT "${rv}" STREQUAL "0")
  message("warning: calling ctest -S script failed with '${rv}'")
endif()

message("=================================")
