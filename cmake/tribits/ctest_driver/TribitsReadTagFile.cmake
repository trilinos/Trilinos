# Set policy here instead of including TribitCTestPolicis.cmake since we want
# this to be a standalone module
cmake_policy(SET CMP0007 NEW)


# @FUNCTION: tribits_read_ctest_tag_file()
#
# Read in the <build>/Testing/TAG file contents
#
# Usage::
#
#   tribits_read_ctest_tag_file( <tagFileIn>
#     <buildStartTimeOut> <cdashGroupOut> <cdashModelOut> )
#
function(tribits_read_ctest_tag_file  tagFileIn
    buildStartTimeOut  cdashGroupOut  cdashModelOut
  )
  file(READ "${tagFileIn}" tagFileStr)
  string(REPLACE "\n" ";" tagFileStrList "${tagFileStr}")
  list(GET tagFileStrList 0 buildStartTime)
  list(GET tagFileStrList 1 cdashGroup)
  list(GET tagFileStrList 2 cdashModel)
  set(${buildStartTimeOut} "${buildStartTime}" PARENT_SCOPE)
  set(${cdashGroupOut} "${cdashGroup}" PARENT_SCOPE)
  set(${cdashModelOut} "${cdashModel}" PARENT_SCOPE)
endfunction()
