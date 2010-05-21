INCLUDE(PrintVar)
INCLUDE(GlobalSet)

#
# Function that wraps the standard CMake/CTest MESSAGE(...) function call in
# order to allow unit testing to intercept the call.
#
FUNCTION(MESSAGE_WRAPPER)
  #MESSAGE("MESSAGE_WRAPPER: ${ARGN}")
  #PRINT_VAR(MESSAGE_WRAPPER_UNIT_TEST_MODE)
  IF (MESSAGE_WRAPPER_UNIT_TEST_MODE)
    GLOBAL_SET(MESSAGE_WRAPPER_INPUT ${ARGN})
  ELSE()
    GLOBAL_SET(MESSAGE_WRAPPER_INPUT "")
    MESSAGE(${ARGN})
  ENDIF()
ENDFUNCTION()

