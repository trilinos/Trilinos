#
# Usage:
#
#   SET_CACHE_VAR(<varName> <varVal> <type>)
#
FUNCTION(SET_CACHE_VAR  VAR_NAME  VAR_VALUE  TYPE_STR)
  IF ("${${VAR_NAME}}" STREQUAL "")
    MESSAGE("-- " "Setting ${VAR_NAME}='${VAR_VALUE}' by default")
    SET(${VAR_NAME} "${VAR_VALUE}" CACHE  ${TYPE_STR}
      "Set in ${CMAKE_CURRENT_LIST_FILE}")
  ENDIF()
ENDFUNCTION()

#
# Usage:
#
#   SET_BOOL_CACHE_VAR(<varName> <varVal>)
#
FUNCTION(SET_BOOL_CACHE_VAR  VAR_NAME  VAR_VALUE)
  SET_CACHE_VAR(${VAR_NAME}  ${VAR_VALUE}  BOOL)
ENDFUNCTION()


function(assert_defined VARS)
  foreach(VAR ${VARS})
    if(NOT DEFINED ${VAR})
      message(SEND_ERROR "Error, the variable ${VAR} is not defined!")
    endif()
  endforeach()
endfunction()
