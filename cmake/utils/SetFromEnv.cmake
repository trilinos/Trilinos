
INCLUDE(SetDefault)
INCLUDE(PrintVar)

#
# Set from the env if defined on the env.  Otherwise, just leave it alone.
#
MACRO(SET_FROM_ENV VAR)
  
  SET(ENV_${VAR} $ENV{${VAR}})
  IF (NOT "${ENV_${VAR}}" STREQUAL "")
    PRINT_VAR(ENV_${VAR})
    SET(${VAR} ${ENV_${VAR}})
  ENDIF()

  PRINT_VAR(${VAR})

ENDMACRO()
