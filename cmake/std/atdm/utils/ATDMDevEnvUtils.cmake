#
# ATDM_SET_CACHE(<varName> <defaultVal> CACHE <dataType>)
#
# This silently sets a cache var (but will be verbose if
# ATDM_CONFIG_VERBOSE_DEBUG=ON).  However, this set statement will only change
# anything if the cache var <varName> is not already set.

#
MACRO(ATDM_SET_CACHE VAR_NAME DEFAULT_VAL CACHE_STR DATA_TYPE)
  IF (ATDM_CONFIG_VERBOSE_DEBUG)
    MESSAGE("-- SET(${VAR_NAME} \"${DEFAULT_VAL}\" ${CACHE_STR} ${DATA_TYPE})")
  ENDIF()
  SET(${VAR_NAME} "${DEFAULT_VAL}"  ${CACHE_STR} ${DATA_TYPE}
    "Set in ${CMAKE_CURRENT_LIST_FILE}")
ENDMACRO()

#
# ATDM_SET_CACHE_FORCE(<varName> <defaultVal> CACHE <dataType>)
#
# This silently force sets a cache var (but will be verbose if
# ATDM_CONFIG_VERBOSE_DEBUG=ON.  This will override any existing value that is
# currently set in the cache.
#
MACRO(ATDM_SET_CACHE_FORCE VAR_NAME DEFAULT_VAL CACHE_STR DATA_TYPE)
  IF (ATDM_CONFIG_VERBOSE_DEBUG)
    MESSAGE("-- SET(${VAR_NAME} \"${DEFAULT_VAL}\" ${CACHE_STR} ${DATA_TYPE} <doc> FORCE)")
  ENDIF()
  SET(${VAR_NAME} "${DEFAULT_VAL}"  ${CACHE_STR} ${DATA_TYPE}
    "Set in ${CMAKE_CURRENT_LIST_FILE}"  FORCE)
ENDMACRO()

#
# ATDM_SET_ENABLE(<varName> <defaultVal>)
#
# This is most appropriate for enable vars that can have an unspecified empty
# state "" and we only wnat to set this default if has an empty state.
#
# This also prints a message about the default value being set but will always
# print the final value if `ATDM_CONFIG_VERBOSE_DEBUG` is set to `TRUE`.
#
MACRO(ATDM_SET_ENABLE  VAR_NAME  VAR_DEFAULT)
  IF ("${${VAR_NAME}}" STREQUAL "")
    MESSAGE("-- " "Setting default ${VAR_NAME}=${VAR_DEFAULT}")
    SET(${VAR_NAME} "${VAR_DEFAULT}" CACHE BOOL
      "Set in ${CMAKE_CURRENT_LIST_FILE}")
  ENDIF()
  IF (ATDM_CONFIG_VERBOSE_DEBUG)
    PRINT_VAR(${VAR_NAME})
  ENDIF()
ENDMACRO()

#
# ATDM_SET_ATDM_VAR_FROM_ENV_AND_DEFAULT(<varBaseName> <defaultVal>)
#
# Sets the var 'ATDM_${VAR_BASE_NAME}' from
# '$ENV{ATDM_CONFIG_${VAR_BASE_NAME}}' or if not set gives a default.
#
# This does **Not** set a cache var.  It only sets a local var with a devalue
# value that can be overrridden as an env var.
#
MACRO(ATDM_SET_ATDM_VAR_FROM_ENV_AND_DEFAULT VAR_BASE_NAME VAR_DEFAULT_VAL)
  SET(ATDM_${VAR_BASE_NAME} "$ENV{ATDM_CONFIG_${VAR_BASE_NAME}}")
  IF ("${ATDM_${VAR_BASE_NAME}}" STREQUAL "")
    SET(ATDM_${VAR_BASE_NAME} "${VAR_DEFAULT_VAL}")
  ENDIF()
  IF (ATDM_CONFIG_VERBOSE_DEBUG)
    PRINT_VAR(ATDM_${VAR_BASE_NAME})
  ENDIF()
ENDMACRO()
