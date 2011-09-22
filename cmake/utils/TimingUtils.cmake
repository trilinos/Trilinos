#
# Timing utilities
#
# Warning: Depends on calling 'date' program so will not be portable to all
# platforms so call with care.


#
# Return the raw time in seconds since epoch, i.e., since 1970-01-01 00:00:00
# UTC
#
FUNCTION(TIMER_GET_RAW_SECONDS   SECONDS_RAW_OUT)
  EXECUTE_PROCESS(COMMAND date "+%s" OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE SECONDS_RAW)
  SET(${SECONDS_RAW_OUT} ${SECONDS_RAW} PARENT_SCOPE)
ENDFUNCTION()


#
# Return the relative time between start and stop seconds
#
FUNCTION(TIMER_GET_REL_SECONDS  SECONDS_RAW_START
  SECONDS_RAW_END  SECONDS_REL_OUT
  )
  MATH(EXPR SECONDS_REL "${SECONDS_RAW_END} - ${SECONDS_RAW_START}")
  SET(${SECONDS_REL_OUT} ${SECONDS_REL} PARENT_SCOPE)
ENDFUNCTION()


#
# Print the relative time in minutes
#
FUNCTION(TIMER_PRINT_REL_TIME  SECONDS_RAW_START   SECONDS_RAW_END
  MESSAGE_STR
  )
  TIMER_GET_REL_SECONDS(${SECONDS_RAW_START}  ${SECONDS_RAW_END}
     SECONDS_REL)
  # CMake does not support floating point so I need to do this manually
  MATH(EXPR  MINUTES_REL  "${SECONDS_REL}/60")
  MATH(EXPR  SECONDS_REL_REMAINING  "${SECONDS_REL} - 60*${MINUTES_REL}")
  MESSAGE("${MESSAGE_STR}: ${MINUTES_REL}m${SECONDS_REL_REMAINING}s")
ENDFUNCTION()
