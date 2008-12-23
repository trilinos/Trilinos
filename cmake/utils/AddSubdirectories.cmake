#
# Macro that adds a list of subdirectories
#

MACRO(ADD_SUBDIRECTORIES)
  FOREACH(DIR ${ARGV})
    ADD_SUBDIRECTORY(${DIR})
  ENDFOREACH()
ENDMACRO()
