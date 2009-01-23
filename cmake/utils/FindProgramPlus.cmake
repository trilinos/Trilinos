INCLUDE(ParseVariableArguments)
INCLUDE(PrintNonemptyVar)

#
# Finds the absolute path for a program given optionally just the program name
#
#  Input keyword arguments:
#
#    NAMES name1 name2 ...
#     List of default names that will be search
#
#    PATHS path1 path2 ...
#      List of paths
#
#    DOC docStr
#      Documentation string
#

FUNCTION(FIND_PROGRAM_PLUS PROG_VAR)

  PARSE_ARGUMENTS(
    PARSE
    "NAMES;PATHS;DOC"
    ""
    ${ARGN}
    )

  PRINT_NONEMPTY_VAR(${PROG_VAR})

  IF (IS_ABSOLUTE ${PROG_VAR})
    #MESSAGE(STATUS "Is Absoute")
    SET(NAMES_ARGS ${PARSE_NAMES})
  ELSE()
    #MESSAGE(STATUS "Is Not Absolute")
    SET(NAMES_ARGS ${${PROG_VAR}} ${PARSE_NAMES}) 
    SET(${PROG_VAR} "${PROG_VAR}-NOTFOUND" CACHE FILEPATH "" FORCE)
  ENDIF()
  #PRINT_VAR(NAMES_ARGS)

  SET(DOC "${PARSE_DOC}  Can be full path or just exec name.")

  # Look for program in given paths first!
  FIND_PROGRAM( ${PROG_VAR}
    NAMES ${NAMES_ARGS}
    PATHS ${PARSE_PATHS}
    DOC ${DOC}
    NO_DEFAULT_PATH
    )
  FIND_PROGRAM( ${PROG_VAR}
    NAMES ${NAMES_ARGS}
    DOC ${DOC}
    )
  MARK_AS_ADVANCED(${PROG_VAR})

  PRINT_NONEMPTY_VAR(${PROG_VAR})

ENDFUNCTION()
