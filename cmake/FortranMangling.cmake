# Function that defines variables describing the Fortran name mangling
# convention
#
# Sets the following outputs on success:
#
#  FC_FN_CASE
#    "UPPER" if names are translated to upper-case, 
#    "LOWER" otherwise.
#
#  FC_FN_UNDERSCORE
#    "NO_UNDER" if nothing is appended, "UNDER" if
#    one underscore is appended, and "SECOND_UNDER"
#    if a function with an underscore has a second
#    appended.
#
# FC_FUNC_DEFAULT
#    The default mange mangling for Fortran functions
#    that do not contain an underscore.
#
# FC_FUNC__DEFAULT
#    The default mange mangling for Fortran functions
#    that do contain an underscore.
#
#  The Fortran 2003 name binding facilities and ISO_C_BINDING module
#  should be preferred over cpp macro trickery whenever possible.
#
FUNCTION(FORTRAN_MANGLING)

  IF(NOT DEFINED FC_FN_CASE)

    IF (VERBOSE_CONFIGURE)
      MESSAGE("FORTRAN_MANGLING: Testing name Mangling Schemes!\n")
    ENDIF()

    FIND_FILE(_fcmakelists fmangle/ ${CMAKE_MODULE_PATH})
    IF (NOT _fcmakelists)
      MESSAGE(STATUS "Error, the file fmangle could not be found so we can not determine Fortran name mangling!")
      RETURN()
    ENDIF()

    SET(_fcmangledir ${PROJECT_BINARY_DIR}/CMakeFiles/CMakeTmp/fmangle)
    FILE(MAKE_DIRECTORY ${_fcmangledir})

    FOREACH(cdef LOWER UPPER)

      FOREACH(udef UNDER NO_UNDER SECOND_UNDER)

        IF (VERBOSE_CONFIGURE)
          MESSAGE("FORTRAN_MANGLING: Testing ${cdef} ${udef}\n\n")
        ENDIF()

        SET(_fcmangledir_case "${_fcmangledir}/${cdef}/${udef}")
        FILE(MAKE_DIRECTORY "${_fcmangledir}/${cdef}")
        FILE(MAKE_DIRECTORY ${_fcmangledir_case})

        TRY_COMPILE(_fcmngl ${_fcmangledir_case} ${_fcmakelists} fmangle
          CMAKE_FLAGS "-DMANGLE_FLAGS:STRING=-DFC_FN_${cdef};-DFC_FN_${udef}"
          OUTPUT_VARIABLE _fcmngl_output
          )
        IF (VERBOSE_CONFIGURE)
          MESSAGE("${_fcmngl_output}\n\n")
        ENDIF()

        IF(_fcmngl)
          IF (VERBOSE_CONFIGURE)
            MESSAGE("FORTRAN_MANGLING: Bingo!  ${cdef} ${udef} is the correct fortran name mangling!\n")
          ENDIF()
          SET(FC_FN_CASE ${cdef} CACHE INTERNAL
            "Case used by Fortran functions" FORCE)
          SET(FC_FN_UNDERSCORE ${udef} CACHE INTERNAL
            "Underscore convention used by Fortran functions" FORCE)
          BREAK()
        ENDIF()

      ENDFOREACH()

      IF(_fcmngl)
        BREAK()
      ENDIF()

    ENDFOREACH()

    IF(_fcmngl)
      MESSAGE(STATUS "Fortran name mangling: ${FC_FN_CASE} ${FC_FN_UNDERSCORE}")
    ELSE()
      MESSAGE(STATUS "Warning, cannot automatically determine Fortran mangling.")
    ENDIF()

  ENDIF()

  IF (FC_FN_CASE STREQUAL LOWER)
    SET(FC_NAME_NAME name)
  ELSEIF (FC_FN_CASE STREQUAL UPPER)
    SET(FC_NAME_NAME NAME)
  ENDIF()

  IF (FC_FN_UNDERSCORE)
    IF(FC_FN_UNDERSCORE STREQUAL "UNDER")
      SET(FC_FUNC_DEFAULT "(name,NAME) ${FC_NAME_NAME} ## _" CACHE INTERNAL "")
      SET(FC_FUNC__DEFAULT "(name,NAME) ${FC_NAME_NAME} ## _" CACHE INTERNAL "")
    ELSEIF(FC_FN_UNDERSCORE STREQUAL "SECOND_UNDER")
      SET(FC_FUNC_DEFAULT "(name,NAME) ${FC_NAME_NAME} ## _" CACHE INTERNAL "")
      SET(FC_FUNC__DEFAULT "(name,NAME) ${FC_NAME_NAME} ## __" CACHE INTERNAL "")
    ELSE()
      SET(FC_FUNC_DEFAULT "(name,NAME) ${FC_NAME_NAME}" CACHE INTERNAL "")
      SET(FC_FUNC__DEFAULT "(name,NAME) ${FC_NAME_NAME}" CACHE INTERNAL "")
    ENDIF()
  ENDIF()

ENDFUNCTION()


# 2008/10/26: rabartl: Below, these were macros that were also present in the
# file that I got off of the web but I did not need them so they are disbled
# (for now)

## - Guess if the Fortran compiler returns REAL in C doubles.
##  FORTRAN_FLOATRET()
##
## If the Fortran compiler follows the f2c convention of
## returning REALs in C doubles, routines like SLAMCH can
## cause very difficult to find stack corruptions.  This
## test is not perfect; it just checks if C->Fortran->C
## twice in a row does not crash.  If the tests do crash,
## FC_FN_FLOATRET is set to a true value, otherwise it
## is set to a false value.
##
## The REAL kinds in Fortran 2003's ISO_C_BINDING module should
## be used instead of this test whenever possible.
##
#MACRO(FORTRAN_FLOATRET)
#  IF(NOT DEFINED FC_FN_FLOATRET)
#    # Find inputs
#    FIND_FILE(_fcindir floatret/ ${CMAKE_MODULE_PATH})
#    SET(_fcdir ${PROJECT_BINARY_DIR}/CMakeFiles/CMakeTmp/floatret)
#    FILE(MAKE_DIRECTORY ${_fcdir})
#    TRY_COMPILE(_fccrv ${_fcdir} ${_fcindir} floatret
#      CMAKE_FLAGS "-DFC_FN_DEFS:STRING=${FC_FN_DEFS}")
#    IF(_fccrv)
#      EXECUTE_PROCESS(COMMAND ${_fcdir}/ctst
#        WORKING_DIRECTORY ${_fcdir}
#        RESULT_VARIABLE _fcrrv)
#    ENDIF(_fccrv)
#    IF(_fcrrv EQUAL 0)
#      SET(_fc_fn_floatret 0)
#    ELSE(_fcrrv EQUAL 0)
#      SET(_fc_fn_floatret 1)
#    ENDIF(_fcrrv EQUAL 0)
#    SET(FC_FN_FLOATRET ${_fc_fn_floatret} CACHE BOOL
#      "Fortran returns REAL in double.")
#    MESSAGE(STATUS "Fortran returns REAL in double: ${FC_FN_FLOATRET}")
#  ENDIF(NOT DEFINED FC_FN_FLOATRET)
#ENDMACRO()
#
#
## - Guess the convention for passing strings from C to Fortran.
##  FORTRAN_STRINGARG()
##
## If string lengths are directly appended to each variable, e.g.
## CALL FOO('bar', 1.0) becomes foo({'b','a','r'}, 3, 1.0), then
## FC_FN_STRINGARG is set to PAIRED.  If the lengths are appended
## to the call, e.g. foo({'b','a','r'}, 1.0, 3), FC_FN_STRINGARG
## is set to TRAILING.
## 
## This macro does not currently check for older Cray and VMS
## conventions that require conversion functions.  It also assumes
## that the length is passed as the "natural" size type, C's size_t.
##
## The string kinds in Fortran 2003's ISO_C_BINDING module should
## be used instead of these conventions whenever possible.
##
#SET(FC_FN_STRINGARG_TYPES "PAIRED" "TRAILING")
#MACRO(FORTRAN_STRINGARG)
#  IF(NOT DEFINED FC_FN_STRINGARG)
#    # Find inputs
#    FIND_FILE(_fcstrindir fstrings/ ${CMAKE_MODULE_PATH})
#    SET(_fcstrdir ${PROJECT_BINARY_DIR}/CMakeFiles/CMakeTmp/fstrings)
#    FOREACH(argtype ${FC_FN_STRINGARG_TYPES})
#      FILE(MAKE_DIRECTORY ${_fcstrdir})
#      TRY_COMPILE(_fcstrcrv ${_fcstrdir} ${_fcstrindir} fstrings
#        CMAKE_FLAGS "-DFC_FN_DEFS:STRING=${FC_FN_DEFS}" "-Dargtype:STRING=${argtype}")
#      IF(_fcstrcrv)
#        EXECUTE_PROCESS(COMMAND ${_fcstrdir}/ccheck
#          WORKING_DIRECTORY ${_fcstrdir}
#          RESULT_VARIABLE _fcstrrrv)
#        IF(_fcstrrrv EQUAL 0)
#          SET(_fcstr ${argtype})
#          BREAK()
#        ENDIF(_fcstrrrv EQUAL 0)
#      ENDIF(_fcstrcrv)
#      FILE(REMOVE_RECURSE ${_fcstrdir})
#    ENDFOREACH(argtype)
#    IF(DEFINED _fcstr)
#      SET(FC_FN_STRINGARG ${_fcstr} CACHE STRING
#        "How Fortran accepts string arguments.")
#      MESSAGE(STATUS "Fortran string passing: ${FC_FN_STRINGARG}")
#    ELSE(DEFINED _fcstr)
#      MESSAGE(STATUS "Cannot determine Fortran string passing.")
#    ENDIF(DEFINED _fcstr)
#  ENDIF(NOT DEFINED FC_FN_STRINGARG)
#ENDMACRO()
