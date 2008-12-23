# Function that defines variables describing the BLAS name mangling convention
#
# Sets the following outputs on success:
#
#  BLAS_FN_CASE
#    "UPPER" if names are translated to upper-case, 
#    "LOWER" otherwise.
#
#  BLAS_FN_UNDERSCORE
#    "NO_UNDER" if nothing is appended, "UNDER" if
#    one underscore is appended, and "SECOND_UNDER"
#    if a function with an underscore has a second
#    appended.
#
# BLAS_FUNC_DEFAULT
#    The default mange mangling for Fortran functions
#    that do not contain an underscore.
#
#  The Fortran 2003 name binding facilities and ISO_C_BINDING module
#  should be preferred over cpp macro trickery whenever possible.
#
FUNCTION(BLAS_MANGLING)

  IF(NOT DEFINED BLAS_FN_CASE)

    IF (Trilinos_VERBOSE_CONFIGURE)
      MESSAGE("BLAS_MANGLING: Testing name Mangling Schemes!\n")
    ENDIF()

    FIND_FILE(_blascmakelists blasmangle/ ${CMAKE_MODULE_PATH})
    IF (NOT _blascmakelists)
      MESSAGE(STATUS "Error, the file blasmangle could not be found so we can not determine Fortran name mangling!")
      RETURN()
    ENDIF()

    SET(_fcmangledir ${PROJECT_BINARY_DIR}/CMakeFiles/CMakeTmp/blasmangle)
    FILE(MAKE_DIRECTORY ${_fcmangledir})

    FOREACH(cdef LOWER UPPER)

      FOREACH(udef UNDER NO_UNDER)

        IF (Trilinos_VERBOSE_CONFIGURE)
          MESSAGE("BLAS_MANGLING: Testing ${cdef} ${udef}\n\n")
        ENDIF()

        SET(_fcmangledir_case "${_fcmangledir}/${cdef}/${udef}")
        FILE(MAKE_DIRECTORY "${_fcmangledir}/${cdef}")
        FILE(MAKE_DIRECTORY ${_fcmangledir_case})

        TRY_COMPILE(_blasmngl ${_fcmangledir_case} ${_blascmakelists} blasmangle
          CMAKE_FLAGS "-DMANGLE_FLAGS:STRING=-DFC_FN_${cdef};-DFC_FN_${udef}"
          OUTPUT_VARIABLE _blasmngl_output
          )
        IF (Trilinos_VERBOSE_CONFIGURE)
          MESSAGE("${_blasmngl_output}\n\n")
        ENDIF()

        IF(_blasmngl)
          IF (Trilinos_VERBOSE_CONFIGURE)
            MESSAGE("BLAS_MANGLING: Bingo!  ${cdef} ${udef} is the correct BLAS name mangling!\n")
          ENDIF()
          SET(BLAS_FN_CASE ${cdef} CACHE INTERNAL
            "Case used by Fortran functions" FORCE)
          SET(BLAS_FN_UNDERSCORE ${udef} CACHE INTERNAL
            "Underscore convention used by Fortran functions" FORCE)
          BREAK()
        ENDIF()

      ENDFOREACH()

      IF(_blasmngl)
        BREAK()
      ENDIF()

    ENDFOREACH()

    IF(_blasmngl)
      MESSAGE(STATUS "BLAS name mangling: ${BLAS_FN_CASE} ${BLAS_FN_UNDERSCORE}")
    ELSE()
      MESSAGE(STATUS "Warning, cannot automatically determine BLAS mangling.")
    ENDIF()

  ENDIF()

  IF (BLAS_FN_CASE STREQUAL LOWER)
    SET(BLAS_NAME_NAME name)
  ELSEIF (BLAS_FN_CASE STREQUAL UPPER)
    SET(BLAS_NAME_NAME NAME)
  ENDIF()

  IF (BLAS_FN_UNDERSCORE)
    IF(BLAS_FN_UNDERSCORE STREQUAL "UNDER")
      SET(BLAS_FUNC_DEFAULT "(name,NAME) ${BLAS_NAME_NAME} ## _" CACHE INTERNAL "")
    ELSE()
      SET(BLAS_FUNC_DEFAULT "(name,NAME) ${BLAS_NAME_NAME}" CACHE INTERNAL "")
    ENDIF()
  ENDIF()

ENDFUNCTION()
