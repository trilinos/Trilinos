INCLUDE(Multiline_Set)
INCLUDE(Global_Set)
INCLUDE(Append_Set)
INCLUDE(Assert_Defined)
INCLUDE(Parse_Variable_Arguments)


MACRO(TPL_DECLARE_LIBRARIES TPL_NAME)

  # Make sure the right name is used
  ASSERT_DEFINED(TPL_ENABLE_${TPL_NAME})

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "REQUIRED_HEADERS;REQUIRED_LIBS_NAMES"
     #options
     ""
     ${ARGN}
     )

  IF (Trilinos_VERBOSE_CONFIGURE)
    MESSAGE("TPL_DECLARE_LIBRARIES: ${TPL_NAME}")
    PRINT_VAR(PARSE_REQUIRED_HEADERS)
    PRINT_VAR(PARSE_REQUIRED_LIBS_NAMES)
  ENDIF()



  #
  # User options
  #

  # Library directories
  
  MULTILINE_SET(DOCSTR
    "List of semi-colon separated paths to look for the TPL ${TPL_NAME}"
    " libraries.  This list of paths will be passed into a FIND_LIBRARY(...)"
    " command to find the libraries listed in ${TPL_NAME}_LIBRARIES."
    "  Note that this set of paths is also the default value used for"
    " ${TPL_NAME}_LIBRARY_DIRS.  Therefore, if the headers exist in the"
    " same directories as the library, you do not need to set"
    " ${TPL_NAME}_LIBRARY_DIRS."
    )
  ADVANCED_SET(${TPL_NAME}_LIBRARY_DIRS "" CACHE STRING ${DOCSTR})

  # Libraries

  MULTILINE_SET(DOCSTR
    "List of semi-colon separated names of libraries needed to link to for"
    " the TPL ${TPL_NAME}.  This list of libraries will be search for in"
    " FIND_LIBRARY(...) calls along with the directories specified with"
    " ${TPL_NAME}_LIBRARY_DIRS.  NOTE: This is not the final list of libraries"
    " used for linking.  That is specified by TPL_${TPL_NAME}_LIBRARIES!"
    )
  ADVANCED_SET(${TPL_NAME}_LIBRARIES ${PARSE_REQUIRED_LIBS_NAMES} 
    CACHE STRING ${DOCSTR})

  # Include directories

  IF (PARSE_REQUIRED_HEADERS)

    MULTILINE_SET(DOCSTR
      "List of semi-colon separated paths to look for the TPL ${TPL_NAME}"
      " headers.  This list of paths will be passed into a FIND_PATH(...)"
      " command to find the headers for ${TPL_NAME} (which are known in advance)."
      )
    ADVANCED_SET(${TPL_NAME}_INCLUDE_DIRS ${${TPL_NAME}_LIBRARY_DIRS}
      CACHE STRING ${DOCSTR})
  
    IF (Trilinos_VERBOSE_CONFIGURE)
      PRINT_VAR(${TPL_NAME}_LIBRARY_DIRS)
      PRINT_VAR(${TPL_NAME}_LIBRARIES)
      PRINT_VAR(${TPL_NAME}_INCLUDE_DIRS)
    ENDIF()

  ENDIF()

  #
  # Direct build options
  #

  # Libraries

  IF (NOT TPL_${TPL_NAME}_LIBRARIES)

    SET(LIBRARIES_FOUND)

    FOREACH(LIBNAME_SET ${${TPL_NAME}_LIBRARIES})

      IF (Trilinos_VERBOSE_CONFIGURE)
        PRINT_VAR(LIBNAME_SET)
      ENDIF()

      SET(LIBNAME_LIST ${LIBNAME_SET})
      SEPARATE_ARGUMENTS(LIBNAME_LIST)

      SET(LIBNAME_SET_LIB)

      FOREACH(LIBNAME ${LIBNAME_LIST})

        IF (Trilinos_VERBOSE_CONFIGURE)
          PRINT_VAR(LIBNAME)
        ENDIF()
  
        IF (${TPL_NAME}_LIBRARY_DIRS)
          SET(PATHS_ARG PATHS ${${TPL_NAME}_LIBRARY_DIRS})
        ELSE()
          SET(PATHS_ARG PATHS)
        ENDIF()
  
        FIND_LIBRARY( _${TPL_NAME}_${LIBNAME}_LIBRARY
          NAMES ${LIBNAME}
          ${PATHS_ARG}
          )
        MARK_AS_ADVANCED(_${TPL_NAME}_${LIBNAME}_LIBRARY)

        IF (Trilinos_VERBOSE_CONFIGURE)
          PRINT_VAR(_${TPL_NAME}_${LIBNAME}_LIBRARY)
        ENDIF()
  
        IF (_${TPL_NAME}_${LIBNAME}_LIBRARY)
          MESSAGE(STATUS "  Found ${TPL_NAME} TPL library: ${_${TPL_NAME}_${LIBNAME}_LIBRARY}")
          SET(LIBNAME_SET_LIB ${_${TPL_NAME}_${LIBNAME}_LIBRARY})
          BREAK()
        ENDIF()
  
      ENDFOREACH()

      IF (NOT LIBNAME_SET_LIB)
        MULTILINE_SET(ERRMSG
          "Warning: Could not find a library in the set \"${LIBNAME_SET}\" for"
          " the TPL ${TPL_NAME}!  Please manually set"
          " ${TPL_NAME}_LIBRARY_DIRS and/or ${TPL_NAME}_LIBRARIES or just"
          " TPL_${TPL_NAME}_LIBRARIES to point to the ${TPL_NAME} libraries!")
        MESSAGE(STATUS ${ERRMSG})
      ENDIF()

      APPEND_SET(LIBRARIES_FOUND ${LIBNAME_SET_LIB})

    ENDFOREACH()

    MULTILINE_SET(DOCSTR
      "List of semi-colon separated full paths to the libraries for the TPL"
      " ${TPL_NAME}.  This is the final variable that is used in the link"
      " commands.  The user variable ${TPL_NAME}_LIBRARY_DIRS is used to look"
      " for the know library names but but is just a suggestion."
      " This varible, however, is the final value and will not be touched."
      )
    ADVANCED_SET( TPL_${TPL_NAME}_LIBRARIES ${LIBRARIES_FOUND}
      CACHE PATH ${DOCSTR} )
  
    IF (NOT TPL_${TPL_NAME}_LIBRARIES)
      MULTILINE_SET(ERRMSG
        "Error, could not find the ${TPL_NAME} Library!  Please manually set"
        " ${TPL_NAME}_LIBRARY_DIRS and/or ${TPL_NAME}_LIBRARIES or just"
        " TPL_${TPL_NAME}_LIBRARIES to point to the ${TPL_NAME} libraries!")
      MESSAGE(FATAL_ERROR ${ERRMSG})
    ENDIF()

  ENDIF()

  # Include directories

  IF (PARSE_REQUIRED_HEADERS AND NOT TPL_${TPL_NAME}_INCLUDE_DIRS)

    MULTILINE_SET(DOCSTR
      "List of semi-colon separated paths to append to the compile invocations"
      " to find the headers for the TPL ${TPL_NAME}.  This is the final variable"
      " that is used in the build commands.  The user variable ${TPL_NAME}_INCLUDE_DIRS"
      " is used to look for the given headers first but is just a suggestion."
      " This varible, however, is the final value and will not be touched."
      )

    FIND_PATH( TPL_${TPL_NAME}_INCLUDE_DIRS NAMES
      NAMES ${PARSE_REQUIRED_HEADERS}
      PATHS ${${TPL_NAME}_INCLUDE_DIRS}
      DOC ${DOCSTR} )
    MARK_AS_ADVANCED(TPL_${TPL_NAME}_LIBRARIES)
  
    IF (NOT TPL_${TPL_NAME}_INCLUDE_DIRS)
      MULTILINE_SET(ERRMSG
        "Error, could not find the ${TPL_NAME} headers include directory!"
        " Please manually set ${TPL_NAME}_INCUDE_DIRS and/or"
        " ${TPL_NAME}_LIBRARY_DIRS or TPL_${TPL_NAME}_INCLUDE_DIRS to point"
        " to the ${TPL_NAME} headers!")
      MESSAGE(FATAL_ERROR ${ERRMSG})
    ENDIF()

    # 2008/12/02: rabartl: ToDo: Above: Put in a check to see that all of the
    # headers that have been specified have indeed been found!

  ELSE()

   # Library has not header files
   GLOBAL_NULL_SET(TPL_${TPL_NAME}_INCLUDE_DIRS)

  ENDIF()

  # Library directories?

  GLOBAL_NULL_SET(TPL_${TPL_NAME}_LIBRARY_DIRS)
  # Not used for anything, just for consistency!

  IF (Trilinos_VERBOSE_CONFIGURE)
    PRINT_VAR(TPL_${TPL_NAME}_LIBRARIES)
    PRINT_VAR(TPL_${TPL_NAME}_INCLUDE_DIRS)
    PRINT_VAR(TPL_${TPL_NAME}_LIBRARY_DIRS)
  ENDIF()

ENDMACRO()
