INCLUDE(MultilineSet)
INCLUDE(GlobalSet)
INCLUDE(AppendSet)
INCLUDE(AssertDefined)
INCLUDE(SetNotFound)
INCLUDE(DualScopeSet)
INCLUDE(ParseVariableArguments)


#
# Function that sets up cache variables for users to specify where to
# find a TPL's headers and libraries.
#
# This function only sets global varibles as a way to return state
# so it can be called from anywhere in the call stack.
#
# This function can set up a with header files and/or libraries.
#
# The following cache variables defined that are intended for the user
# to set:
#
#   ${TPL_NAME}_INCLUDE_DIRS:  A list of common-separated directory paths
#       that will be searched ...
#
# ToDO: Finish this documentation.
#

FUNCTION(TPL_DECLARE_LIBRARIES TPL_NAME)

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

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("TPL_DECLARE_LIBRARIES: ${TPL_NAME}")
    PRINT_VAR(PARSE_REQUIRED_HEADERS)
    PRINT_VAR(PARSE_REQUIRED_LIBS_NAMES)
  ENDIF()

  #
  # User options
  #

  IF (PARSE_REQUIRED_LIBS_NAMES)

    # Library directories
  
    MULTILINE_SET(DOCSTR
      "List of semi-colon separated paths to look for the TPL ${TPL_NAME}"
      " libraries.  This list of paths will be passed into a FIND_LIBRARY(...)"
      " command to find the libraries listed in ${TPL_NAME}_LIBRARY_NAMES."
      "  Note that this set of paths is also the default value used for"
      " ${TPL_NAME}_LIBRARY_DIRS.  Therefore, if the headers exist in the"
      " same directories as the library, you do not need to set"
      " ${TPL_NAME}_LIBRARY_DIRS."
      )
    ADVANCED_SET(${TPL_NAME}_LIBRARY_DIRS "" CACHE PATH ${DOCSTR})

    # Libraries
  
    MULTILINE_SET(DOCSTR
      "List of semi-colon separated names of libraries needed to link to for"
      " the TPL ${TPL_NAME}.  This list of libraries will be search for in"
      " FIND_LIBRARY(...) calls along with the directories specified with"
      " ${TPL_NAME}_LIBRARY_DIRS.  NOTE: This is not the final list of libraries"
      " used for linking.  That is specified by TPL_${TPL_NAME}_LIBRARIES!"
      )
    ADVANCED_SET(${TPL_NAME}_LIBRARY_NAMES ${PARSE_REQUIRED_LIBS_NAMES} 
      CACHE STRING ${DOCSTR})

    # Let the user override what the names of the libraries which might
    # actually mean that no libraies are searched for.
    SET(PARSE_REQUIRED_LIBS_NAMES ${${TPL_NAME}_LIBRARY_NAMES})

  ELSE()

    SET(${TPL_NAME}_LIBRARY_DIRS) # Just to ignore below!
  
  ENDIF()

  # Include directories

  IF (PARSE_REQUIRED_HEADERS)

    MULTILINE_SET(DOCSTR
      "List of semi-colon separated paths to look for the TPL ${TPL_NAME}"
      " headers.  This list of paths will be passed into a FIND_PATH(...)"
      " command to find the headers for ${TPL_NAME} (which are known in advance)."
      )
    ADVANCED_SET(${TPL_NAME}_INCLUDE_DIRS ${${TPL_NAME}_LIBRARY_DIRS}
      CACHE PATH ${DOCSTR})
  
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(${TPL_NAME}_LIBRARY_DIRS)
      PRINT_VAR(${TPL_NAME}_LIBRARY_NAMES)
      PRINT_VAR(${TPL_NAME}_INCLUDE_DIRS)
    ENDIF()

  ENDIF()

  #
  # Direct build options
  #

  IF (PARSE_REQUIRED_LIBS_NAMES)

    # Libraries
  
    IF (NOT TPL_${TPL_NAME}_LIBRARIES)

      IF (${TPL_NAME}_LIBRARY_DIRS)
        MESSAGE(STATUS "  ${TPL_NAME}_LIBRARY_DIRS='${${TPL_NAME}_LIBRARY_DIRS}'")
      ENDIF()
  
      SET(LIBRARIES_FOUND)
  
      FOREACH(LIBNAME_SET ${${TPL_NAME}_LIBRARY_NAMES})
  
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          PRINT_VAR(LIBNAME_SET)
        ENDIF()
  
        SET(LIBNAME_LIST ${LIBNAME_SET})
        SEPARATE_ARGUMENTS(LIBNAME_LIST)
  
        SET(LIBNAME_SET_LIB)
  
        FOREACH(LIBNAME ${LIBNAME_LIST})

          MESSAGE(STATUS "  Searching for library '${LIBNAME}' ...")
    
          IF (${TPL_NAME}_LIBRARY_DIRS)
            SET(PATHS_ARG PATHS ${${TPL_NAME}_LIBRARY_DIRS})
          ELSE()
            SET(PATHS_ARG PATHS)
          ENDIF()
    
          SET_NOTFOUND(_${TPL_NAME}_${LIBNAME}_LIBRARY)
          FIND_LIBRARY( _${TPL_NAME}_${LIBNAME}_LIBRARY
            NAMES ${LIBNAME}
            ${PATHS_ARG} NO_DEFAULT_PATH )
          FIND_LIBRARY( _${TPL_NAME}_${LIBNAME}_LIBRARY
            NAMES ${LIBNAME} )
          MARK_AS_ADVANCED(_${TPL_NAME}_${LIBNAME}_LIBRARY)
  
          IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
            PRINT_VAR(_${TPL_NAME}_${LIBNAME}_LIBRARY)
          ENDIF()
    
          IF (_${TPL_NAME}_${LIBNAME}_LIBRARY)
            MESSAGE(STATUS "    Found ${TPL_NAME} TPL library: ${_${TPL_NAME}_${LIBNAME}_LIBRARY}")
            SET(LIBNAME_SET_LIB ${_${TPL_NAME}_${LIBNAME}_LIBRARY})
            BREAK()
          ENDIF()
    
        ENDFOREACH()
  
        IF (NOT LIBNAME_SET_LIB)
          MULTILINE_SET(ERRMSG
            "Warning: Could not find a library in the set \"${LIBNAME_SET}\" for"
            " the TPL ${TPL_NAME}!  Please manually set"
            " ${TPL_NAME}_LIBRARY_DIRS and/or ${TPL_NAME}_LIBRARY_NAMES or just"
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
        " This variable, however, is the final value and will not be touched."
        )
      ADVANCED_SET( TPL_${TPL_NAME}_LIBRARIES ${LIBRARIES_FOUND}
        CACHE FILEPATH ${DOCSTR} )
    
      IF (NOT TPL_${TPL_NAME}_LIBRARIES)
        MULTILINE_SET(ERRMSG
          "Error, could not find the ${TPL_NAME} Library!  Please manually set"
          " ${TPL_NAME}_LIBRARY_DIRS and/or ${TPL_NAME}_LIBRARY_NAMES or just"
          " TPL_${TPL_NAME}_LIBRARIES to point to the ${TPL_NAME} libraries!")
        MESSAGE(FATAL_ERROR ${ERRMSG})
      ENDIF()
  
    ENDIF()
  
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(TPL_${TPL_NAME}_LIBRARIES)
    ENDIF()

  ELSE()
  
    # There are no libraries so set the libraries to null but don't
    # change the cache which should not even have this varaible in it.
    # This set command is only to follow the standards for the package
    # support CMake code.
    GLOBAL_NULL_SET(TPL_${TPL_NAME}_LIBRARIES)
        
  ENDIF()

  # Include directories
  
  IF (PARSE_REQUIRED_HEADERS)

    IF (NOT TPL_${TPL_NAME}_INCLUDE_DIRS)
    
      FOREACH(INCLUDE_FILE_SET ${PARSE_REQUIRED_HEADERS})
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          PRINT_VAR(INCLUDE_FILE_SET)
        ENDIF()
        
        SET(INCLUDE_FILE_LIST ${INCLUDE_FILE_SET})
        SEPARATE_ARGUMENTS(INCLUDE_FILE_LIST)
        
        SET(INCLUDE_PATHS)
        
        FOREACH(INCLUDE_FILE ${INCLUDE_FILE_LIST})
          IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
            PRINT_VAR(INCLUDE_FILE)
          ENDIF()
          
          SET_NOTFOUND(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
          FIND_PATH( _${TPL_NAME}_${INCLUDE_FILE}_PATH
            NAMES ${INCLUDE_FILE}
            PATHS ${${TPL_NAME}_INCLUDE_DIRS}
            NO_DEFAULT_PATH)
          FIND_PATH( _${TPL_NAME}_${INCLUDE_FILE}_PATH
            NAMES ${INCLUDE_FILE} )
          MARK_AS_ADVANCED(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
          
          IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
            PRINT_VAR(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
          ENDIF()
          
          IF(_${TPL_NAME}_${INCLUDE_FILE}_PATH)
            MESSAGE(STATUS "  Found ${TPL_NAME} TPL header: ${_${TPL_NAME}_${INCLUDE_FILE}_PATH}")
            SET(INCLUDE_FILE_SET_PATH ${_${TPL_NAME}_${INCLUDE_FILE}_PATH})
            BREAK()
          ENDIF()
        ENDFOREACH()
        
        IF(NOT INCLUDE_FILE_SET_PATH)
          MULTILINE_SET(ERRMSG
            "Warning: Could not find a header in the set \"${INCLUDE_FILE_SET}\" for"
            " the TPL ${TPL_NAME}!  Please manually set"
            " ${TPL_NAME}_INCLUDE_DIRS and or just"
            " TPL_${TPL_NAME}_INCLUDE_DIRS to point to the ${TPL_NAME} includes!")
        ENDIF()
        
        APPEND_SET(INCLUDES_FOUND ${INCLUDE_FILE_SET_PATH})

      ENDFOREACH()
    
      MULTILINE_SET(DOCSTR
        "List of semi-colon separated paths to append to the compile invocations"
        " to find the headers for the TPL ${TPL_NAME}.  This is the final variable"
        " that is used in the build commands.  The user variable ${TPL_NAME}_INCLUDE_DIRS"
        " is used to look for the given headers first but is just a suggestion."
        " This variable, however, is the final value and will not be touched."
        )
      
      ADVANCED_SET(TPL_${TPL_NAME}_INCLUDE_DIRS ${INCLUDES_FOUND}
        CACHE PATH ${DOCSTR})
  
      IF (NOT TPL_${TPL_NAME}_INCLUDE_DIRS)
        MULTILINE_SET(ERRMSG
          "Error, could not find the ${TPL_NAME} headers include directory!"
          " Please manually set ${TPL_NAME}_INCLUDE_DIRS and/or"
          " ${TPL_NAME}_LIBRARY_DIRS or TPL_${TPL_NAME}_INCLUDE_DIRS to point"
          " to the ${TPL_NAME} headers!")
        MESSAGE(FATAL_ERROR ${ERRMSG})
      ENDIF()
    
      IF (TPL_${TPL_NAME}_INCLUDE_DIRS)
        MESSAGE(STATUS "  Found ${TPL_NAME} TPL header path: ${TPL_${TPL_NAME}_INCLUDE_DIRS}")
      ENDIF()
    ELSE()

      # TPL_${TPL_NAME}_INCLUDE_DIRS is already in the cache so leave it alone!

    ENDIF()

  ELSE()

    # Library has no header files so just set them to null
    GLOBAL_NULL_SET(TPL_${TPL_NAME}_INCLUDE_DIRS)

  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(TPL_${TPL_NAME}_INCLUDE_DIRS)
  ENDIF()

  # Set library directories to null always.  We do this because
  # the package support code expects this variable and it is used
  # for package dependencies.  Therefore, we need it to allow
  # TPLs and internal packages to be treated in the same way.
  GLOBAL_NULL_SET(TPL_${TPL_NAME}_LIBRARY_DIRS)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(TPL_${TPL_NAME}_LIBRARY_DIRS)
  ENDIF()

ENDFUNCTION()
