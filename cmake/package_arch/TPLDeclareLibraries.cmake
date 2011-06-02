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
# This function can set up a with header files and/or libraries.
#
# The input arguments to this function are:
#
#   REQUIRED_HEADERS:  List of header files that are searched for the TPL
#     using FIND_PATH(...).
#
#   MUST_FIND_ALL_HEADERS:  If set, then all of the header files listed in
#     REQUIRED_HEADERS must be found in order for TPL_${TPL_NAME}_INCLUDE_DIRS
#     to be defined.
#
#   REQUIRED_LIBS_NAMES: List of libraries that are searched for when
#     looked for the TPLs libraries with FIND_LIBRARY(...).
#
#   MUST_FIND_ALL_LIBS:  If set, then all of the library files listed in
#     REQUIRED_LIBS_NAMES must be found or the TPL is considered not
#     found!
#
#   NO_PRINT_ENABLE_SUCCESS_FAIL: If set, then the final success/fail
#     will not be printed
#
# The input cmake cache variables that this funciton uses (if defined) are:
#
#   ${TPL_NAME}_INCLUDE_DIRS:PATH: List of paths to search first for header
#      files defined in REQUIRED_HEADERS.
#
#   ${TPL_NAME}_LIBRARY_DIRS:PATH: The list of directories to search first
#      for libraies defined in REQUIRED_LIBS_NAMES.
#
#   ${TPL_NAME}_LIBRARY_NAMES:STIRNG: List of library names to be looked for
#      instead of what is specified in REQUIRED_LIBS_NAMES.
#
# This function only sets global varibles as a way to return state so it can
# be called from anywhere in the call stack.  The following cache variables
# defined that are intended for the user to set and/or use:
#
#   TPL_${TPL_NAME}_INCLUDE_DIRS:  A list of common-separated full directory paths
#     that contain the TPLs headers.  If this varible is set before calling
#     this function, then no headers are searched for and this variable will
#     be assumed to have the correct list of header paths.
#
#   TPL_${TPL_NAME}_LIBRARIES:  A list of commons-seprated full library
#     names (output from FIND_LIBRARY(...)) for all of the libraries found
#     for the TPL.  IF this varible is set before calling this function,
#     no libraries are searched for and this varaible will be assumed to
#     have the correct list of libraries to link to.
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
     "MUST_FIND_ALL_LIBS;MUST_FIND_ALL_HEADERS;NO_PRINT_ENABLE_SUCCESS_FAIL"
     ${ARGN}
     )

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("TPL_DECLARE_LIBRARIES: ${TPL_NAME}")
    PRINT_VAR(PARSE_REQUIRED_HEADERS)
    PRINT_VAR(PARSE_REQUIRED_LIBS_NAMES)
    PRINT_VAR(TPL_${TPL_NAME}_INCLUDE_DIRS)
    PRINT_VAR(TPL_${TPL_NAME}_LIBRARIES)
  ENDIF()

  IF (TPL_TENTATIVE_ENABLE_${TPL_NAME})
    MESSAGE(STATUS "  Attempting to enable tentatively enabled TPL '${TPL_NAME}' ...")
    SET(ERROR_MSG_MODE STATUS)
  ELSE()
    SET(ERROR_MSG_MODE FATAL_ERROR)
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
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(${TPL_NAME}_LIBRARY_DIRS)
    ENDIF()

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

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(${TPL_NAME}_LIBRARY_NAMES)
      PRINT_VAR(PARSE_REQUIRED_LIBS_NAMES)
    ENDIF()

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
  # Set the lib extentions to find
  #

  # Save the default the first time through
  IF (NOT CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT)
   SET(TPL_CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT ${CMAKE_FIND_LIBRARY_SUFFIXES})
   #PRINT_VAR(TPL_CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT)
  ENDIF()

  #PRINT_VAR(TPL_FIND_SHARED_LIBS)
  #PRINT_VAR(CMAKE_FIND_LIBRARY_SUFFIXES)
  # Set libraries to find
  IF (TPL_FIND_SHARED_LIBS)
    # The default should be to find shared libs first
    SET(TPL_CMAKE_FIND_LIBRARY_SUFFIXES ${TPL_CMAKE_FIND_LIBRARY_SUFFIXES_DEFAULT})
  ELSE()
    if(WIN32)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a)
    else()
      SET(CMAKE_FIND_LIBRARY_SUFFIXES .a )
    endif()
  ENDIF()
  #PRINT_VAR(CMAKE_FIND_LIBRARY_SUFFIXES)

  #
  # Direct build options
  #

  SET(_${TPL_NAME}_ENABLE_SUCCESS TRUE)

  IF (PARSE_REQUIRED_LIBS_NAMES)

    # Libraries
  
    IF (NOT TPL_${TPL_NAME}_LIBRARIES)

      IF (PARSE_MUST_FIND_ALL_LIBS)
        MESSAGE(STATUS "  Must find all libraries in \"${PARSE_REQUIRED_LIBS_NAMES}\"")
      ENDIF()

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
            "Could not find a library in the set \"${LIBNAME_SET}\" for"
            " the TPL ${TPL_NAME}!  Please manually set"
            " ${TPL_NAME}_LIBRARY_DIRS and/or ${TPL_NAME}_LIBRARY_NAMES or just"
            " TPL_${TPL_NAME}_LIBRARIES to point to the ${TPL_NAME} libraries!")
          IF (PARSE_MUST_FIND_ALL_LIBS)
            SET(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
            MESSAGE(${ERROR_MSG_MODE} "  Error: ${ERRMSG}")
          ELSE()
            MESSAGE(STATUS "  Warning: ${ERRMSG}")
          ENDIF()
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
          "Could not find the ${TPL_NAME} Library!  Please manually set"
          " ${TPL_NAME}_LIBRARY_DIRS and/or ${TPL_NAME}_LIBRARY_NAMES or just"
          " TPL_${TPL_NAME}_LIBRARIES to point to the ${TPL_NAME} libraries!")
        SET(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
        MESSAGE(${ERROR_MSG_MODE} ${ERRMSG})
        PRINT_VAR(_${TPL_NAME}_ENABLE_SUCCESS)
      ENDIF()
  
    ENDIF()
  
    # Print the final value to be used *always*
    MESSAGE(STATUS "  TPL_${TPL_NAME}_LIBRARIES='${TPL_${TPL_NAME}_LIBRARIES}'")

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

      IF (PARSE_MUST_FIND_ALL_HEADERS)
        MESSAGE(STATUS "  Must find all headers in \"${PARSE_REQUIRED_HEADERS}\"")
      ENDIF()
    
      FOREACH(INCLUDE_FILE_SET ${PARSE_REQUIRED_HEADERS})

        IF (PARSE_MUST_FIND_ALL_HEADERS)
          MESSAGE(STATUS "  Searching for headers \"${INCLUDE_FILE_SET}\"")
        ENDIF()
        
        SET(INCLUDE_FILE_LIST ${INCLUDE_FILE_SET})
        SEPARATE_ARGUMENTS(INCLUDE_FILE_LIST)
        SET(INCLUDE_FILE_SET_PATH) # Start out as empty list
        
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
            MESSAGE(STATUS "    Found ${TPL_NAME} TPL header: ${_${TPL_NAME}_${INCLUDE_FILE}_PATH}/${INCLUDE_FILE}")
            APPEND_SET(INCLUDE_FILE_SET_PATH ${_${TPL_NAME}_${INCLUDE_FILE}_PATH})
            IF(NOT PARSE_MUST_FIND_ALL_HEADERS)
              BREAK()
            ENDIF()
          ELSE()
            SET(USERMSG "    Did not find ${TPL_NAME} TPL header: ${INCLUDE_FILE}")
            IF(PARSE_MUST_FIND_ALL_HEADERS)
              SET(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
              MESSAGE(${ERROR_MSG_MODE} ${USERMSG})
              BREAK()
            ELSE()
              MESSAGE(STATUS ${USERMSG})
            ENDIF()
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

      ENDFOREACH(INCLUDE_FILE_SET ${PARSE_REQUIRED_HEADERS})
    
      IF (INCLUDES_FOUND)
        LIST(REMOVE_DUPLICATES INCLUDES_FOUND)
      ENDIF()

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
          "Could not find the ${TPL_NAME} headers include directory!"
          " Please manually set ${TPL_NAME}_INCLUDE_DIRS and/or"
          " ${TPL_NAME}_LIBRARY_DIRS or TPL_${TPL_NAME}_INCLUDE_DIRS to point"
          " to the ${TPL_NAME} headers!")
        SET(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
        MESSAGE(${ERROR_MSG_MODE} ${ERRMSG})
        PRINT_VAR(_${TPL_NAME}_ENABLE_SUCCESS)
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

  # Print the final value to be used *always*
  MESSAGE(STATUS "  TPL_${TPL_NAME}_INCLUDE_DIRS='${TPL_${TPL_NAME}_INCLUDE_DIRS}'")

  # Set library directories to null always.  We do this because
  # the package support code expects this variable and it is used
  # for package dependencies.  Therefore, we need it to allow
  # TPLs and internal packages to be treated in the same way.
  GLOBAL_NULL_SET(TPL_${TPL_NAME}_LIBRARY_DIRS)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(TPL_${TPL_NAME}_LIBRARY_DIRS)
  ENDIF()
  # 2011/05/09: rabartl: ToDo: Remove this above variable from everywhere!

  #PRINT_VAR(TPL_TENTATIVE_ENABLE_${TPL_NAME})
  #PRINT_VAR(_${TPL_NAME}_ENABLE_SUCCESS)
  IF (TPL_TENTATIVE_ENABLE_${TPL_NAME})
    IF (_${TPL_NAME}_ENABLE_SUCCESS)
      IF (NOT PARSE_NO_PRINT_ENABLE_SUCCESS_FAIL)
        MESSAGE(STATUS "  Attempt to enable tentatively enabled TPL '${TPL_NAME}' passed!")
      ENDIF()
    ELSE()
      IF (NOT PARSE_NO_PRINT_ENABLE_SUCCESS_FAIL)
        MESSAGE(STATUS "  Attempt to enable tentatively enabled TPL '${TPL_NAME}' failed!  Setting TPL_ENABLE_${TPL_NAME}=OFF")
      ENDIF()
      SET(TPL_ENABLE_${TPL_NAME} OFF CACHE STRING "autoset" FORCE)
    ENDIF()
  ENDIF()

ENDFUNCTION()


#
# Function that sets up for an optionally enabled TPL that is attempted to be
# enabled but will be disabled if all of the parts are not found.
#

FUNCTION(TPL_TENTATIVELY_ENABLE TPL_NAME)

  IF ("${TPL_ENABLE_${TPL_NAME}}" STREQUAL "")
    # The TPL's enable status has not been set so we will tentatively enable
    # it.
    SET(TPL_ENABLE_${TPL_NAME} ON CACHE STRING "autoset" FORCE)
    ADVANCED_SET(TPL_TENTATIVE_ENABLE_${TPL_NAME} ON CACHE STRING "autoset" FORCE)
  ELSE()
    # The TPL's enable status has already be hard set to be ON or OFF so we
    # will leave it alone.
  ENDIF()

ENDFUNCTION()
