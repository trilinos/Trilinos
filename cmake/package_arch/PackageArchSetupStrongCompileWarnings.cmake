INCLUDE(CheckCSourceCompiles)
INCLUDE(CheckCXXSourceCompiles)

INCLUDE(AdvancedSet)
INCLUDE(AssertDefined)
INCLUDE(MultilineSet)


FUNCTION(CHECK_C_COMPILER_FLAGS FLAGS VAR)
  IF (NOT ${VAR} STREQUAL OFF)
    SET(CMAKE_REQUIRED_FLAGS ${FLAGS})
    CHECK_C_SOURCE_COMPILES("int main() {return 0;}" ${VAR})
    MARK_AS_ADVANCED(${VAR})
  ENDIF()
ENDFUNCTION()


FUNCTION(CHECK_CXX_COMPILER_FLAGS FLAGS VAR)
  IF (NOT ${VAR} STREQUAL OFF)
    SET(CMAKE_REQUIRED_FLAGS ${FLAGS})
    CHECK_CXX_SOURCE_COMPILES("int main() {return 0;}"  ${VAR})
    MARK_AS_ADVANCED(${VAR})
  ENDIF()
ENDFUNCTION()

#
# Function that sets up strong compile options for the primary
# development platform (i.e. gcc)
#
# NOTE: The compiler flags in the cache, which may have been set by
# the user, are not disturbed in this function.  Instead variables in
# the parent base scope are set.
#

FUNCTION(PACKAGE_ARCH_SETUP_STRONG_COMPILE_WARNINGS)

  #
  # Setup and general flags
  #

  SET(CMAKE_BUILD_TYPES_LIST
    DEBUG MINSIZEREL RELEASE RELWITHDEBINFO)

  SET(GENERAL_DEBUG_FLAGS "-g -O0")

  SET(GENERAL_RELEASE_FLAGS "-O3")

  MULTILINE_SET(C_STRONG_COMPILE_WARNING_FLAGS
    " -ansi" # Check for C89 or C++98 standard code
    " -pedantic" # Adds more strick checking to remove non-ANSI GNU extensions
    " -Wall " # Do a bunch of default warnings (turns on a lot of the warnings above)
    " -Wwrite-strings" # Checks for non-const char * copy of string constants
    " -fexceptions" # Make sure that exceptions can be propogated through C code
    " -Wno-long-long" # Allow long long int since it is used by MPI, SWIG etc
    )
  
  MULTILINE_SET(CXX_STRONG_COMPILE_WARNING_FLAGS
    ${C_STRONG_COMPILE_WARNING_FLAGS}
    " -Wshadow" # Warn about general shadowing
    " -Woverloaded-virtual" # Warn about hiding virtual functions
    )

  ADVANCED_SET( ${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS "-Werror"
    CACHE STRING "Flags for treating warnings as errors.  To turn off set to ''")

  #
  # C compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_C)
  IF (${PROJECT_NAME}_ENABLE_C)
    
    IF (NOT CMAKE_C_FLAGS_DEBUG)
      CHECK_C_COMPILER_FLAGS(${GENERAL_DEBUG_FLAGS}
        ${PROJECT_NAME}_ENABLE_C_DEBUG_COMPILE_FLAGS )
      IF (${PROJECT_NAME}_ENABLE_C_DEBUG_COMPILE_FLAGS)
        SET(CMAKE_C_FLAGS_DEBUG ${GENERAL_DEBUG_FLAGS} PARENT_SCOPE)
      ENDIF()
    ENDIF()
    
    IF (NOT CMAKE_C_FLAGS_RELEASE)
      CHECK_C_COMPILER_FLAGS(${GENERAL_RELEASE_FLAGS}
        ${PROJECT_NAME}_ENABLE_C_RELEASE_COMPILE_FLAGS )
      IF (${PROJECT_NAME}_ENABLE_C_RELEASE_COMPILE_FLAGS)
        SET(CMAKE_C_FLAGS_RELEASE ${GENERAL_RELEASE_FLAGS} PARENT_SCOPE)
      ENDIF()
    ENDIF()

    CHECK_C_COMPILER_FLAGS( ${C_STRONG_COMPILE_WARNING_FLAGS}
      ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS )
    IF (${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS)
      FOREACH(BUILD_TYPE ${CMAKE_BUILD_TYPES_LIST})
        SET(CMAKE_C_FLAGS_${BUILD_TYPE}
         "${C_STRONG_COMPILE_WARNING_FLAGS} ${CMAKE_C_FLAGS_${BUILD_TYPE}}"
          PARENT_SCOPE)
      ENDFOREACH()
    ENDIF()
  
  ENDIF()

  #
  # C++ compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX)
  IF (${PROJECT_NAME}_ENABLE_CXX)
    
    IF (NOT CMAKE_CXX_FLAGS_DEBUG)
      CHECK_CXX_COMPILER_FLAGS(${GENERAL_DEBUG_FLAGS}
        ${PROJECT_NAME}_ENABLE_CXX_DEBUG_COMPILE_FLAGS )
      IF (${PROJECT_NAME}_ENABLE_CXX_DEBUG_COMPILE_FLAGS)
        SET(CMAKE_CXX_FLAGS_DEBUG ${GENERAL_DEBUG_FLAGS} PARENT_SCOPE)
      ENDIF()
    ENDIF()
    
    IF (NOT CMAKE_CXX_FLAGS_RELEASE)
      CHECK_CXX_COMPILER_FLAGS(${GENERAL_RELEASE_FLAGS}
        ${PROJECT_NAME}_ENABLE_CXX_RELEASE_COMPILE_FLAGS )
      IF (${PROJECT_NAME}_ENABLE_CXX_RELEASE_COMPILE_FLAGS)
        SET(CMAKE_CXX_FLAGS_RELEASE ${GENERAL_RELEASE_FLAGS} PARENT_SCOPE)
      ENDIF()
    ENDIF()
    
    CHECK_CXX_COMPILER_FLAGS( ${CXX_STRONG_COMPILE_WARNING_FLAGS}
      ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS )
    IF (${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS)
      FOREACH(BUILD_TYPE ${CMAKE_BUILD_TYPES_LIST})
        SET(CMAKE_CXX_FLAGS_${BUILD_TYPE}
         "${CXX_STRONG_COMPILE_WARNING_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TYPE}}"
         PARENT_SCOPE)
      ENDFOREACH()
    ENDIF()

  ENDIF()
  
ENDFUNCTION()
