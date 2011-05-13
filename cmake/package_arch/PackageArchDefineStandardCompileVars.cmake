INCLUDE(CMakeBuildTypesList)

INCLUDE(AdvancedSet)
INCLUDE(AssertDefined)
INCLUDE(MultilineSet)
INCLUDE(DualScopeSet)
INCLUDE(PrependCmndlineArgs)


#
# Macro that just defines the basic flags
#

MACRO(PACKAGE_ARCH_DEFINE_STANDARD_COMPILE_FLAGS_VARS  ENABLE_SHADOWING_WARNINGS)

  #
  # Setup and general flags
  #

  SET(GENERAL_BUILD_FLAGS) # Applies to all builds, period
 
  IF (${PROJECT_NAME}_ENABLE_CHECKED_STL)
    PREPEND_CMNDLINE_ARGS(GENERAL_BUILD_FLAGS "-D_GLIBCXX_DEBUG")
  ENDIF()

  SET(DEBUG_SYMBOLS_FLAGS "-g")

  SET(GENERAL_DEBUG_FLAGS "${DEBUG_SYMBOLS_FLAGS} -O0")
 
  SET(GENERAL_RELEASE_FLAGS "-O3")

  IF (${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS)
    PREPEND_CMNDLINE_ARGS(GENERAL_BUILD_FLAGS "${DEBUG_SYMBOLS_FLAGS}")
  ENDIF()

  MULTILINE_SET(C_STRONG_COMPILE_WARNING_FLAGS
    "-ansi" # Check for C89 or C++98 standard code
    " -pedantic" # Adds more static checking to remove non-ANSI GNU extensions
    " -Wall" # Enable a bunch of default warnings
    " -Wno-long-long" # Allow long long int since it is used by MPI, SWIG, etc.
    )
  
  MULTILINE_SET(CXX_STRONG_COMPILE_WARNING_FLAGS
    ${C_STRONG_COMPILE_WARNING_FLAGS}
    " -Wwrite-strings" # Checks for non-const char * copy of string constants
    )

  IF (${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS)
    SET(LOCAL_ENABLE_SHADOWING_WARNINGS ON)
  ELSEIF (${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS STREQUAL "OFF")
    SET(LOCAL_ENABLE_SHADOWING_WARNINGS OFF)
  ELSE()
    SET(LOCAL_ENABLE_SHADOWING_WARNINGS ${ENABLE_SHADOWING_WARNINGS})
  ENDIF()
  
  IF (LOCAL_ENABLE_SHADOWING_WARNINGS)
    MULTILINE_SET(CXX_STRONG_COMPILE_WARNING_FLAGS
      ${CXX_STRONG_COMPILE_WARNING_FLAGS}
      " -Wshadow" # Warn about general shadowing issues
      " -Woverloaded-virtual" # Warn about hiding virtual functions
      )
  ENDIF()

ENDMACRO()
