##################################################################################
#
#                    Set defaults for XSDK CMake projects
#
##################################################################################

#
# This module implements standard behavior for XSDK CMake projects.  The main
# thing it does in XSDK mode (i.e. USE_XSDK_DEFAULTS=TRUE) is to print out
# when the env vars CC, CXX, FC and compiler flags CFLAGS, CXXFLAGS, and
# FFLAGS/FCFLAGS are used to select the compilers and compiler flags (raw
# CMake does this silently) and to set BUILD_SHARED_LIBS=TRUE and
# CMAKE_BUILD_TYPE=DEBUG by default.  It does not implement *all* of the
# standard XSDK configuration parameters.  The parent CMake project must do
# that.
#
# Note that when USE_XSDK_DEFAULTS=TRUE, then the Fortran flags will be read
# from either of the env vars FFLAGS or FCFLAGS.  If both are set, but are the
# same, then FFLAGS it used (which is the same as FCFLAGS).  However, if both
# are set but are not equal, then a FATAL_ERROR is raised and CMake configure
# processing is stopped.
#
# To be used in a parent project, this module must be included after
#
#   project(${PROJECT_NAME}  NONE)
#
# is called but before the compilers are defined and processed using:
#
#   enable_language(<LANG>)
#
# For example, one would do:
#
#   project(${PROJECT_NAME}  NONE)
#   ...
#   set(USE_XSDK_DEFAULTS_DEFAULT TRUE) # Set to false if desired
#   include("${CMAKE_CURRENT_SOURCE_DIR}/stdk/XSDKDefaults.cmake")
#   ...
#   enable_language(C)
#   enable_language(C++)
#   enable_language(Fortran)
#
# The variable `USE_XSDK_DEFAULTS_DEFAULT` is used as the default for the
# cache var `USE_XSDK_DEFAULTS`.  That way, a project can decide if it wants
# XSDK defaults turned on or off by default and users can independently decide
# if they want the CMake project to use standard XSDK behavior or raw CMake
# behavior.
#
# By default, the XSDKDefaults.cmake module assumes that the project will need
# C, C++, and Fortran.  If any language is not needed then, set
# XSDK_ENABLE_C=OFF, XSDK_ENABLE_CXX=OFF, or XSDK_ENABLE_Fortran=OFF *before*
# including this module.  Note, these variables are *not* cache vars because a
# project either does or does not have C, C++ or Fortran source files, the
# user has nothing to do with this so there is no need for cache vars.  The
# parent CMake project just needs to tell XSDKDefault.cmake what languages is
# needs or does not need.
#
# For example, if the parent CMake project only needs C, then it would do:
#
#   project(${PROJECT_NAME}  NONE)'
#   ...
#   set(USE_XSDK_DEFAULTS_DEFAULT TRUE)
#   set(XSDK_ENABLE_CXX OFF)
#   set(XSDK_ENABLE_Fortran OFF)
#   include("${CMAKE_CURRENT_SOURCE_DIR}/stdk/XSDKDefaults.cmake")
#   ...
#   enable_langauge(C)
#
# This module code will announce when it sets any variables.
#

#
# Helper functions
#

if (NOT COMMAND PRINT_VAR)
  function(print_var  VAR_NAME)
    message("${VAR_NAME} = '${${VAR_NAME}}'")
  endfunction()
endif()

if (NOT COMMAND SET_DEFAULT)
  macro(set_default VAR)
    if ("${${VAR}}" STREQUAL "")
      set(${VAR} ${ARGN})
    endif()
  endmacro()
endif()

#
# XSDKDefaults.cmake control variables
#

# USE_XSDK_DEFAULTS
if ("${USE_XSDK_DEFAULTS_DEFAULT}" STREQUAL "")
  set(USE_XSDK_DEFAULTS_DEFAULT  FALSE)
endif()
set(USE_XSDK_DEFAULTS  ${USE_XSDK_DEFAULTS_DEFAULT}  CACHE  BOOL
  "Use XSDK defaults and behavior.")
print_var(USE_XSDK_DEFAULTS)

set_default(XSDK_ENABLE_C  TRUE)
set_default(XSDK_ENABLE_CXX  TRUE)
set_default(XSDK_ENABLE_Fortran  TRUE)

# Handle the compiler and flags for a language
macro(xsdk_handle_lang_defaults  CMAKE_LANG_NAME  ENV_LANG_NAME
  ENV_LANG_FLAGS_NAMES
  )

  # Announce using env var ${ENV_LANG_NAME}
  if (NOT "$ENV{${ENV_LANG_NAME}}" STREQUAL "" AND
    "${CMAKE_${CMAKE_LANG_NAME}_COMPILER}" STREQUAL ""
    )
    message("-- " "XSDK: Setting CMAKE_${CMAKE_LANG_NAME}_COMPILER from env var"
      " ${ENV_LANG_NAME}='$ENV{${ENV_LANG_NAME}}'!")
    set(CMAKE_${CMAKE_LANG_NAME}_COMPILER "$ENV{${ENV_LANG_NAME}}" CACHE FILEPATH
      "XSDK: Set by default from env var ${ENV_LANG_NAME}")
  endif()

  # Announce using env var ${ENV_LANG_FLAGS_NAME}
  foreach(ENV_LANG_FLAGS_NAME  ${ENV_LANG_FLAGS_NAMES})
    if (NOT "$ENV{${ENV_LANG_FLAGS_NAME}}" STREQUAL "" AND
      "${CMAKE_${CMAKE_LANG_NAME}_FLAGS}" STREQUAL ""
      )
      message("-- " "XSDK: Setting CMAKE_${CMAKE_LANG_NAME}_FLAGS from env var"
        " ${ENV_LANG_FLAGS_NAME}='$ENV{${ENV_LANG_FLAGS_NAME}}'!")
      set(CMAKE_${CMAKE_LANG_NAME}_FLAGS "$ENV{${ENV_LANG_FLAGS_NAME}} " CACHE  STRING
        "XSDK: Set by default from env var ${ENV_LANG_FLAGS_NAME}")
      # NOTE: CMake adds the space after $ENV{${ENV_LANG_FLAGS_NAME}} so we
      # duplicate that here!
    endif()
  endforeach()

endmacro()


#
# Set XSDK Defaults
#

# Set default compilers and flags
if (USE_XSDK_DEFAULTS)

  # Handle env vars for languages C, C++, and Fortran

  if (XSDK_ENABLE_C)
    xsdk_handle_lang_defaults(C  CC  CFLAGS)
  endif()

  if (XSDK_ENABLE_CXX)
    xsdk_handle_lang_defaults(CXX  CXX  CXXFLAGS)
  endif()

  if (XSDK_ENABLE_Fortran)
    set(ENV_FFLAGS "$ENV{FFLAGS}")
    set(ENV_FCFLAGS "$ENV{FCFLAGS}")
    if (
      (NOT "${ENV_FFLAGS}" STREQUAL "") AND (NOT "${ENV_FCFLAGS}" STREQUAL "")
      AND
      ("${CMAKE_Fortran_FLAGS}" STREQUAL "")
      )
      if (NOT "${ENV_FFLAGS}" STREQUAL "${ENV_FCFLAGS}")
        message(FATAL_ERROR "Error, env vars FFLAGS='${ENV_FFLAGS}' and"
          " FCFLAGS='${ENV_FCFLAGS}' are both set in the env but are not equal!")
      endif()
    endif()
    xsdk_handle_lang_defaults(Fortran  FC  "FFLAGS;FCFLAGS")
  endif()
  
  # Set XSDK defaults for other CMake variables
  
  if ("${BUILD_SHARED_LIBS}"  STREQUAL  "")
    message("-- " "XSDK: Setting default BUILD_SHARED_LIBS=TRUE")
    set(BUILD_SHARED_LIBS  TRUE  CACHE  BOOL  "Set by default in XSDK mode")
  endif()
  
  if ("${CMAKE_BUILD_TYPE}"  STREQUAL  "")
    message("-- " "XSDK: Setting default CMAKE_BUILD_TYPE=DEBUG")
    set(CMAKE_BUILD_TYPE  DEBUG  CACHE  STRING  "Set by default in XSDK mode")
  endif()

endif()
