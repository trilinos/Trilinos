# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(CMakeBuildTypesList)

include(AdvancedSet)
include(AssertDefined)
include(MultilineSet)
include(DualScopeSet)
include(PrependCmndlineArgs)


#
# Macro that just defines the basic flags
#

macro(tribits_define_standard_compile_flags_vars  ENABLE_SHADOWING_WARNINGS)

  #
  # Setup and general flags
  #

  set(GENERAL_BUILD_FLAGS) # Applies to all builds, period

  if (${PROJECT_NAME}_ENABLE_CHECKED_STL)
    prepend_cmndline_args(GENERAL_BUILD_FLAGS "-D_GLIBCXX_DEBUG")
  endif()

  set(DEBUG_SYMBOLS_FLAGS "-g")

  set(GENERAL_DEBUG_FLAGS "${DEBUG_SYMBOLS_FLAGS} -O0")

  set(C_DEBUG_FLAGS "${GENERAL_DEBUG_FLAGS}")
  set(CXX_DEBUG_FLAGS "${GENERAL_DEBUG_FLAGS}")
  set(Fortran_DEBUG_FLAGS "${GENERAL_DEBUG_FLAGS}")

  set(GENERAL_RELEASE_FLAGS "-O3")

  set(C_RELEASE_FLAGS "${GENERAL_RELEASE_FLAGS} -DNDEBUG")
  set(CXX_RELEASE_FLAGS "${GENERAL_RELEASE_FLAGS} -DNDEBUG")
  set(Fortran_RELEASE_FLAGS "${GENERAL_RELEASE_FLAGS}")

  if (${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS)
    prepend_cmndline_args(GENERAL_BUILD_FLAGS "${DEBUG_SYMBOLS_FLAGS}")
  endif()

  if (${PROJECT_NAME}_COMMON_STRONG_COMPILE_WARNING_FLAGS)
    set(COMMON_STRONG_COMPILE_WARNING_FLAGS
       ${${PROJECT_NAME}_COMMON_STRONG_COMPILE_WARNING_FLAGS})
  else()
    multiline_set(COMMON_STRONG_COMPILE_WARNING_FLAGS
      " -pedantic" # Adds more static checking to remove non-ANSI GNU extensions
      " -Wall" # Enable a bunch of default warnings
      " -Wno-long-long" # Allow long long int since it is used by MPI, SWIG, etc.
      )
   endif()

  if (${PROJECT_NAME}_C_STRONG_COMPILE_WARNING_FLAGS)
    set(C_STRONG_COMPILE_WARNING_FLAGS ${${PROJECT_NAME}_C_STRONG_COMPILE_WARNING_FLAGS})
  else()
    if ("${${PROJECT_NAME}_C_Standard}" STREQUAL "")
      set(${PROJECT_NAME}_C_Standard_USED c99)
    else()
      set(${PROJECT_NAME}_C_Standard_USED ${${PROJECT_NAME}_C_Standard})
    endif()
    multiline_set(C_STRONG_COMPILE_WARNING_FLAGS
      ${COMMON_STRONG_COMPILE_WARNING_FLAGS}
      " -std=${${PROJECT_NAME}_C_Standard_USED}"
      )
  endif()

  if (${PROJECT_NAME}_CXX_STRONG_COMPILE_WARNING_FLAGS)
    set(CXX_STRONG_COMPILE_WARNING_FLAGS ${${PROJECT_NAME}_CXX_STRONG_COMPILE_WARNING_FLAGS})
  else()
    multiline_set(CXX_STRONG_COMPILE_WARNING_FLAGS
      ${COMMON_STRONG_COMPILE_WARNING_FLAGS}
      " -Wwrite-strings" # Checks for non-const char * copy of string constants
      )
  endif()

  if (${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS)
    set(LOCAL_ENABLE_SHADOWING_WARNINGS ON)
  elseif (${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS STREQUAL "OFF")
    set(LOCAL_ENABLE_SHADOWING_WARNINGS OFF)
  else()
    set(LOCAL_ENABLE_SHADOWING_WARNINGS ${ENABLE_SHADOWING_WARNINGS})
  endif()

  if (LOCAL_ENABLE_SHADOWING_WARNINGS)
    multiline_set(CXX_STRONG_COMPILE_WARNING_FLAGS
      ${CXX_STRONG_COMPILE_WARNING_FLAGS}
      " -Wshadow" # Warn about general shadowing issues
      " -Woverloaded-virtual" # Warn about hiding virtual functions
      )
  endif()

endmacro()
