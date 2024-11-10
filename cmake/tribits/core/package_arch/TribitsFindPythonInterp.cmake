# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# TriBITS Wrapper for finding Python (or not) for a TriBITS project.
macro(tribits_find_python_interp)
  if (${PROJECT_NAME}_REQUIRES_PYTHON)
    set(${PROJECT_NAME}_USES_PYTHON  TRUE)
  endif()
  if ("${${PROJECT_NAME}_USES_PYTHON}" STREQUAL "")
    # Unless the project tells us they can use Python or not, let's go ahead
    # and look for Python in case some packages want to use it.
    set(${PROJECT_NAME}_USES_PYTHON  TRUE)
  endif()
  if (${PROJECT_NAME}_USES_PYTHON)
    tribits_find_python()
    print_var(Python3_EXECUTABLE)
    if (${PROJECT_NAME}_REQUIRES_PYTHON  AND  Python3_EXECUTABLE  STREQUAL "")
      message_wrapper(FATAL_ERROR "Error, Python3_EXECUTABLE='' but"
        " ${PROJECT_NAME}_REQUIRES_PYTHON=${${PROJECT_NAME}_REQUIRES_PYTHON}!" )
    endif()
  else()
    message_wrapper("-- " "NOTE: Skipping check for Python because"
      " ${PROJECT_NAME}_USES_PYTHON='${${PROJECT_NAME}_USES_PYTHON}'")
  endif()
endmacro()


# Find Python executable which is needed for dependency file building
macro(tribits_find_python)
  tribits_find_python_set_python3_find_version()
  tribits_find_python_backward_compatible_python_executable()
  tribits_find_python_find_python3()
endmacro()


macro(tribits_find_python_set_python3_find_version)
  # Get minimum version of Python to find
  set(${PROJECT_NAME}_Python3_FIND_VERSION_MIN "3.6")
  if ("${${PROJECT_NAME}_Python3_FIND_VERSION_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_Python3_FIND_VERSION_DEFAULT
      "${${PROJECT_NAME}_Python3_FIND_VERSION_MIN}")
  endif()
  advanced_set(${PROJECT_NAME}_Python3_FIND_VERSION
    ${${PROJECT_NAME}_Python3_FIND_VERSION_DEFAULT}
    CACHE  STRING
    "Default version of Python to find (must be ${${PROJECT_NAME}_Python3_FIND_VERSION_DEFAULT} or greater")
  if (${PROJECT_NAME}_Python3_FIND_VERSION  VERSION_LESS
      "${${PROJECT_NAME}_Python3_FIND_VERSION_MIN}"
    )
    message_wrapper(FATAL_ERROR  "Error,"
      " ${PROJECT_NAME}_Python3_FIND_VERSION=${${PROJECT_NAME}_Python3_FIND_VERSION} < ${${PROJECT_NAME}_Python3_FIND_VERSION_MIN}"
      " is not allowed!" )
  endif()
endmacro()


macro(tribits_find_python_backward_compatible_python_executable)
  # Provide backward compatibility for user setting PYTHON_EXECUTABLE
  if ((NOT "${PYTHON_EXECUTABLE}" STREQUAL "") AND ("${Python3_EXECUTABLE}" STREQUAL ""))
    tribits_deprecated("Python3_EXECUTABLE being set by default to PYTHON_EXECUTABLE = '${PYTHON_EXECUTABLE}' is deprecated!")
    set(Python3_EXECUTABLE "${PYTHON_EXECUTABLE}" CACHE FILEPATH
      "Set by default to PYTHON_EXECUTABLE!")
  endif()
endmacro()


macro(tribits_find_python_find_python3)
  # Find Python3
  if (${PROJECT_NAME}_REQUIRES_PYTHON)
    set(Python3_REQUIRED_ARG "REQUIRED")
  else()
    set(Python3_REQUIRED_ARG "")
  endif()
  set(FIND_Python3_ARGS
    Python3 ${${PROJECT_NAME}_Python3_FIND_VERSION} ${Python3_REQUIRED_ARG})
  if(DEFINED Python3_EXECUTABLE)
    # Already defined (even if it is set to empty), so no need to call anything!
  elseif (TRIBITS_FIND_PYTHON_UNITTEST)
    set(Python3_EXECUTABLE  ${Python3_EXECUTABLE_UNITTEST_VAL})
  else()
    find_package(${FIND_Python3_ARGS})
  endif()
endmacro()



