# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Find Python executable which is needed for dependency file building
macro(tribits_find_python)
  set(PythonInterp_FIND_VERSION_MIN "2.6")
  if ("${PythonInterp_FIND_VERSION_DEFAULT}" STREQUAL "")
    set(PythonInterp_FIND_VERSION_DEFAULT "${PythonInterp_FIND_VERSION_MIN}")
  endif()
  advanced_set(PythonInterp_FIND_VERSION  ${PythonInterp_FIND_VERSION_DEFAULT}
    CACHE  STRING
    "Default version of Python to find (must be ${PythonInterp_FIND_VERSION_DEFAULT} or greater")
  if (PythonInterp_FIND_VERSION  VERSION_LESS  "${PythonInterp_FIND_VERSION_MIN}")
    message_wrapper(FATAL_ERROR  "Error,"
      " PythonInterp_FIND_VERSION=${PythonInterp_FIND_VERSION} < ${PythonInterp_FIND_VERSION_MIN}"
      " is not allowed!" )
  endif()
  advanced_set(PythonInterp_MUST_BE_FOUND FALSE CACHE BOOL "Require Python to be found or not.")
  if (${PROJECT_NAME}_REQUIRES_PYTHON)
    set(PythonInterp_REQUIRED_ARG "REQUIRED")
  else()
    set(PythonInterp_REQUIRED_ARG "")
  endif()
  set(FIND_PythonInterp_ARGS PythonInterp ${PythonInterp_REQUIRED_ARG}) 
  if (TRIBITS_FIND_PYTHON_UNITTEST)
    set(PYTHON_EXECUTABLE  ${PYTHON_EXECUTABLE_UNITTEST_VAL})
  else()
    find_package(${FIND_PythonInterp_ARGS})
  endif()
endmacro()


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
    print_var(PYTHON_EXECUTABLE)
    if (${PROJECT_NAME}_REQUIRES_PYTHON  AND  PYTHON_EXECUTABLE  STREQUAL "")
      message_wrapper(FATAL_ERROR "Error, PYTHON_EXECUTABLE='' but"
        " ${PROJECT_NAME}_REQUIRES_PYTHON=${${PROJECT_NAME}_REQUIRES_PYTHON}!" )
    endif()
  else()
    message_wrapper("-- " "NOTE: Skipping check for Python because"
      " ${PROJECT_NAME}_USES_PYTHON='${${PROJECT_NAME}_USES_PYTHON}'")
  endif()
endmacro()
