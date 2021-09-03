# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
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
