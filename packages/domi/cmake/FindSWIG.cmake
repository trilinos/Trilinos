# @HEADER
# ***********************************************************************
#
#     Domi: Multi-dimensional Distributed Linear Algebra Services
#                 Copyright (2014) Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia
# Corporation, the U.S. Government retains certain rights in this
# software.
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
# Questions? Contact William F. Spotz (wfspotz@sandia.gov)
#
# ***********************************************************************
# @HEADER

# This is a modification of the standard SWIG module that adds version
# checking.

# - Find SWIG
# This module finds an installed SWIG.  It sets the following variables:
#  SWIG_FOUND - set to true if SWIG is found
#  SWIG_DIR - the directory where swig is installed
#  SWIG_EXECUTABLE - the path to the swig executable
#  SWIG_VERSION   - the version number of the swig executable
#
# All informations are collected from the SWIG_EXECUTABLE so the
# version to be found can be changed from the command line by
# means of setting SWIG_EXECUTABLE
#

SET(SWIG_FOUND FALSE)

FIND_PROGRAM(SWIG_EXECUTABLE swig)

IF(SWIG_EXECUTABLE)
  EXECUTE_PROCESS(COMMAND ${SWIG_EXECUTABLE} -swiglib
    OUTPUT_VARIABLE SWIG_swiglib_output
    ERROR_VARIABLE SWIG_swiglib_error
    RESULT_VARIABLE SWIG_swiglib_result)

  IF(SWIG_swiglib_result) 
    IF(SWIG_FIND_REQUIRED)
      MESSAGE(SEND_ERROR "Command \"${SWIG_EXECUTABLE} -swiglib\" failed with output:\n${SWIG_swiglib_error}")
    ELSE(SWIG_FIND_REQUIRED)
      MESSAGE(STATUS "Command \"${SWIG_EXECUTABLE} -swiglib\" failed with output:\n${SWIG_swiglib_error}")
    ENDIF(SWIG_FIND_REQUIRED)
  ELSE(SWIG_swiglib_result)
    STRING(REGEX REPLACE "[\n\r]+" ";" SWIG_swiglib_output ${SWIG_swiglib_output})
    # force the path to be computed each time in case SWIG_EXECUTABLE has changed.
    SET(SWIG_DIR SWIG_DIR-NOTFOUND)
    FIND_PATH(SWIG_DIR swig.swg PATHS ${SWIG_swiglib_output})
    IF(SWIG_DIR)
      SET(SWIG_FOUND 1)
      SET(SWIG_USE_FILE ${CMAKE_ROOT}/Modules/UseSWIG.cmake)
      EXECUTE_PROCESS(COMMAND ${SWIG_EXECUTABLE} -version
        OUTPUT_VARIABLE SWIG_version_output
        ERROR_VARIABLE SWIG_version_output
        RESULT_VARIABLE SWIG_version_result)
      IF(SWIG_version_result)
        MESSAGE(SEND_ERROR "Command \"${SWIG_EXECUTABLE} -version\" failed with output:\n${SWIG_version_output}")
      ELSE(SWIG_version_result)
        STRING(REGEX REPLACE ".*SWIG Version[^0-9.]*\([0-9.]+\).*" "\\1"
          SWIG_version_output "${SWIG_version_output}")
        SET(SWIG_VERSION ${SWIG_version_output} CACHE STRING "Swig version" FORCE)
	# Begin local modification
	IF(SWIG_FIND_VERSION)
	  IF(${SWIG_VERSION} VERSION_LESS ${SWIG_FIND_VERSION})
	    MESSAGE(FATAL_ERROR
	      "SWIG version " ${SWIG_VERSION}
	      " is less than required version " ${SWIG_FIND_VERSION}
	      )
	  ENDIF(${SWIG_VERSION} VERSION_LESS ${SWIG_FIND_VERSION})
	ENDIF(SWIG_FIND_VERSION)
	# End local modification
      ENDIF(SWIG_version_result)
    ENDIF(SWIG_DIR)
  ENDIF(SWIG_swiglib_result)
ENDIF(SWIG_EXECUTABLE)

IF(NOT SWIG_FOUND)
  IF(NOT SWIG_FIND_QUIETLY)
    IF(SWIG_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "SWIG was not found. Please specify Swig executable location")
    ELSE(SWIG_FIND_REQUIRED)
      MESSAGE(STATUS "SWIG was not found. Please specify Swig executable location")
    ENDIF(SWIG_FIND_REQUIRED)
  ENDIF(NOT SWIG_FIND_QUIETLY)
ENDIF(NOT SWIG_FOUND)
