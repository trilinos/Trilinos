# @HEADER
# ************************************************************************
#
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2010) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
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
