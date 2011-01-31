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

# - Find the python include path

# This module determines where the python include file Python.h is,
# based on the current value of PYTHON_EXECUTABLE. This code sets the
# following variable:
#
#  PYTHON_INCLUDE_PATH  = path to where Python.h is found
#

IF(PYTHON_EXECUTABLE)
  # Obtain the candidate path for python include
  EXECUTE_PROCESS(COMMAND
    ${PYTHON_EXECUTABLE} -c 
    "import sys; print sys.prefix + '/include/python' + sys.version[:3]"
    OUTPUT_VARIABLE CANDIDATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  # Verify theat Python.h exists within the candidate path
  IF(EXISTS "${CANDIDATE}")
    IF(EXISTS ${CANDIDATE}/Python.h)
      SET(PYTHON_INCLUDE_FOUND TRUE)
      SET(PYTHON_INCLUDE_PATH ${CANDIDATE})
    ENDIF(EXISTS ${CANDIDATE}/Python.h)
  ENDIF(EXISTS "${CANDIDATE}")
  # Obtain the candidate path for python library
  EXECUTE_PROCESS(COMMAND
    ${PYTHON_EXECUTABLE} -c 
    "import sys; print 'python' + sys.version[:3]"
    OUTPUT_VARIABLE PYVERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  EXECUTE_PROCESS(COMMAND
    ${PYTHON_EXECUTABLE} -c 
    "import sys; print sys.prefix + '/lib/python' + sys.version[:3] + '/config'"
    OUTPUT_VARIABLE CANDIDATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  # Verify theat Python.h exists within the candidate path
  IF(EXISTS "${CANDIDATE}")
    FOREACH(SUFFIX dll dylib so a)
      SET(LIB_CANDIDATE ${CANDIDATE}/lib${PYVERSION}.${SUFFIX})
      IF(EXISTS ${LIB_CANDIDATE})
	SET(PYTHON_LIB_FOUND TRUE)
	SET(PYTHON_LIBRARIES ${LIB_CANDIDATE})
	BREAK()
      ENDIF(EXISTS ${LIB_CANDIDATE})
    ENDFOREACH(SUFFIX)
  ENDIF(EXISTS "${CANDIDATE}")
ENDIF(PYTHON_EXECUTABLE)

# Set any python variables that did not get set by the above logic
IF(PYTHON_INCLUDE_FOUND AND PYTHON_LIB_FOUND)
  # We're good
ELSE(PYTHON_INCLUDE_FOUND AND PYTHON_LIB_FOUND)
  FIND_PACKAGE(PythonLibs)
ENDIF(PYTHON_INCLUDE_FOUND AND PYTHON_LIB_FOUND)
