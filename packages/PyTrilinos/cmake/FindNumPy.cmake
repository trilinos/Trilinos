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
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
# @HEADER

# - Find the python numpy module

# This cmake module determines whether the python module numpy exists
# for the current PYTHON_EXECUTABLE. This code sets the following
# variable:
#
#  NumPy_FOUND       = Set to TRUE if numpy is found
#  NumPy_VERSION     = NumPy version number
#  NumPy_INCLUDE_DIR = Path to numpy include files
#

#
# If NumPy is required and python executable does not exist, then send
# an error
IF(NOT PYTHON_EXECUTABLE)
  IF(NumPy_FIND_REQUIRED)
    MESSAGE(SEND_ERROR
      "Python executable not found, so required NumPy module not found"
      )
  ENDIF(NumPy_FIND_REQUIRED)
#
# Continue processing if python executable is known
ELSE(NOT PYTHON_EXECUTABLE)

  # Retrieve the NumPy version
  EXECUTE_PROCESS(COMMAND
    ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.__version__"
    OUTPUT_VARIABLE NumPy_VERSION
    ERROR_VARIABLE  NumPy_VERSION_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # If NumPy_VERSION_ERROR does not exist, then we know NumPy exists;
  # now look for the NumPy include directory
  IF(NOT NumPy_VERSION_ERROR)
    EXECUTE_PROCESS(COMMAND
      ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
      OUTPUT_VARIABLE NumPy_INCLUDE_DIR
      ERROR_VARIABLE  NumPy_INCLUDE_ERROR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )

    # If there is a NumPy include error, it is because NumPy is older
    # and the wrong function name was called
    IF(NumPy_INCLUDE_ERROR)
      EXECUTE_PROCESS(COMMAND
	${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_numpy_include()"
        OUTPUT_VARIABLE NumPy_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
	)
    ENDIF(NumPy_INCLUDE_ERROR)

    # Handle the QUIETLY and REQUIRED arguments and set NumPy_FOUND to
    # TRUE if all listed variables are TRUE
    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(NumPy
      FOUND_VAR     NumPy_FOUND
      REQUIRED_VARS NumPy_VERSION NumPy_INCLUDE_DIR
      VERSION_VAR   NumPy_VERSION)

    #
    # Version checking: If a version check is requested, compare
    # NumPy_VERSION to the requested version
    IF(NumPy_FIND_VERSION)
      IF(${NumPy_VERSION} VERSION_LESS ${NumPy_FIND_VERSION})
	MESSAGE(FATAL_ERROR
	  "NumPy version " ${NumPy_VERSION}
	  " is less than required version " ${NumPy_FIND_VERSION}
	  )
      ENDIF(${NumPy_VERSION} VERSION_LESS ${NumPy_FIND_VERSION})
    ENDIF(NumPy_FIND_VERSION)

  #
  # A NumPy version error means that NumPy was not found
  ELSE(NOT NumPy_VERSION_ERROR)
    IF(NumPy_FIND_REQUIRED)
      MESSAGE(SEND_ERROR
	"Required NumPy python module not found"
	)
    ENDIF(NumPy_FIND_REQUIRED)

  ENDIF(NOT NumPy_VERSION_ERROR)

ENDIF(NOT PYTHON_EXECUTABLE)
