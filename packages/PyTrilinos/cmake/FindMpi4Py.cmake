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

# - Find the python mpi4py module

# This cmake module determines whether the python module mpi4py exists
# for the current PYTHON_EXECUTABLE. This code sets the following
# variables:
#
#  Mpi4Py_FOUND       = Set to TRUE if mpi4py is found
#  Mpi4Py_VERSION     = Mpi4Py version number
#  Mpi4Py_INCLUDE_DIR = Path to mpi4py include files
#

#
# If Mpi4Py is required and python executable does not exist, then send
# an error
IF(NOT PYTHON_EXECUTABLE)
  IF(Mpi4Py_FIND_REQUIRED)
    MESSAGE(SEND_ERROR
      "Python executable not found, so required Mpi4Py module not found"
      )
  ENDIF(Mpi4Py_FIND_REQUIRED)
#
# Continue processing if python executable is known
ELSE(NOT PYTHON_EXECUTABLE)

  # Retrieve the Mpi4Py version
  EXECUTE_PROCESS(COMMAND
    ${PYTHON_EXECUTABLE} -c "import mpi4py; print mpi4py.__version__"
    OUTPUT_VARIABLE Mpi4Py_VERSION
    ERROR_VARIABLE  Mpi4Py_VERSION_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # If Mpi4Py_VERSION_ERROR does not exist, then we know Mpi4Py exists;
  # now look for the Mpi4Py include directory
  IF(NOT Mpi4Py_VERSION_ERROR)
    EXECUTE_PROCESS(COMMAND
      ${PYTHON_EXECUTABLE} -c "import mpi4py; print mpi4py.get_include()"
      OUTPUT_VARIABLE Mpi4Py_INCLUDE_DIR
      ERROR_VARIABLE  Mpi4Py_INCLUDE_ERROR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )

    # Handle the QUIETLY and REQUIRED arguments and set Mpi4Py_FOUND to
    # TRUE if all listed variables are TRUE
    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(Mpi4Py
      FOUND_VAR     Mpi4Py_FOUND
      REQUIRED_VARS Mpi4Py_VERSION Mpi4Py_INCLUDE_DIR
      VERSION_VAR   Mpi4Py_VERSION)

    # Version checking: If a version check is requested, compare
    # Mpi4Py_VERSION to the requested version
    IF(Mpi4Py_FIND_VERSION)
      IF(${Mpi4Py_VERSION} VERSION_LESS ${Mpi4Py_FIND_VERSION})
	MESSAGE(FATAL_ERROR
	  "Mpi4Py version " ${Mpi4Py_VERSION}
	  " is less than required version " ${Mpi4Py_FIND_VERSION}
	  )
      ENDIF(${Mpi4Py_VERSION} VERSION_LESS ${Mpi4Py_FIND_VERSION})
    ENDIF(Mpi4Py_FIND_VERSION)

  #
  # A Mpi4Py version error means that Mpi4Py was not found
  ELSE(NOT Mpi4Py_VERSION_ERROR)
    IF(Mpi4Py_FIND_REQUIRED)
      MESSAGE(SEND_ERROR
	"Required Mpi4Py python module not found"
	)
    ENDIF(Mpi4Py_FIND_REQUIRED)

  ENDIF(NOT Mpi4Py_VERSION_ERROR)

ENDIF(NOT PYTHON_EXECUTABLE)
