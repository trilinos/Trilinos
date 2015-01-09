# @HEADER
# ***********************************************************************
#
#          PyTrilinos: Python Interfaces to Trilinos Packages
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
