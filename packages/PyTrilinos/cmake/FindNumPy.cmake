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
