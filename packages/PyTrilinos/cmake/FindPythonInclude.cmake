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

# - Find the python include path

# This module determines where the python include file Python.h is,
# based on the current value of Python3_EXECUTABLE. This code sets the
# following variable:
#
#  PYTHON_INCLUDE_PATH  = path to where Python.h is found
#

IF(Python3_EXECUTABLE)

  # Obtain the Python version string
  EXECUTE_PROCESS(COMMAND
    ${Python3_EXECUTABLE} -c 
    "import sys; print('python' + sys.version[:3])"
    OUTPUT_VARIABLE PYVERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # Obtain the Python prefix path
  EXECUTE_PROCESS(COMMAND
    ${Python3_EXECUTABLE} -c 
    "import sys; print(sys.prefix)"
    OUTPUT_VARIABLE PYPREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # Obtain the Python include path
  EXECUTE_PROCESS(COMMAND
    ${Python3_EXECUTABLE} -c
    "import sysconfig; print(sysconfig.get_paths()['include'])"
    OUTPUT_VARIABLE PYTHON_INCLUDE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # # Check ${PYPREFIX}/include/${PYVERSION}
  # SET(CANDIDATE1 ${PYPREFIX}/include/${PYVERSION})

  # # Verify theat Python.h exists within the candidate path
  # IF(EXISTS "${CANDIDATE1}")
  #   IF(EXISTS ${CANDIDATE1}/Python.h)
  #     SET(PYTHON_INCLUDE_FOUND TRUE)
  #     SET(PYTHON_INCLUDE_PATH ${CANDIDATE1})
  #   ENDIF(EXISTS ${CANDIDATE1}/Python.h)
  # ELSE(EXISTS "${CANDIDATE1}")

  #   # Check ${PYPREFIX}/include/${PYVERSION}
  #   SET(CANDIDATE2 ${PYPREFIX}/Headers)

  #   # Verify theat Python.h exists within the candidate path
  #   IF(EXISTS "${CANDIDATE2}")
  #     IF(EXISTS ${CANDIDATE2}/Python.h)
  #       SET(PYTHON_INCLUDE_FOUND TRUE)
  #       SET(PYTHON_INCLUDE_PATH ${CANDIDATE2})
  #     ENDIF(EXISTS ${CANDIDATE2}/Python.h)
  #   ENDIF(EXISTS "${CANDIDATE2}")
  # ENDIF(EXISTS "${CANDIDATE1}")

  MESSAGE(STATUS "PYTHON_INCLUDE_PATH is ${PYTHON_INCLUDE_PATH}")

ENDIF(Python3_EXECUTABLE)
