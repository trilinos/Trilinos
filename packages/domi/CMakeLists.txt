# -*- cmake -*-

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

INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# Define the package
TRIBITS_PACKAGE(Domi DISABLE_STRONG_WARNINGS)

#
# Standardized debugging options
TRIBITS_ADD_DEBUG_OPTION()

#
# Set up package-specific options
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_ABC
  HAVE_DOMI_ARRAY_BOUNDSCHECK
  "Enable Domi feature for array bounds checking."
  ${${PROJECT_NAME}_ENABLE_DEBUG}
  )

########################################################################

# Python interpreter, library, and header files are required.
#SET(CMAKE_MODULE_PATH
#  ${CMAKE_PROJECT_DIR}/packages/PyTrilinos/cmake/FindPythonInclude.cmake
#  ${CMAKE_MODULE_PATH})
# SET(CMAKE_PREFIX_PATH
#   ${CMAKE_PROJECT_DIR}/packages/PyTrilinos/cmake
#   ${CMAKE_PREFIX_PATH})
# FIND_PACKAGE(PythonInterp 2.2 REQUIRED)
# FIND_PACKAGE(PythonLibs 2.2 REQUIRED)
# FIND_PACKAGE(PythonInclude REQUIRED)

# Python module numpy is required
# FIND_PACKAGE(NumPy 1.0 REQUIRED)

# Simple wrapper interface generator (SWIG) is required
# SET(SWIG_VERSION_FLOOR   2.0.0)
# SET(SWIG_VERSION_CEILING 2.0.8)
# FIND_PACKAGE(SWIG REQUIRED)
# IF(${SWIG_VERSION} VERSION_LESS ${SWIG_VERSION_FLOOR})
#   MESSAGE(FATAL_ERROR
#     "SWIG version " ${SWIG_VERSION}
#     " is less than minimum version " ${SWIG_VERSION_FLOOR}
#     )
# ENDIF()
# IF(${SWIG_VERSION} VERSION_GREATER ${SWIG_VERSION_CEILING})
#   MESSAGE(FATAL_ERROR
#     "SWIG version " ${SWIG_VERSION}
#     " is greater than maximum version " ${SWIG_VERSION_CEILING}
#     )
# ENDIF()
# ${SWIG_USE_FILE} is hard-coded to be the installed UseSWIG.cmake
# file.  We want to include the local, modified version.
# INCLUDE(UseSWIG)

########################################################################

#
# Add the libraries, tests, and examples
ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Add the Parameter List documentation
TRIBITS_ADD_TEST_DIRECTORIES(doc/parameterList)

#
# Do standard postprocessing
TRIBITS_PACKAGE_POSTPROCESS()
