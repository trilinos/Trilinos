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

#
# Include the TriBITS library macros

#
# Package-specific configuration options
IF(Domi_ENABLE_Epetra)
  SET(HAVE_EPETRA ON)
ENDIF()

IF(Domi_ENABLE_Tpetra)
  SET(HAVE_TPETRA ON)
ENDIF()

#
# Here we set the Domi_ORDINAL_TYPE to be int.  This should be
# sufficient for all conceivable situations, even for global ordinals
# for distributed MDVectors.  We do not run into the same upper bound
# on problem size that, say, Epetra does, because Domi is a
# multi-dimensional package, and each dimension can be the maximum
# value supported by int.
SET(${PACKAGE_NAME}_ORDINAL_TYPE int
  CACHE STRING
  "Select the data type for Domi_Ordinal (Domi::Ordinal)."
  )

#
# Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# Set the include directories
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# IF (Domi_ENABLE_PyTrilinos)
#   INCLUDE_DIRECTORIES(PYTHON_INCLUDE_PATH)
#   INCLUDE_DIRECTORIES(NumPy_INCLUDE_DIR)
# ENDIF (Domi_ENABLE_PyTrilinos)

# ADD_SUBDIRECTORY(PyTrilinos)

#
# Define the header and source files (and directories)
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

APPEND_SET(HEADERS
  Domi_ConfigDefs.hpp
  Domi_Version.hpp
  Domi_Utils.hpp
  Domi_Exceptions.hpp
  Domi_DefaultNode.hpp
  Domi_Slice.hpp
  Domi_MDIterator.hpp
  Domi_MDRevIterator.hpp
  Domi_MDArray.hpp
  Domi_MDArrayView.hpp
  Domi_MDArrayRCP.hpp
  Domi_MDComm.hpp
  Domi_MDMap.hpp
  Domi_MDVector.hpp
  Domi_getValidParameters.hpp
  )

APPEND_SET(SOURCES
  Domi_Version.cpp
  Domi_Utils.cpp
  Domi_Exceptions.cpp
  Domi_DefaultNode.cpp
  Domi_Slice.cpp
  Domi_MDComm.cpp
  Domi_getValidParameters.cpp
  )

#
# Define the targets for package's library(s)
TRIBITS_ADD_LIBRARY(
  domi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
