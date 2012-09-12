# @HEADER
#
# ***********************************************************************
#
#           Amesos2: Templated Direct Sparse Solver Package
#                  Copyright 2011 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ***********************************************************************
#
# @HEADER


INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

# If we use Trilinos_version then "make install" will fail.
#APPEND_SET(HEADERS
#${PROJECT_SOURCE_DIR}/Trilinos_version.h
#  )

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so that MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Amesos2.hpp
  Amesos2_AbstractConcreteMatrixAdapter_def.hpp
  Amesos2_Control.hpp
  Amesos2_Factory.hpp
  Amesos2_FunctionMap.hpp
  Amesos2_TpetraRowMatrix_AbstractMatrixAdapter.hpp
  Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_decl.hpp
  Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_def.hpp
  Amesos2_TpetraCrsMatrix_MatrixAdapter.hpp
  Amesos2_TpetraCrsMatrix_MatrixAdapter_decl.hpp
  Amesos2_TpetraCrsMatrix_MatrixAdapter_def.hpp
  Amesos2_MatrixAdapter.hpp
  Amesos2_MatrixAdapter_decl.hpp
  Amesos2_MatrixAdapter_def.hpp
  Amesos2_AbstractConcreteMatrixAdapter.hpp
  Amesos2_ctassert.hpp
  Amesos2_ConcreteMatrixAdapter.hpp
  Amesos2_ConcreteMatrixAdapter_decl.hpp
  Amesos2_ConcreteMatrixAdapter_def.hpp
  Amesos2_MatrixHelper.hpp
  Amesos2_MatrixTraits.hpp
  Amesos2_Meta.hpp
  Amesos2_MultiVecAdapter.hpp
  Amesos2_MultiVecAdapter_decl.hpp
  Amesos2_MultiVecAdapter_def.hpp
  Amesos2_Solver.hpp
  Amesos2_Solver_decl.hpp
  Amesos2_SolverCore.hpp
  Amesos2_SolverCore_decl.hpp
  Amesos2_SolverCore_def.hpp
  Amesos2_SolverTraits.hpp
  Amesos2_Status.hpp
  Amesos2_Timers.hpp
  Amesos2_TpetraMultiVecAdapter.hpp
  Amesos2_TpetraMultiVecAdapter_decl.hpp
  Amesos2_TpetraMultiVecAdapter_def.hpp
  Amesos2_TypeMap.hpp
  Amesos2_TypeDecl.hpp
  Amesos2_Util.hpp
  Amesos2_Version.hpp
  )

# Set sources

APPEND_SET(SOURCES
  Amesos2_Control.cpp
  Amesos2_Factory.cpp
  Amesos2_MatrixAdapter.cpp
  Amesos2_ConcreteMatrixAdapter.cpp
  Amesos2_MultiVecAdapter.cpp
  Amesos2_SolverCore.cpp
  Amesos2_Solver.cpp
  Amesos2_TpetraRowMatrix_AbstractMatrixAdapter.cpp
  Amesos2_TpetraCrsMatrix_MatrixAdapter.cpp
  Amesos2_TpetraMultiVecAdapter.cpp
  Amesos2_Util.cpp
  Amesos2_Version.cpp
  )

# Set additional headers and sources based on package-specific
# configuration options.  In Amesos2, each solver interface should be
# added only if that solver has been enabled by the user.

IF (${PACKAGE_NAME}_ENABLE_SuperLU)
  APPEND_SET(HEADERS
    Amesos2_Superlu.hpp
    Amesos2_Superlu_decl.hpp
    Amesos2_Superlu_def.hpp
    Amesos2_Superlu_FunctionMap.hpp
    Amesos2_Superlu_TypeMap.hpp
    )

  APPEND_SET(SOURCES
    Amesos2_Superlu.cpp
    Amesos2_Superlu_TypeMap.cpp
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_SuperLUMT)
  APPEND_SET(HEADERS
    Amesos2_Superlumt.hpp
    Amesos2_Superlumt_FunctionMap.hpp
    Amesos2_Superlumt_TypeMap.hpp
    Amesos2_Superlumt_decl.hpp
    Amesos2_Superlumt_def.hpp
    )

  APPEND_SET(SOURCES
    Amesos2_Superlumt.cpp
    Amesos2_Superlumt_TypeMap.cpp
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_SuperLUDist) # This naming syntax is consistent with Trilinos
  APPEND_SET(HEADERS
    Amesos2_Superludist.hpp
    Amesos2_Superludist_FunctionMap.hpp
    Amesos2_Superludist_TypeMap.hpp
    Amesos2_Superludist_decl.hpp
    Amesos2_Superludist_def.hpp
    )

  APPEND_SET(SOURCES
    Amesos2_Superludist.cpp
#    Amesos2_Superludist_TypeMap.cpp
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_PARDISO_MKL)
  APPEND_SET(HEADER
    Amesos2_PardisoMKL.hpp
    Amesos2_PardisoMKL_FunctionMap.hpp
    Amesos2_PardisoMKL_TypeMap.hpp
    Amesos2_PardisoMKL_decl.hpp
    Amesos2_PardisoMKL_def.hpp
    )

  APPEND_SET(SOURCES
    Amesos2_PardisoMKL.cpp
  )
ENDIF()


IF (${PACKAGE_NAME}_ENABLE_Epetra)
  APPEND_SET(HEADERS
    Amesos2_EpetraRowMatrix_AbstractMatrixAdapter.hpp
    Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_decl.hpp
    Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_def.hpp
    Amesos2_EpetraCrsMatrix_MatrixAdapter.hpp
    Amesos2_EpetraCrsMatrix_MatrixAdapter_decl.hpp
    Amesos2_EpetraCrsMatrix_MatrixAdapter_def.hpp
    Amesos2_EpetraMultiVecAdapter.hpp
    Amesos2_EpetraMultiVecAdapter_decl.hpp
    Amesos2_EpetraMultiVecAdapter_def.hpp
    )

  APPEND_SET(SOURCES
    # Amesos2_EpetraCrsMatrixAdapter.cpp
    # Amesos2_EpetraJadMatrixAdapter.cpp
    # Amesos2_EpetraMsrMatrixAdapter.cpp
    # Amesos2_EpetraRowMatrixAdapter.cpp
    # Amesos2_EpetraVbrMatrixAdapter.cpp
    Amesos2_EpetraRowMatrix_AbstractMatrixAdapter.cpp
    Amesos2_EpetraCrsMatrix_MatrixAdapter.cpp
    Amesos2_EpetraMultiVecAdapter.cpp
    )
ENDIF()

ADD_SUBDIRECTORY(SuiteSparse/KLU2/Source)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  amesos2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
