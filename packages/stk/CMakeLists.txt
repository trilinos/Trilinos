# Copyright (c) 2013, Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Governement retains certain rights in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#     * Neither the name of Sandia Corporation nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE_DECL( STK  )

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# B) Set up package-specific options
#

IF (Trilinos_ENABLE_KokkosCore)
  SET(HAVE_STK_KokkosCore ON)
ENDIF()

IF (Trilinos_ENABLE_STKIO)
  SET(HAVE_STKIO ON)
ENDIF()

IF (TPL_ENABLE_CUDA)
  SET(HAVE_STK_CUDA ON)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Trilinos
  HAVE_STK_Trilinos
  "Building as a Trilinos package"
  ON
  )

# set fortran mangling options

IF("${FC_FN_UNDERSCORE}" STREQUAL "NO_UNDER")
  SET(FORTRAN_NO_UNDERSCORE ON)
ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "UNDER")
  SET(FORTRAN_ONE_UNDERSCORE ON)
ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "SECOND_UNDER")
  SET(FORTRAN_TWO_UNDERSCORES ON)
ELSE()
  MESSAGE("Could not determine the Fortran mangling defaulting to no underscore.")
  SET(FORTRAN_NO_UNDERSCORE ON)
ENDIF()


#
# C) Add the libraries, tests, and examples
TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

#need to find out what to exclude
#TRIBITS_EXCLUDE_FILES(
#  )
  
#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()

