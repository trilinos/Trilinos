# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
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
# ************************************************************************
# @HEADER

INCLUDE(TribitsPackageMacros)


#
# @MACRO: TRIBITS_SUBPACKAGE()
#
# Forward declare a `TriBITS Subpackage`_ called at the top of the
# subpackage's `<packageDir>/<spkgDir>/CMakeLists.txt`_ file.
#
# Usage::
#
#   TRIBITS_SUBPACKAGE(<spkgName>)
#
# Once called, the following local variables are in scope:
#
#   ``PARENT_PACKAGE_NAME``
#
#     The name of the parent package.
#
#   ``SUBPACKAGE_NAME``
#
#     The local name of the subpackage (does not contain
#     the parent package name).
#
#   ``SUBPACKAGE_FULLNAME``
#
#     The full project-level name of the subpackage (which includes the parent
#     package name at the beginning,
#     ``${PARENT_PACKAGE_NAME}${SUBPACKAGE_NAME}``).
#
#   ``PACKAGE_NAME``
#
#     Inside the subpackage, the same as ``SUBPACKAGE_FULLNAME``.
#
MACRO(TRIBITS_SUBPACKAGE SUBPACKAGE_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nSUBPACKAGE: ${SUBPACKAGE_NAME_IN}")
  ENDIF()

  # check that this is not being called from a package
  IF (NOT CURRENTLY_PROCESSING_SUBPACKAGE)
  # we are in a package

    MESSAGE(FATAL_ERROR "Cannot call TRIBITS_SUBPACKAGE() from a package."
    " Use TRIBITS_PACKAGE() instead"
    " ${CURRENT_PACKAGE_CMAKELIST_FILE}")

  ELSE()
  # We are in a subpackage

    # check to see if postprocess is called before subpackage
    IF(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      MESSAGE(FATAL_ERROR "TRIBITS_SUBPACKAGE_POSTPROCESS() called before TRIBITS_SUBPACKAGE()")
    ENDIF()

    # check to see if we have already called this macro
    IF(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED)
      MESSAGE(FATAL_ERROR "Already called TRIBITS_SUBPACKGE() for the"
	"${PARENT_PACKAGE_NAME} subpackage ${TRIBITS_SUBPACKAGE}")
    ENDIF()

    # make sure the name in the macro call matches the name in the packages cmake file
    IF (NOT ${SUBPACKAGE_NAME_IN} STREQUAL ${SUBPACKAGE_NAME})
      MESSAGE(FATAL_ERROR "Error, the package-defined subpackage name"
	" '${SUBPACKAGE_NAME_IN}' is not the same as the subpackage name"
	" '${SUBPACKAGE_NAME}' defined in the parent packages's"
	" Dependencies.cmake file")
    ENDIF()
  ENDIF()


  # To provide context for various macros
  SET(PACKAGE_NAME ${SUBPACKAGE_FULLNAME})

  SET(PARENT_PACKAGE_SOURCE_DIR "${PACKAGE_SOURCE_DIR}")
  SET(PARENT_PACKAGE_BINARY_DIR "${PACKAGE_BINARY_DIR}")

  # Now override the package-like variables
  TRIBITS_SET_COMMON_VARS(${SUBPACKAGE_FULLNAME})
  TRIBITS_DEFINE_LINKAGE_VARS(${SUBPACKAGE_FULLNAME})

  # Set flags that are used  to check that macros are called in the correct order
  SET(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED TRUE)
  SET(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED FALSE)

ENDMACRO()


#
# @MACRO: TRIBITS_SUBPACKAGE_POSTPROCESS()
#
# Macro that performs standard post-processing after defining a `TriBITS
# Subpackage`_ which is called at the bottom of a subpackage's
# `<packageDir>/<spkgDir>/CMakeLists.txt`_ file.
#
# Usage::
#
#   TRIBITS_SUBPACKAGE_POSTPROCESS()
#
# NOTE: It is unfortunate that a Subpackages's CMakeLists.txt file must call
# this macro but limitations of the CMake language make it necessary to do so.
#
MACRO(TRIBITS_SUBPACKAGE_POSTPROCESS)

  # check that this is not being called from a package
  IF (NOT CURRENTLY_PROCESSING_SUBPACKAGE)

  # This is being called from a package

    MESSAGE(FATAL_ERROR "Cannot call TRIBITS_SUBPACKAGE_POSTPROCESS() from a package."
    " Use TRIBITS_PACKAGE_POSTPROCESS() instead"
    " ${CURRENT_PACKAGE_CMAKELIST_FILE}")

  ELSE()
  # This is being caleld from a subpackage

    # check to make sure this has not already been called
    IF (${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      MESSAGE(FATAL_ERROR "Already called TRIBITS_SUBPACKGE_POSTPROCESS() for the"
        "${PARENT_PACKAGE_NAME} subpackage ${TRIBITS_SUBPACKAGE}")
    ENDIF()
  
    # make sure subpackage is called prior to subpackage postprocess
    IF(NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED)
      MESSAGE(FATAL_ERROR "TRIBITS_SUBPACKAGE() must be called before TRIBITS_SUBPACKAGE_POSTPROCESS()"
        "for the ${PARENT_PACKAGE_NAME} subpackage ${TRIBITS_SUBPACKAGE}")
    ENDIF()

  ENDIF()

  # Set flags that are used  to check that macros are called in the correct order
  DUAL_SCOPE_SET(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED TRUE)

  TRIBITS_PACKAGE_POSTPROCESS_COMMON()

ENDMACRO()
