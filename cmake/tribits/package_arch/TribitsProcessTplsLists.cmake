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


INCLUDE(TribitsConstants)
INCLUDE(TribitsListHelpers)

INCLUDE(PrintVar)
INCLUDE(Split)

#
# @MACRO: TRIBITS_REPOSITORY_DEFINE_TPLS()
#
# Define the list of `TriBITS TPLs`_ for a given `TriBITS Repository`_ which
# includes the TPL name, find module, and classification .  This macro is
# typically called from inside of the repository's `<repoDir>/TPLsList.cmake`_
# file.
#
# Usage::
#
#   TRIBITS_REPOSITORY_DEFINE_TPLS(
#     <tpl0_name>   <tpl0_findmod>  <tpl0_classif>
#     <tpl1_name>   <tpl1_findmod>  <tpl1_classif>
#     ...
#     )
#
# This macro sets up a 2D array of ``NumTPLS`` by ``NumColumns`` listing out
# the `TriBITS TPLs`_ for a `TriBITS Repository`_.  Each row (with 3 entries)
# specifies a TPL which contains the columns (ordered 0-2):
#
# 0. **TPL** (``<tpli_name>``): The name of the TriBITS TPL ``<tplName>``.
#    This name must be unique across all other TriBITS TPLs in this or any
#    other TriBITS repo that might be combined into a single TriBITS project
#    meta-build (see `Globally unique TriBITS TPL names`_).  However, a TPL
#    can be redefined from an upstream repo (see below).  The name should be a
#    valid identifier (e.g. matches the regex ``[a-zA-Z_][a-zA-Z0-9_]*``).
#    TPL names typically use mixed case (e.g. ``SomeTpl`` and not
#    ``SOMETPL``).
#
# 1. **FINDMOD** (``<tpli_findmod>``): The relative path for the find module,
#    usually with the name ``FindTPL<tplName>.cmake``.  This path is relative
#    to the repository base directory.  If just the base path for the find
#    module is given, ending with ``"/"`` (e.g. ``"cmake/tpls/"``), then the
#    find module will be assumed to be under that this directory with the
#    standard name (e.g. ``cmake/tpls/FindTPL<tplName>.cmake``).  A standard
#    way to write a ``FindTPL<tplName>.cmake`` module is to use the function
#    `TRIBITS_TPL_DECLARE_LIBRARIES()`_.
#
# 2. **CLASSIFICATION** (``<pkgi_classif>``): Gives the `SE Package Test
#    Group`_ `PT`_, `ST`_, or `EX`_ and the maturity level ``EP``, ``RS``,
#    ``PG``, ``PM``, ``GRS``, ``GPG``, ``GPM``, ``UM``.  These are separated
#    by a coma with no space in between such as ``"RS,PT"`` for a "Research
#    Stable", "Primary Tested" package.  No spaces are allowed so that CMake
#    treats this a one field in the array.  The maturity level can be left off
#    in which case it is assumed to be ``UM`` for "Unspecified Maturity".
#
# A TPL defined in a upstream repo can listed again in a downstream repo,
# which allows redefining the find module that is used to specify the TPL.
# This allows downstream repos to add additional requirements for a given TPL
# (i.e. add more libraries, headers, etc.).  However, the downstream repo's
# find module file must find the TPL components that are fully compatible with
# the upstream's find module.
#
# This macro just sets the variable::
#
#   ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS
#
# in the current scope.  The advantages of using this macro instead of
# directly setting this variable are that the macro:
#
# * Asserts that the variable ``REPOSITORY_NAME`` is defined and set
#
# * Avoids having to hard-code the assumed repository name
#   ``${REPOSITORY_NAME}``.  This provides more flexibility for how other
#   TriBITS projects choose to name a given TriBITS repo (i.e. the name of
#   repo subdirs).
#
# * Avoids misspelling the name of the variable
#   ``${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS``.  If one misspells
#   the name of a macro, it is an immediate error in CMake.
#
MACRO(TRIBITS_REPOSITORY_DEFINE_TPLS)
  ASSERT_DEFINED(REPOSITORY_NAME)
  SET(${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS "${ARGN}")
ENDMACRO()


#
# Macro that processes the list of TPLs
#
# This macro reads from the variable
# ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS for a given repository and
# and fills the variables ${PROJECT_NAME}_TPLS, ${PROJECT_NAME}_NUM_TPLS,
# ${PROJECT_NAME}_REVERSE_TPLS.  For each TPL, it also sets the variable
# ${TPL_NAME}_FINDMOD and ${TPL_NAME}_TESTGROUP.
#

MACRO(TRIBITS_PROCESS_TPLS_LISTS  REPOSITORY_NAME  REPOSITORY_DIR)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("TRIBITS_PROCESS_TPLS_LISTS:  '${REPOSITORY_NAME}'  '${REPOSITORY_DIR}'")
  ENDIF()

  #SET(TRIBITS_PROCESS_TPLS_LISTS_DEBUG ON)
  SET(TRIBITS_PROCESS_TPLS_LISTS_DEBUG OFF)

  SET(TPL_NAME_OFFSET 0)
  SET(TPL_FINDMOD_OFFSET 1)
  SET(TPL_CLASSIFICATION_OFFSET 2)
  SET(TPL_NUM_COLUMNS 3)

  LIST(LENGTH ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS
    ${REPOSITORY_NAME}_CURR_NUM_TPLS_FULL)
  MATH(EXPR ${REPOSITORY_NAME}_CURR_NUM_TPLS
    "${${REPOSITORY_NAME}_CURR_NUM_TPLS_FULL}/${TPL_NUM_COLUMNS}")

  IF (${REPOSITORY_NAME}_CURR_NUM_TPLS GREATER 0)

    MATH(EXPR ${REPOSITORY_NAME}_LAST_TPL_IDX
      "${${REPOSITORY_NAME}_CURR_NUM_TPLS}-1")

    FOREACH(TPL_IDX RANGE ${${REPOSITORY_NAME}_LAST_TPL_IDX})

      IF (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        PRINT_VAR(TPL_IDX)
      ENDIF()

      # Get fields for this TPL

      MATH(EXPR TPL_NAME_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_NAME_OFFSET}")
      LIST(GET ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_NAME_IDX}
        TPL_NAME)
      IF (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        PRINT_VAR(TPL_NAME)
      ENDIF()

      MATH(EXPR TPL_FINDMOD_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_FINDMOD_OFFSET}")
      LIST(GET ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_FINDMOD_IDX}
        TPL_FINDMOD)
      IF (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        PRINT_VAR(TPL_FINDMOD)
      ENDIF()

      MATH(EXPR TPL_CLASSIFICATION_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_CLASSIFICATION_OFFSET}")
      LIST(GET ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_CLASSIFICATION_IDX}
        TPL_CLASSIFICATION)
      IF (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        PRINT_VAR(TPL_CLASSIFICATION)
      ENDIF()

      # ToDo: Parse out TESTGROUP and MATURITYLEVEL (Trilinos #6042)
      SET(TPL_TESTGROUP ${TPL_CLASSIFICATION})

      TRIBITS_UPDATE_PS_PT_SS_ST(TPL  ${TPL_NAME}  TPL_TESTGROUP)

      # Update TPLS list (unless the TPL already exists)

      IF (${TPL_NAME}_FINDMOD)
        # If the variable ${TPL_NAME}_FINDMOD already exists, then this TPL
        # has already been defined in a previous repository.  In this case, we
        # will just leave the TPL in its current position.
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE("-- " "NOTE: The TPL ${TPL_NAME} has already been defined so leaving it"
            " in the same location and not adding it again!")
        ENDIF()
      ELSE()
        LIST(APPEND ${PROJECT_NAME}_TPLS ${TPL_NAME})
      ENDIF()

      # Set ${TPL_NAME}_TESTGROUP

      IF (TPL_TESTGROUP STREQUAL PT
        OR TPL_TESTGROUP STREQUAL ST
        OR TPL_TESTGROUP STREQUAL TT
        OR TPL_TESTGROUP STREQUAL EX
        )
      ELSE()
        MESSAGE(FATAL_ERROR "Error the TPL classification '${TPL_TESTGROUP}'"
          " for the TPL ${TPL_NAME} is not a valid classification." )
      ENDIF()

      IF (NOT ${TPL_NAME}_TESTGROUP) # Allow for testing override
        SET(${TPL_NAME}_TESTGROUP ${TPL_TESTGROUP})
      ENDIF()

      # Set ${TPL_NAME}_FINDMOD

      #PRINT_VAR(REPOSITORY_DIR)

      IF ("${REPOSITORY_DIR}" STREQUAL "." OR IS_ABSOLUTE ${TPL_FINDMOD})
        SET(REPOSITORY_DIR_AND_SEP "")
      ELSE()
        SET(REPOSITORY_DIR_AND_SEP "${REPOSITORY_DIR}/")
      ENDIF()
      #PRINT_VAR(REPOSITORY_DIR_AND_SEP)

      SET(TPL_FINDMOD "${REPOSITORY_DIR_AND_SEP}${TPL_FINDMOD}")
      #PRINT_VAR(TPL_FINDMOD)

      SET(TPL_FINDMOD_STD_NAME "FindTPL${TPL_NAME}.cmake")

      IF (TPL_FINDMOD)
        STRING(REGEX MATCH ".+/$" FINDMOD_IS_DIR "${TPL_FINDMOD}")
        #PRINT_VAR(FINDMOD_IS_DIR)
        IF (FINDMOD_IS_DIR)
          SET(${TPL_NAME}_FINDMOD "${TPL_FINDMOD}${TPL_FINDMOD_STD_NAME}")
        ELSE()
          SET(${TPL_NAME}_FINDMOD ${TPL_FINDMOD})
        ENDIF()
      ELSE()
        SET(${TPL_NAME}_FINDMOD ${TPL_FINDMOD_STD_NAME})
      ENDIF()

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(${TPL_NAME}_FINDMOD)
      ENDIF()

      # Set the enable cache variable for ${TPL_NAME}

      MULTILINE_SET(DOCSTR
        "Enable support for the TPL ${TPL_NAME} in all supported ${PROJECT_NAME} packages."
        "  This can be set to 'ON', 'OFF', or left empty ''."
        )
      SET_CACHE_ON_OFF_EMPTY( TPL_ENABLE_${TPL_NAME} "" ${DOCSTR} )

      # 2008/11/25: rabartl: Above, we use the prefix TPL_ instead of
      # ${PROJECT_NAME}_ in order to make it clear that external TPLs are
      # different from packages so users don't get confused and
      # think that the project actually includes some TPL when it does not!

    ENDFOREACH()

  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PROJECT_NAME}_TPLS)
  ENDIF()

  # Get the final length

  LIST(LENGTH ${PROJECT_NAME}_TPLS ${PROJECT_NAME}_NUM_TPLS)
  PRINT_VAR(${PROJECT_NAME}_NUM_TPLS)

  # Create a reverse list for later use

  IF (${PROJECT_NAME}_TPLS)
    SET(${PROJECT_NAME}_REVERSE_TPLS ${${PROJECT_NAME}_TPLS})
    LIST(REVERSE ${PROJECT_NAME}_REVERSE_TPLS)
  ELSE()
    SET(${PROJECT_NAME}_REVERSE_TPLS)
  ENDIF()

ENDMACRO()
