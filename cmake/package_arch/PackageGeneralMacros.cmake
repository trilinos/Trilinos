# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

INCLUDE(AppendSet)
INCLUDE(AssertDefined)


#
# Function that extracts all of the required and optional
# items for a given class of package lists
#

FUNCTION(PACKAGE_GATHER_ENABLED_ITEMS PACKAGE_NAME LISTTYPE_PREFIX LISTTYPE_POSTFIX
  GATHERED_ITEMS_LIST
  )

  SET(GATHERED_ITEMS_LIST_TMP
    ${${PACKAGE_NAME}_${LISTTYPE_PREFIX}_REQUIRED_DEP_${LISTTYPE_POSTFIX}}
    )

  FOREACH(ITEM
    ${${PACKAGE_NAME}_${LISTTYPE_PREFIX}_OPTIONAL_DEP_${LISTTYPE_POSTFIX}}
    )
    ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_${ITEM})
    IF (${PACKAGE_NAME}_ENABLE_${ITEM})
      APPEND_SET(GATHERED_ITEMS_LIST_TMP ${ITEM})
    ENDIF()
  ENDFOREACH()

  SET(${GATHERED_ITEMS_LIST} ${GATHERED_ITEMS_LIST_TMP} PARENT_SCOPE)

ENDFUNCTION()


#
# Function that sorts a list of TPLs
# into reverse order for link order
# purposes
#

FUNCTION(PACKAGE_SORT_LIST MASTER_LIST LIST_VAR)

  #MESSAGE("PACKAGE_SORT_LIST:")
  #PRINT_VAR(MASTER_LIST)
  #PRINT_VAR(LIST_VAR)
  #PRINT_VAR(${LIST_VAR})

  SET(SORTED_LIST)

  FOREACH(TPL ${MASTER_LIST})
    LIST(FIND ${LIST_VAR} ${TPL} TPL_IDX)
    #PRINT_VAR(TPL)
    #PRINT_VAR(TPL_IDX)
    IF (NOT TPL_IDX EQUAL -1)
      APPEND_SET(SORTED_LIST ${TPL})
    ENDIF()
  ENDFOREACH()

  #PRINT_VAR(SORTED_LIST)

  SET(${LIST_VAR} ${SORTED_LIST} PARENT_SCOPE)

ENDFUNCTION()


#
# Function that appends the TPL paths and libraries
#

FUNCTION(PACKAGE_APPEND_PATHS_LIBS PREFIX LIST EXTRA_DEP_LIBS_ARG)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE STREQUAL "MAX")
    MESSAGE("\nPACKAGE_APPEND_PATHS_LIBS: ${PREFIX} ${LIST} ${EXTRA_DEP_LIBS_ARG}")
  ENDIF()
  SET(EXTRA_DEP_LIBS_ARG_TMP ${${EXTRA_DEP_LIBS_ARG}})
  FOREACH(ITEM ${LIST})
    ASSERT_DEFINED(${PREFIX}${ITEM}_LIBRARIES)
    ASSERT_DEFINED(${PREFIX}${ITEM}_INCLUDE_DIRS)
    ASSERT_DEFINED(${PREFIX}${ITEM}_LIBRARY_DIRS)
    APPEND_SET(EXTRA_DEP_LIBS_ARG_TMP ${${PREFIX}${ITEM}_LIBRARIES})
    INCLUDE_DIRECTORIES(${${PREFIX}${ITEM}_INCLUDE_DIRS})
    IF(PREFIX)
      # TODO: Is there a better way to know if we need this?
      # We want LINK_DIRECTORIES for TPLs but not packages.
      LINK_DIRECTORIES(${${PREFIX}${ITEM}_LIBRARY_DIRS})
    ENDIF()
    SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS ${${PREFIX}${ITEM}_LIBRARY_DIRS})
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE STREQUAL "MAX")
      PRINT_VAR(${PREFIX}${ITEM}_LIBRARIES)
      PRINT_VAR(${PREFIX}${ITEM}_INCLUDE_DIRS)
      PRINT_VAR(${PREFIX}${ITEM}_LIBRARY_DIRS)
    ENDIF()
  ENDFOREACH()
  SET(${EXTRA_DEP_LIBS_ARG} ${EXTRA_DEP_LIBS_ARG_TMP} PARENT_SCOPE)
ENDFUNCTION()


#
# Function that sorts and appends all the items in a dependency list
# for TPLs or packages
#

FUNCTION(PACKAGE_SORT_AND_APPEND_PATHS_LIBS MASTER_SORT_LIST LIST PREFIX
  EXTRA_DEP_LIBS_ARG
  )

  #MESSAGE("PACKAGE_SORT_AND_APPEND_PATHS_LIBS:")
  #PRINT_VAR(MASTER_SORT_LIST)
  #PRINT_VAR(LIST)
  #PRINT_VAR(PREFIX)
  #PRINT_VAR(EXTRA_DEP_LIBS_ARG)

  SET(LOCAL_LIST ${LIST})
  #PRINT_VAR(LOCAL_LIST)

  IF (LOCAL_LIST)

    PACKAGE_SORT_LIST("${MASTER_SORT_LIST}" LOCAL_LIST)
    #PRINT_VAR(LOCAL_LIST)

    SET(EXTRA_DEP_LIBS_ARG_TMP ${${EXTRA_DEP_LIBS_ARG}})
    PACKAGE_APPEND_PATHS_LIBS("${PREFIX}" "${LOCAL_LIST}" EXTRA_DEP_LIBS_ARG_TMP)
    SET(${EXTRA_DEP_LIBS_ARG} ${EXTRA_DEP_LIBS_ARG_TMP} PARENT_SCOPE)

  ENDIF()

ENDFUNCTION()
