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

FUNCTION( TRIBITS_GATHER_ENABLED_ITEMS  PACKAGE_NAME  LISTTYPE_PREFIX 
  LISTTYPE_POSTFIX  GATHERED_ITEMS_LIST_OUT
  )

  #MESSAGE("TRIBITS_GATHER_ENABLED_ITEMS:  '${PACKAGE_NAME}'  '${LISTTYPE_PREFIX}'" 
  #  "  '${LISTTYPE_POSTFIX}'  '${GATHERED_ITEMS_LIST_OUT}'")

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

  #MESSAGE("TRIBITS_GATHER_ENABLED_ITEMS:"
  #  "  ${GATHERED_ITEMS_LIST_OUT} = ${GATHERED_ITEMS_LIST_TMP}")

  SET(${GATHERED_ITEMS_LIST_OUT} ${GATHERED_ITEMS_LIST_TMP} PARENT_SCOPE)

ENDFUNCTION()


#
# Function that does an in-place sort of a list of items according to the
# ordering in a master list
#
# NOTE: This function has wost-case N^2 complexity as the number of packages N
# or TPLs increases.  It actually has N * n complexity where N is the total
# number of packages/TPLs and n is the number of passed-in packages/TPLs.
# However, since N is not likely to ever be more than a few hundred, this is
# likely not going to be a big performance problem.  If this does become a
# performance problem, LIST(SORT ...) could be used but would require some
# work to build up the datastructures to make this very efficient.
#

FUNCTION(TRIBITS_SORT_LIST MASTER_LIST LIST_VAR)

  #MESSAGE("TRIBITS_SORT_LIST:")
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
# Function that appends the Package/TPL include and library paths for given
# list of enabled Packages/TPLs
#
# As a side effect of calling this function, INCLUDE_DIRECTORIES(...) to set
# all of the include directories for a given set of enabled Packages/TPLs.
#
# NOTE: The Packages/TPLs should be sorted in decending dependency order
# before calling this function.
#

# NOTE: Because this function may be called in cases where a package's
# required subpackages are not actually enabled (e.g. SEACAS subpackages)

#

FUNCTION( TRIBITS_APPEND_INCLUDE_AND_LINK_DIRS  PREFIX  LIST  EXTRA_DEP_LIBS_INOUT)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE STREQUAL "MAX")
    MESSAGE("\nPACKAGE_APPEND_INCLUDE_AND_LINK_DIRS: ${PREFIX} ${LIST} ${EXTRA_DEP_LIBS_INOUT}")
  ENDIF()
  SET(EXTRA_DEP_LIBS_INOUT_TMP ${${EXTRA_DEP_LIBS_INOUT}})
  FOREACH(ITEM ${LIST})
    IF (${PREFIX}${ITEM}_LIBRARIES)
      APPEND_SET(EXTRA_DEP_LIBS_ARG_TMP ${${PREFIX}${ITEM}_LIBRARIES})
    ENDIF()
    IF (${PREFIX}${ITEM}_INCLUDE_DIRS)
      INCLUDE_DIRECTORIES(${${PREFIX}${ITEM}_INCLUDE_DIRS})
    ENDIF()
    IF (${PREFIX}${ITEM}_LIBRARY_DIRS)
      IF (PREFIX)
        # TODO: Is there a better way to know if we need this?
        # We want LINK_DIRECTORIES for TPLs but not packages.
        LINK_DIRECTORIES(${${PREFIX}${ITEM}_LIBRARY_DIRS})
      ENDIF()
      SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS
        ${${PREFIX}${ITEM}_LIBRARY_DIRS})
    ENDIF()
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE STREQUAL "MAX")
      PRINT_VAR(${PREFIX}${ITEM}_LIBRARIES)
      PRINT_VAR(${PREFIX}${ITEM}_INCLUDE_DIRS)
      PRINT_VAR(${PREFIX}${ITEM}_LIBRARY_DIRS)
    ENDIF()
  ENDFOREACH()
  SET(${EXTRA_DEP_LIBS_INOUT} ${EXTRA_DEP_LIBS_ARG_TMP} PARENT_SCOPE)
ENDFUNCTION()


#
# Function that sorts and appends all the items in a dependency list for
# packages or TPLs.
#

FUNCTION( TRIBITS_SORT_AND_APPEND_INCLUDE_AND_LINK_DIRS_AND_LIBS
  MASTER_SORT_LIST  LIST  PREFIX
  EXTRA_DEP_LIBS_INOUT
  )

  #MESSAGE("TRIBITS_SORT_AND_APPEND_INCLUDE_AND_LINK_DIRS_AND_LIBS:")
  #PRINT_VAR(MASTER_SORT_LIST)
  #PRINT_VAR(LIST)
  #PRINT_VAR(PREFIX)
  #PRINT_VAR(EXTRA_DEP_LIBS_INOUT)

  SET(LOCAL_LIST ${LIST})
  #PRINT_VAR(LOCAL_LIST)

  IF (LOCAL_LIST)

    TRIBITS_SORT_LIST("${MASTER_SORT_LIST}"  LOCAL_LIST)
    #PRINT_VAR(LOCAL_LIST)

    SET(EXTRA_DEP_LIBS_ARG_TMP ${${EXTRA_DEP_LIBS_INOUT}})
    TRIBITS_APPEND_INCLUDE_AND_LINK_DIRS("${PREFIX}"
      "${LOCAL_LIST}" EXTRA_DEP_LIBS_ARG_TMP)
    SET(${EXTRA_DEP_LIBS_INOUT} ${EXTRA_DEP_LIBS_ARG_TMP} PARENT_SCOPE)

  ENDIF()

ENDFUNCTION()


#
# Fully process the include and link directories and list of libraries for a
# package's list of dependent packages for use in creating a library or an
# executable
#

FUNCTION( TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS
  PACKAGE_NAME  LIB_OR_TEST_ARG  EXTRA_DEP_LIBS_INOUT
  )

  TRIBITS_GATHER_ENABLED_ITEMS(${PACKAGE_NAME}  ${LIB_OR_TEST_ARG}
    PACKAGES  ALL_DEP_PACKAGES)

  SET(EXTRA_DEP_LIBS_TMP ${${EXTRA_DEP_LIBS_INOUT}})
  TRIBITS_SORT_AND_APPEND_INCLUDE_AND_LINK_DIRS_AND_LIBS(
    "${${PROJECT_NAME}_REVERSE_SE_PACKAGES}"
    "${ALL_DEP_PACKAGES}"  ""  EXTRA_DEP_LIBS_TMP)
  SET(${EXTRA_DEP_LIBS_INOUT} ${EXTRA_DEP_LIBS_TMP} PARENT_SCOPE)

ENDFUNCTION()


#
# Fully process the include and link directories and list of libraries for a
# package's list of dependent TPLs for use in creating a library or an
# executable
#

FUNCTION( TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS
  PACKAGE_NAME  LIB_OR_TEST_ARG  EXTRA_DEP_LIBS_INOUT
  )

  TRIBITS_GATHER_ENABLED_ITEMS(${PACKAGE_NAME}  ${LIB_OR_TEST_ARG}
    TPLS  ALL_TPLS)

  SET(EXTRA_DEP_LIBS_TMP ${${EXTRA_DEP_LIBS_INOUT}})
  TRIBITS_SORT_AND_APPEND_INCLUDE_AND_LINK_DIRS_AND_LIBS(
    "${${PROJECT_NAME}_REVERSE_TPLS}"
    "${ALL_TPLS}"  TPL_  EXTRA_DEP_LIBS_TMP)
  SET(${EXTRA_DEP_LIBS_INOUT} ${EXTRA_DEP_LIBS_TMP} PARENT_SCOPE)

ENDFUNCTION()
