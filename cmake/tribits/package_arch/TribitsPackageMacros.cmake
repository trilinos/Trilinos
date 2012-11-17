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

INCLUDE(TribitsPackageSetupCompilerFlags)
INCLUDE(TribitsPackageWritePackageConfig)
INCLUDE(TribitsGeneralMacros)

INCLUDE(ParseVariableArguments)
INCLUDE(GlobalNullSet)
INCLUDE(AppendGlobalSet)
INCLUDE(PrintVar)
INCLUDE(PrependSet)
INCLUDE(PrependGlobalSet)
INCLUDE(RemoveGlobalDuplicates)
INCLUDE(TribitsAddOptionAndDefine)


#
# Macro that defines the package architecture system varaibles used to link
# different SE packages together
#
# See README.DEPENDENCIES for information on what these varaibles mean and how
# they are used.
#

MACRO(TRIBITS_DEFINE_LINKAGE_VARS PACKAGE_NAME_IN)
  GLOBAL_NULL_SET(${PACKAGE_NAME_IN}_INCLUDE_DIRS)
  GLOBAL_NULL_SET(${PACKAGE_NAME_IN}_LIBRARY_DIRS)
  GLOBAL_NULL_SET(${PACKAGE_NAME_IN}_LIBRARIES)
ENDMACRO()


#
# Macro that defines varaibles that create global targets
#

MACRO(TRIBITS_DEFINE_TARGET_VARS PARENT_PACKAGE_NAME_IN)
  GLOBAL_NULL_SET(${PARENT_PACKAGE_NAME_IN}_LIB_TARGETS)
  GLOBAL_NULL_SET(${PARENT_PACKAGE_NAME_IN}_ALL_TARGETS)
ENDMACRO()

#
# Set up some common varaibles used in the creation of an SE package
#

MACRO(TRIBITS_SET_COMMON_VARS PACKAGE_NAME_IN)

  STRING(TOUPPER ${PACKAGE_NAME_IN} PACKAGE_NAME_UC)

  # Write TRIBITS_PACKAGE versions of common variables
  SET(PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  SET(PACKAGE_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")

  # Get the name of the directory this ${PROJECT_NAME} package is in
  FILE(TO_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR} STANDARD_PACKAGE_SOURCE_DIR)
  STRING(REGEX REPLACE "/.+/(.+)" "\\1" PACKAGE_DIR_NAME "${STANDARD_PACKAGE_SOURCE_DIR}")

ENDMACRO()


#
# TRIBITS_PACKAGE_DECL(...): Macro called at the very beginning of a ${PROJECT_NAME}
# package's top-level CMakeLists.txt file when there are subpackages.
#
# Usage is:
#
#   TRIBITS_PACKAGE_DECL(
#     <packageName>
#     [ENABLE_SHADOWING_WARNINGS]
#     [DISABLE_STRONG_WARNINGS]
#     [CLEANED]
#     [DISABLE_CIRCULAR_REF_DETECTION_FAILURE]
#     )
#
# The arguments are:
#
#   <packageName>
#
#     Gives the name of the Package, mostly just for checking and
#     documentation purposes.
#
#   ENABLE_SHADOWING_WARNINGS
#
#     If specified, then shadowing warnings will be turned on for supported
#     platforms/compilers.  The default is for shadowing warnings to be turned
#     off.  Note that this can be overridden globally by setting the cache
#     variable ${PROJECT_NAME}_ENABLE_SHADOWING_WARNINGS.
#
#   DISABLE_STRONG_WARNINGS
#
#     If specified, then all strong warnings will be turned off, if they are
#     not already turned off by global cache variables.
#
#   CLEANED
#
#     If specified, then warnings will be promoted to errors for all defined
#     warnings.
#
#   DISABLE_CIRCULAR_REF_DETECTION_FAILURE
#
#     If specified, then the standard grep looking for RCPNode circular
#     references that causes tests to fail will be disabled.  Note that if
#     these warnings are being produced then it means that the test is leaking
#     memory and user like may also be leaking memory.
#
MACRO(TRIBITS_PACKAGE_DECL PACKAGE_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_DECL: ${PACKAGE_NAME_IN}")
  ENDIF()
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
    #prefix
    PARSE
    #lists
    ""
    #options
    "CLEANED;ENABLE_SHADOWING_WARNINGS;DISABLE_STRONG_WARNINGS;DISABLE_CIRCULAR_REF_DETECTION_FAILURE"
    ${ARGN}
    )

  #
  # B) Assert that the global and local package names are the same!
  #

  IF (DEFINED PACKAGE_NAME)
    IF (NOT ${PACKAGE_NAME_IN} STREQUAL ${PACKAGE_NAME})
      MESSAGE(FATAL_ERROR "Error, the package-defined package name"
        " '${PACKAGE_NAME_IN}' is not the same as the package name"
        " defined at the global level '${PACKAGE_NAME}'")
    ENDIF()
  ENDIF()

  #
  # C) Set up the CMake support for this ${PROJECT_NAME} package and define some
  # top-level varaibles.
  #

  TRIBITS_SET_COMMON_VARS(${PACKAGE_NAME_IN})

  # Set up the compile flags for the package
  TRIBITS_SETUP_COMPILER_FLASGS()

  # Set up circular reference detection test failure
  IF (PARSE_DISABLE_CIRCULAR_REF_DETECTION_FAILURE)
    SET(${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OFF)
  ELSE()
    SET(${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE ON)
  ENDIF()

  # Set up parent package linkage varaibles
  TRIBITS_DEFINE_TARGET_VARS(${PACKAGE_NAME})

ENDMACRO()


#
# TRIBITS_PACKAGE_DEF(): Macro called after subpackages are processed in order to
# handle the libraries, tests, and examples of the final package.
#

MACRO(TRIBITS_PACKAGE_DEF)

  # Reset since it was changed by the subpackages
  SET(PACKAGE_NAME ${PARENT_PACKAGE_NAME})

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_DEF: ${PACKAGE_NAME}")
  ENDIF()

  IF (NOT ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("\n${PACKAGE_NAME} not enabled so exiting package processing")
    ENDIF()
    RETURN()
  ENDIF()
  
  # Reset in case were changed by subpackages
  TRIBITS_SET_COMMON_VARS(${PACKAGE_NAME})

  # Define package linkage varaibles
  TRIBITS_DEFINE_LINKAGE_VARS(${PACKAGE_NAME})

ENDMACRO()


#
# TRIBITS_PACKAGE(...): Macro called at the very beginning of a ${PROJECT_NAME}
# package's top-level CMakeLists.txt file.
#
# Usage is:
#
#   TRIBITS_PACKAGE(
#     <packageName>
#     [ENABLE_SHADOWING_WARNINGS]
#     [DISABLE_STRONG_WARNINGS]
#     [CLEANED]
#     [DISABLE_CIRCULAR_REF_DETECTION_FAILURE]
#     )
#
# See TRIBITS_PACKAGE_DECL(...) for a description of the arguments.
#
MACRO(TRIBITS_PACKAGE PACKAGE_NAME_IN)
  TRIBITS_PACKAGE_DECL(${PACKAGE_NAME_IN} ${ARGN})
  TRIBITS_PACKAGE_DEF()
ENDMACRO()


#
# Macro called to add a set of test directories for a package
#
# This macro only needs to be called from the top most CMakeList.txt file for
# which all subdirectories area all "tests".
#
# This macro can be called several times within a package and it will have the
# right effect.
#
# This macro defines hooks for inserting certain types of behavior in a
# uniform way.
#
MACRO(TRIBITS_ADD_TEST_DIRECTORIES)

  IF(${PACKAGE_NAME}_ENABLE_TESTS OR ${PARENT_PACKAGE_NAME}_ENABLE_TESTS)
    FOREACH(TEST_DIR ${ARGN})
      ADD_SUBDIRECTORY(${TEST_DIR})
    ENDFOREACH()
  ENDIF()

ENDMACRO()


#
# Common options to add to a package
#


MACRO(TRIBITS_ADD_DEBUG_OPTION)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_DEBUG
    HAVE_${PACKAGE_NAME_UC}_DEBUG
    "Enable a host of runtime debug checking."
    ${${PROJECT_NAME}_ENABLE_DEBUG}
    )
ENDMACRO()


MACRO(TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION)
  OPTION(
    ${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
     "Enable Teuchos time monitors for package ${PACKAGE_NAME}"
    ${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR}
    )
ENDMACRO()


MACRO(TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION)
  ADVANCED_SET(
    ${PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS  ${${PROJECT_NAME}_SHOW_DEPRECATED_WARNINGS}
    CACHE BOOL
    "Show warnings about deprecated code in ${PACKAGE_NAME}"
    )
  ADVANCED_SET(
    ${PACKAGE_NAME}_HIDE_DEPRECATED_CODE  ${${PROJECT_NAME}_HIDE_DEPRECATED_CODE}
    CACHE BOOL
    "Fully exclude deprecated code in ${PACKAGE_NAME}"
    )
ENDMACRO()


MACRO(TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
    HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION
    "Enable the use of explicit template instantiation."
    ${${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION}
    )
ENDMACRO()

MACRO(TRIBITS_ADD_ETI_SUPPORT)
  APPEND_GLOBAL_SET(${PROJECT_NAME}_ETI_PACKAGES ${PACKAGE_NAME})
  GLOBAL_NULL_SET(${PACKAGE_NAME}_ETI_LIBRARYSET)
ENDMACRO()


#
# Macro called to add a set of example directories for a package
#
# This macro only needs to be called from the top most CMakeList.txt file for
# which all subdirectories area all "examples".
#
# This macro can be called several times within a package and it will have the
# right effect.
#
# This macro defines hooks for inserting certain types of behavior in a
# uniform way.
#

MACRO(TRIBITS_ADD_EXAMPLE_DIRECTORIES)

  IF(${PACKAGE_NAME}_ENABLE_EXAMPLES OR ${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES)
    FOREACH(EXAMPLE_DIR ${ARGN})
      ADD_SUBDIRECTORY(${EXAMPLE_DIR})
    ENDFOREACH()
  ENDIF()

ENDMACRO()


#
# Function that sets up package linkage linkage variables in case the package
# has no libraries.
#

FUNCTION(PACAKGE_SETUP_DEPENDENCY_VARS_IF_NO_LIBS)

  IF (${PACKAGE_NAME}_SUBPACKAGES)

    # A package with subpackages should get all of its dependency vars from
    # its enabled subpackages.

    GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS)
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS)
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES)

    SET(SUBPACKAGE_IDX 0)
    FOREACH(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})
  
      SET(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
      SET(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
  
      IF (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
        PREPEND_GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS
          ${${SUBPACKAGE_FULLNAME}_INCLUDE_DIRS})
        PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS
          ${${SUBPACKAGE_FULLNAME}_LIBRARY_DIRS})
        PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES
          ${${SUBPACKAGE_FULLNAME}_LIBRARIES})
      ENDIF()
  
      MATH(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")
  
    ENDFOREACH()
  

  ELSEIF(NOT ${PACKAGE_NAME}_INCLUDE_DIRS)

    # No libraries have been defined for this package so we are going to set
    # them based on this package's dependencies.

    TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  LIB  LINK_LIBS) 

    TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)

    GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT  INCLUDE_DIRECTORIES)
    GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT  PACKAGE_LIBRARY_DIRS)

    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS  ${INCLUDE_DIRS_CURRENT})
    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS  ${LIBRARY_DIRS_CURRENT})
    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES  ${LINK_LIBS})

  ENDIF()

ENDFUNCTION()


#
# Helper macro for [SUB]TRIBITS_PACKAGE_POSTPROCESS()
#

MACRO(TRIBITS_PACKAGE_POSTPROCESS_COMMON)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_POSTPROCESS: ${PACKAGE_NAME}")
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

  IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    # Create the configure file so external projects can find packages with a
    # call to find_package(<package_name>)
    # This also creates the Makefile.export.* files.
    TRIBITS_WRITE_PACKAGE_CONFIG_FILE(${PACKAGE_NAME})
  ENDIF()
  
  SET(${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE TRUE
    CACHE INTERNAL "")

ENDMACRO()


#
# Macro called at the very end of a package's top-level CMakeLists.txt file
#

MACRO(TRIBITS_PACKAGE_POSTPROCESS)

  # Only parent packages have the targets (${PACKAGE_NAME}_libs and
  # (${PACKAGE_NAME}_all
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_POSTPROCESS: ${PACKAGE_NAME}")
    PRINT_VAR(${PACKAGE_NAME}_LIB_TARGETS)
    PRINT_VAR(${PACKAGE_NAME}_ALL_TARGETS)
  ENDIF()
  ADD_CUSTOM_TARGET(${PACKAGE_NAME}_libs DEPENDS ${${PACKAGE_NAME}_LIB_TARGETS})
  ADD_CUSTOM_TARGET(${PACKAGE_NAME}_all DEPENDS ${${PACKAGE_NAME}_ALL_TARGETS})

  PACAKGE_SETUP_DEPENDENCY_VARS_IF_NO_LIBS()
  TRIBITS_PACKAGE_POSTPROCESS_COMMON()

ENDMACRO()


#
# Macro that processes subpackages for packages that have them
#

MACRO(TRIBITS_PROCESS_SUBPACKAGES)

  #MESSAGE("TRIBITS_PROCESS_SUBPACKAGES: ${PARENT_PACKAGE_NAME}")

  SET(SUBPACKAGE_IDX 0)
  FOREACH(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})

    #MESSAGE("")
    #PRINT_VAR(SUBPACKAGE_IDX)
    #PRINT_VAR(TRIBITS_SUBPACKAGE)

    SET(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
    SET(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    #PRINT_VAR(SUBPACKAGE_FULLNAME)

    IF (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
      
      LIST(GET ${PARENT_PACKAGE_NAME}_SUBPACKAGE_DIRS ${SUBPACKAGE_IDX} SUBPACKAGE_DIR)
      #PRINT_VAR(SUBPACKAGE_DIR)

      DUAL_SCOPE_SET(${SUBPACKAGE_FULLNAME}_SOURCE_DIR
        ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/${SUBPACKAGE_DIR})
      DUAL_SCOPE_SET(${SUBPACKAGE_FULLNAME}_BINARY_DIR
        ${${PARENT_PACKAGE_NAME}_BINARY_DIR}/${SUBPACKAGE_DIR})

      ADD_SUBDIRECTORY(${${SUBPACKAGE_FULLNAME}_SOURCE_DIR}
        ${${SUBPACKAGE_FULLNAME}_BINARY_DIR})

    ENDIF()

    MATH(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")

  ENDFOREACH()

ENDMACRO()


#
# Append the local package's cmake directory in order to help pull in 
# configure-time testing macros
#

PREPEND_SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
