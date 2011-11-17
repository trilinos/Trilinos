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



#
# This file contains global-level macros that are specific to Trilinos
#


MACRO(TRIBITS_PROJECT_SETUP_EXTRA_OPTIONS)
  
  ADVANCED_SET(Trilinos_DATA_DIR  NOTFOUND
    CACHE PATH
    "Path TrilinosData directory to find more tests and other stuff" )
    
  IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting Trilinos_ENABLE_ForTrilinos=OFF"
      " because Trilinos_ENABLE_Fortran=OFF!"
      "\n***\n"
      )
    SET(Trilinos_ENABLE_ForTrilinos OFF)
  ENDIF()
    
  # ToDo: What is this and why is it needed?
  SET(TRILINOS_BUILD_SHARED_LIBS "@BUILD_SHARED_LIBS@")

ENDMACRO()


#
# Macro that defines Trilinos testing support
#

MACRO(TRIBITS_PROJECT_SETUP_TESTING_SUPPORT)

  CONFIGURE_FILE(
    ${${PROJECT_NAME}_TRIBITS_DIR}/ctest/CTestCustom.ctest.in
    ${${PROJECT_NAME}_BINARY_DIR}/CTestCustom.ctest
    )

ENDMACRO()


#
# Macro that defines Trilinos packaging options:
#
# ToDo: Seprate this into TriBITS general and Trilinos specific.
#

MACRO(TRIBITS_PROJECT_DEFINE_PACKAGING)

  SET(CPACK_SOURCE_IGNORE_FILES
    ${CPACK_SOURCE_IGNORE_FILES}
    /.git/
    ".gitignore"
    classicMakefile
    ".*.pyc"
    ${Trilinos_SOURCE_DIR}/cmake/tribits/common_tools/git
    ${Trilinos_SOURCE_DIR}/cmake/CMakeKitwareBacklog.txt
    ${Trilinos_SOURCE_DIR}/cmake/TODO
    ${Trilinos_SOURCE_DIR}/packages/ITAPS
    ${Trilinos_SOURCE_DIR}/packages/external
    ${Trilinos_SOURCE_DIR}/packages/jpetra
    ${Trilinos_SOURCE_DIR}/packages/cmmlib
    ${Trilinos_SOURCE_DIR}/packages/configure.ac
    ${Trilinos_SOURCE_DIR}/packages/configure
    ${Trilinos_SOURCE_DIR}/packages/Makefile.am
    ${Trilinos_SOURCE_DIR}/packages/Makefile.in
    ${Trilinos_SOURCE_DIR}/doc/[^b]
    ${Trilinos_SOURCE_DIR}/README_old
    ${Trilinos_SOURCE_DIR}/sampleScripts/old_autotools
    ${Trilinos_SOURCE_DIR}/sampleScripts/git-profiles
    ${Trilinos_SOURCE_DIR}/SIERRA
    ${Trilinos_SOURCE_DIR}/commonTools/test/coverage
    ${Trilinos_SOURCE_DIR}/commonTools/test/harness
    ${Trilinos_SOURCE_DIR}/commonTools/test/utilities/README
    ${Trilinos_SOURCE_DIR}/commonTools/test/utilities/dependencies
    ${Trilinos_SOURCE_DIR}/commonTools/test/utilities/packages
    ${Trilinos_SOURCE_DIR}/commonTools/test/utilities/r.*
    ${Trilinos_SOURCE_DIR}/commonTools/scripts
    ${Trilinos_SOURCE_DIR}/commonTools/release
    ${Trilinos_SOURCE_DIR}/packages/common/DoxyfilePackageTemplate
    ${Trilinos_SOURCE_DIR}/stamp-h.in
    ${Trilinos_SOURCE_DIR}/configure.ac
    ${Trilinos_SOURCE_DIR}/aclocal.m4
    ${Trilinos_SOURCE_DIR}/configure
    ${Trilinos_SOURCE_DIR}/Makefile.am
    ${Trilinos_SOURCE_DIR}/Makefile.in
    ${Trilinos_SOURCE_DIR}/bootstrap
    ${Trilinos_SOURCE_DIR}/config
  )
  
  #removing any packages not enabled from the tarball
  set(ENABLED_FLAG OFF)
  set(INCLUDE_EMPTY TRUE)
  TRIBITS_GET_ENABLED_LIST(${PROJECT_NAME}_PACKAGES ${PROJECT_NAME} ${ENABLED_FLAG} ${INCLUDE_EMPTY} 
    NON_ENABLED_PACKAGES NUM_NON_ENABLED)
  STRING(REPLACE " " ";" NON_ENABLED_PACKAGES "${NON_ENABLED_PACKAGES}")

  FOREACH(TRIBITS_PACKAGE ${NON_ENABLED_PACKAGES})
    #if the package is the TrilinosFramework we do not want to exclude it from the tarball
    #because that would exclude the cmake directory and the entire build system. So as a
    #special case we do not remove the TrilinosFramework from the tarball
    IF(NOT ${TRIBITS_PACKAGE} STREQUAL "TrilinosFramework")
      LIST(FIND ${PROJECT_NAME}_PACKAGES ${TRIBITS_PACKAGE} PACKAGE_IDX)
      LIST(GET ${PROJECT_NAME}_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
      
      #checking if we have a relative path to the package's files. Since the exclude is a
      #regular expression any "../" will be interpretted as <any char><any char>/ which
      #would never match the package's actual directory. There isn't a direct way in cmake
      #to convert a relative path into an absolute path with string operations so as a way
      #of making sure that we get the correct path of the package we use a find_path for the
      #CMakeLists.txt file for the package. Since the package has to have this file to work
      #correctly it should be guaranteed to be there.
      STRING(REGEX MATCH "[.][.]/" IS_RELATIVE_PATH ${PACKAGE_DIR})
      IF("${IS_RELATIVE_PATH}" STREQUAL "")
        SET(CPACK_SOURCE_IGNORE_FILES ${Trilinos_SOURCE_DIR}/packages/${PACKAGE_DIR} ${CPACK_SOURCE_IGNORE_FILES})
      ELSE()
        FIND_PATH(ABSOLUTE_PATH CMakeLists.txt PATHS ${Trilinos_SOURCE_DIR}/packages/${PACKAGE_DIR} NO_DEFAULT_PATH)
        IF("${ABSOLUTE_PATH}" STREQUAL "ABSOLUTE_PATH-NOTFOUND")
          MESSAGE(AUTHOR_WARNING "Relative path found for disabled package ${TRIBITS_PACKAGE} but package was missing a CMakeLists.txt file. This disabled package will likely not be excluded from a source release")
        ENDIF()
        SET(CPACK_SOURCE_IGNORE_FILES ${ABSOLUTE_PATH} ${CPACK_SOURCE_IGNORE_FILES})
      ENDIF()
    ENDIF()
  ENDFOREACH()

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("Exclude files when building source packages")
    FOREACH(item ${CPACK_SOURCE_IGNORE_FILES})
      MESSAGE(${item})
    ENDFOREACH()
  ENDIF()

  # The CPACK_RESOURCE_FILE_[LICENSE|README] files must end in one of
  # .txt .rtf .html. Copy the pertinant file to the binary directory with
  # a .txt extension. This is only the case with the PackageMaker 
  # generator, but it doesn't hurt to do it for other generators as
  # well.
  MACRO(COPY_INSTALLER_RESOURCE _varname _source _destination)
    SET("${_varname}" "${_destination}")
    IF (EXISTS "${_destination}")
      FILE(REMOVE_RECURSE "${_destination}")
    ENDIF ()
    CONFIGURE_FILE(
      "${_source}" 
      "${_destination}" 
      COPYONLY)
  ENDMACRO()
  COPY_INSTALLER_RESOURCE(Trilinos_README
    "${Trilinos_SOURCE_DIR}/README"
    "${Trilinos_BINARY_DIR}/README.txt")
  COPY_INSTALLER_RESOURCE(Trilinos_LICENSE
    "${Trilinos_SOURCE_DIR}/LICENSE"
    "${Trilinos_BINARY_DIR}/LICENSE.txt")

  SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
  SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_NAME "trilinos")
  SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  SET(CPACK_PACKAGE_VERSION "${Trilinos_VERSION}")
  SET(CPACK_RESOURCE_FILE_README "${Trilinos_README}")
  SET(CPACK_RESOURCE_FILE_LICENSE "${Trilinos_LICENSE}")
  SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
  SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")
  SET(CPACK_COMPONENTS_ALL ${Trilinos_PACKAGES} Unspecified)
  
  set(ENABLED_FLAG ON)
  set(INCLUDE_EMPTY FALSE)
  TRIBITS_GET_ENABLED_LIST( Trilinos_PACKAGES Trilinos ${ENABLED_FLAG}
    ${INCLUDE_EMPTY} ENABLED_PACKAGES NUM_ENABLED)
  string(REPLACE " " ";" ENABLED_PACKAGES "${ENABLED_PACKAGES}")
  
  #message("ENABLED PACKAGES: ${ENABLED_PACKAGES} ${NUM_ENABLED}")
  FOREACH(PKG ${ENABLED_PACKAGES})
    IF(NOT "${${PKG}_LIB_REQUIRED_DEP_PACKAGES}" STREQUAL "")
        string(TOUPPER ${PKG} UPPER_PKG)
        #message("${UPPER_PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
        SET(CPACK_COMPONENT_${UPPER_PKG}_DEPENDS ${${PKG}_LIB_REQUIRED_DEP_PACKAGES})
    ENDIF()
    #message("${PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
  ENDFOREACH()

  
  IF(WIN32)
    #Resetting the name to avoid overwriting registery keys when installing
    SET(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-${Trilinos_VERSION}")
    IF (TPL_ENABLE_MPI)
      SET(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-mpi")
    ELSE ()
      SET(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-serial")
    ENDIF()
    SET(CPACK_GENERATOR "NSIS")
    SET(CPACK_NSIS_MODIFY_PATH OFF)
  ENDIF()
  
  INCLUDE(CPack)

ENDMACRO()
