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

# Define the TriBITS minimum required CMake version
SET(TRIBITS_CMAKE_MINIMUM_REQUIRED 2.8.11)

MACRO(TRIBITS_ASESRT_MINIMUM_CMAKE_VERSION)

  IF (CMAKE_VERSION VERSION_LESS ${TRIBITS_CMAKE_MINIMUM_REQUIRED})
    MESSAGE(FATAL_ERROR "Error, TriBiTS must have version"
        " ${TRIBITS_CMAKE_MINIMUM_REQUIRED} or higher!")
  ENDIF()
  
  IF(NOT CMAKE_VERSION VERSION_LESS 3.3)
    CMAKE_POLICY(VERSION 3.3)
  ELSE()
    CMAKE_POLICY(VERSION ${CMAKE_VERSION})
  ENDIF()

ENDMACRO()

# File names for TriBITS system

SET(${PROJECT_NAME}_PACKAGES_FILE_NAME PackagesList.cmake)

SET(${PROJECT_NAME}_TPLS_FILE_NAME TPLsList.cmake)

SET(${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME ExtraRepositoriesList.cmake)

SET(${PROJECT_NAME}_EXTRA_PACKAGES_FILE_NAME PackagesList.cmake)

SET(${PROJECT_NAME}_REPO_VERSION_FILE_NAME ${PROJECT_NAME}RepoVersion.txt)

SET(${PROJECT_NAME}_EXTRA_TPLS_FILE_NAME TPLsList.cmake)

# Directories relative to the TriBITS base directory

SET(TRIBITS_PYTHON_SCRIPTS_DIR "python_utils")

SET(TRIBITS_CI_SUPPORT_DIR "ci_support")

SET(TRIBITS_CTEST_DRIVER_DIR "ctest_driver")

SET(TRIBITS_CMAKE_UTILS_DIR "core/utils")

SET(TRIBITS_CMAKE_PACKAGE_ARCH_DIR "core/package_arch")

SET(TRIBITS_CMAKE_INSTALLATION_FILES_DIR "core/installation")

# Files and directories related to the specific project

SET(${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME ${PROJECT_NAME}PackageDependencies.xml)

SET(${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME CDashSubprojectDependencies.xml)

SET(${PROJECT_NAME}_PACKAGE_DEPS_TABLE_HTML_FILE_NAME ${PROJECT_NAME}PackageDependenciesTable.html)

SET(${PROJECT_NAME}_PACKAGE_DEPS_FILES_DIR "cmake/dependencies")

# Other stuff

IF(WIN32)
  #Apparently FIND_PROGRAM looks for an exact match of the file name.
  #So even though "git clone ..." is valid to use on windows we need to give the
  #full name of the command we want to run.
  SET(GIT_NAME git.cmd)
ELSE(WIN32)
  SET(GIT_NAME git)
ENDIF(WIN32)
