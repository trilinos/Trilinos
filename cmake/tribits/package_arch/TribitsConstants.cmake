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

# Define the TriBITS minimum required CMake version
SET(TRIBITS_CMAKE_MINIMUM_REQUIRED 2.7)

# File names for TriBITS system

SET(${PROJECT_NAME}_PACKAGES_FILE_NAME PackagesList.cmake)

SET(${PROJECT_NAME}_TPLS_FILE_NAME TPLsList.cmake)

SET(${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME ExtraRepositoriesList.cmake)

SET(${PROJECT_NAME}_EXTRA_PACKAGES_FILE_NAME PackagesList.cmake)

SET(${PROJECT_NAME}_EXTRA_TPLS_FILE_NAME TPLsList.cmake)

# Directories relative to the TriBITS base directory

SET(TRIBITS_PYTHON_SCRIPTS_DIR "python")

SET(TRIBITS_CMAKE_UTILS_DIR "utils")

SET(TRIBITS_CMAKE_PACKAGE_ARCH_DIR "package_arch")

SET(TRIBITS_CMAKE_INSTALLATION_FILES_DIR "installation")

SET(TRIBITS_MOCK_TRILINOS_DIR "package_arch/UnitTests/MockTrilinos")

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
