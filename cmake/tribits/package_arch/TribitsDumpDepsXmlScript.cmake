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
# CMake script file to generate an xml dependencies file for PROJECT_NAME
#
# This script is designed to be run as:
#
#   $ cmake -P -DPACKAGE_NAME=<packageName> [options] \
#       SOME_BASE_DIR/TribitsDumpDepsXmlScript.cmake
#

# A) Echo input options (must be specified with -D arguments to CMake command)

MESSAGE("PROJECT_NAME = ${PROJECT_NAME}")
MESSAGE("PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = ${${PROJECT_NAME}_TRIBITS_DIR}")
MESSAGE("${PROJECT_NAME}_EXTRA_REPOSITORIES = ${${PROJECT_NAME}_EXTRA_REPOSITORIES}")
MESSAGE("${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR = ${${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR}")
# Get the utility CMake code that can determine the dependencies

SET( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/package_arch"
  )

INCLUDE(TribitsGlobalMacros)

# Generate the dependencies file

SET(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES FALSE)
SET(${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK TRUE)
SET(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES FALSE)
TRIBITS_READ_IN_NATIVE_REPOSITORIES()
  SET(${PROJECT_NAME}_ALL_REPOSITORIES ${${PROJECT_NAME}_NATIVE_REPOSITORIES}
    ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()
