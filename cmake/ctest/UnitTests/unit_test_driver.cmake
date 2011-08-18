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
# General setup code:
#
# Do not modify any of this directly, use use environment variables instead!
#

#
# Include some CMake/CTest code files
#

#MESSAGE("CTEST_SCRIPT_DIRECTORY = ${CTEST_SCRIPT_DIRECTORY}")

SET( CMAKE_MODULE_PATH
  "${CTEST_SCRIPT_DIRECTORY}/.."
  "${CTEST_SCRIPT_DIRECTORY}/../../utils"
  )


#
# Includes
#

INCLUDE(TrilinosCTestDriverCore)
INCLUDE(GetLastDirName)


#
# Override some configuration variables
#

# All these can be changed by env vars
SET(CTEST_TEST_TYPE Experimental)
#SET(CTEST_DO_UPDATES FALSE)
SET(Trilinos_WARNINGS_AS_ERRORS_FLAGS "-Werror")

# Don't change these in the env!
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)
SET(CTEST_GENERATE_DEPS_XML_OUTPUT_FILE TRUE)
SET(CTEST_WIPE_CACHE FALSE)

SET(CTEST_SOURCE_DIRECTORY "${CTEST_SCRIPT_DIRECTORY}/../../..")
SET(Trilinos_DEPS_HOME_DIR "${CTEST_SOURCE_DIRECTORY}/cmake/DependencyUnitTests/MockTrilinos")
GET_FILENAME_COMPONENT(PWD . REALPATH)
SET(CTEST_BINARY_DIRECTORY "${PWD}")

SET_DEFAULT_AND_FROM_ENV(Trilinos_CTEST_COMMAND ctest)
SET(CTEST_COMMAND ${Trilinos_CTEST_COMMAND})


#
# Run the build/test/submit driver
#

TRILINOS_CTEST_DRIVER()
