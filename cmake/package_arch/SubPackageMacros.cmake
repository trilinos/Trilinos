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

INCLUDE(PackageMacros)


#
# Define a subpackage.
#
# Once called, the following varibles are in scope:
#
#   PARENT_PACKAGE_NAME: The name of the parent package
#
#   SUBPACKAGE_NAME: The local name of the subpackage (does not contain the
#   parent package name).
#
#   SUBPACKAGE_FULLNAME: The full project-level name of the subpackage (which
#   includes the parent package name at the beginning).
#
#   PACKAGE_NAME: Inside the subpackage, the same as SUBPACKAGE_FULLNAME.
#
MACRO(SUBPACKAGE SUBPACKAGE_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nSUBPACKAGE: ${SUBPACKAGE_NAME_IN}")
  ENDIF()

  IF (NOT ${SUBPACKAGE_NAME_IN} STREQUAL ${SUBPACKAGE_NAME})
    MESSAGE(FATAL_ERROR "Error, the package-defined subpackage name"
      " '${SUBPACKAGE_NAME_IN}' is not the same as the subpackage name"
      " '${SUBPACKAGE_NAME}' defined in the parent packages's"
      " Dependencies.cmake file")
  ENDIF()

  # To provide context for various macros
  SET(PACKAGE_NAME ${SUBPACKAGE_FULLNAME})

  SET(PARENT_PACKAGE_SOURCE_DIR "${PACKAGE_SOURCE_DIR}")
  SET(PARENT_PACKAGE_BINARY_DIR "${PACKAGE_BINARY_DIR}")

  # Now override the package-like varibles
  PACKAGE_SET_COMMON_VARS(${SUBPACKAGE_FULLNAME})
  PACKAGE_DEFINE_LINKAGE_VARS(${SUBPACKAGE_FULLNAME})

ENDMACRO()


#
# Postprocess after defining a subpackage
#

MACRO(SUBPACKAGE_POSTPROCESS)
  PACKAGE_POSTPROCESS()
ENDMACRO()
