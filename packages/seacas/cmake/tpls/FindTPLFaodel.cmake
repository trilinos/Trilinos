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
# First, set up the variables for the (backward-compatible) TriBITS way of
# finding Faodel.  These are used in case FIND_PACKAGE(Faodel ...) is not
# called or does not find Faodel.  Also, these variables need to be non-null
# in order to trigger the right behavior in the function
# TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES().
#
SET(REQUIRED_HEADERS "kelpie/Kelpie.hh")
SET(REQUIRED_LIBS_NAMES kelpie dirman opbox lunasa nnti faodel-services whookie faodel-common sbl tcmalloc spinlock)

#
# Second, search for Faodel components (if allowed) using the standard
# FIND_PACKAGE(Faodel ...).
#
TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE(Faodel Faodel_ALLOW_PREFIND)
IF (Faodel_ALLOW_PREFIND)

  MESSAGE(STATUS "Using FIND_PACKAGE(Faodel...) ...")

  SET(CMAKE_MODULE_PATH
    "${CMAKE_MODULE_PATH}"
    "${CMAKE_CURRENT_LIST_DIR}/find_modules"
    "${CMAKE_CURRENT_LIST_DIR}/utils"
    )

  find_package(Faodel REQUIRED CONFIG
    COMPONENTS
    tcmalloc spinlock sbl common whookie services nnti lunasa opbox dirman kelpie
    )

  IF (Faodel_FOUND)
    set(TPL_Faodel_INCLUDE_DIRS ${Faodel_INCLUDE_DIR} CACHE PATH
      "List of semi-colon separated list of directories containing Faodel header files")
    set(TPL_Faodel_LIBRARY_DIRS ${Faodel_LIBRARY_DIR} CACHE PATH
      "List of semi-colon separated list of directories containing Faodel library files")
    set(TPL_Faodel_LIBRARIES Faodel::kelpie CACHE PATH
      "List of semi-colon separated (full) paths to the Faodel libraries")
  ENDIF()

ENDIF()

#
# Third, call TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()
#
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Faodel
  REQUIRED_HEADERS ${REQUIRED_HEADERS}
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES})
