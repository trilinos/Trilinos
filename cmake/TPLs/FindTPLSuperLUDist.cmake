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


TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( SuperLUDist
  REQUIRED_HEADERS "superlu_defs.h superludefs.h" supermatrix.h
  REQUIRED_LIBS_NAMES "superludist superlu_dist superlu_dist_2.0 superlu_dist_2.5"
  )

include(CheckCSourceCompiles)
include(MultilineSet)

# Versions 3.0 and later of SuperLU_DIST namespace the IterRefine_t
# enum values with "SLU_" (e.g. "SLU_DOUBLE" versus "DOUBLE" in
# previous versions).  Check which style is used so that packages like
# Amesos and Amesos2 can use the correct enum value.
FUNCTION(CHECK_SUPERLUDIST_ENUM_NAMESPACE  VARNAME)
  SET(SOURCE
  "
#include <superlu_enum_consts.h>

int main()
{
  IterRefine_t refine = SLU_DOUBLE;
  return 0;
}
"
  )

  SET(CMAKE_REQUIRED_INCLUDES ${TPL_SuperLUDist_INCLUDE_DIRS})
  CHECK_C_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()

# Version 4.0 of SuperLU_DIST changed the calling parameters of the
# LUstructInit function.  Check which is used here.
FUNCTION(CHECK_SUPERLUDIST_LUSTRUCTINIT  VARNAME)
  SET(SOURCE
  "
#include <superlu_ddefs.h>

int main()
{
  LUstruct_t lu;
  /* This will fail to compile if the 3-arg version is declared. */
  LUstructInit(10, &lu);
}
"
  )

  SET(CMAKE_REQUIRED_INCLUDES ${TPL_SuperLUDist_INCLUDE_DIRS})
  SET(CMAKE_REQUIRED_INCLUDES ${TPL_SuperLUDist_LIBRARIES})
  CHECK_C_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()

IF (TPL_ENABLE_SuperLUDist)
  CHECK_SUPERLUDIST_ENUM_NAMESPACE(HAVE_SUPERLUDIST_ENUM_NAMESPACE)
  CHECK_SUPERLUDIST_LUSTRUCTINIT(HAVE_SUPERLUDIST_LUSTRUCTINIT_2ARG)
ENDIF(TPL_ENABLE_SuperLUDist)
