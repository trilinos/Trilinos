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
# finding LAPACK.  These are used in case find_package(LAPACK ...) is not called
# or does not find LAPACK.  Also, these variables need to be non-null in order
# to trigger the right behavior in the function
# tribits_tpl_find_include_dirs_and_libraries().
#
set(REQUIRED_LIBS_NAMES "lapack lapack_win32")

#
# Second, search for LAPACK components (if allowed) using the standard
# find_package(LAPACK ...).
#
tribits_tpl_allow_pre_find_package(LAPACK  LAPACK_ALLOW_PREFIND)
if (LAPACK_ALLOW_PREFIND)

  message("-- Using find_package(LAPACK ...) ...")

  find_package(LAPACK)

  if (LAPACK_FOUND)
    # Tell TriBITS that we found LAPACK and there no need to look any further!
    set(TPL_LAPACK_INCLUDE_DIRS "" CACHE PATH
      "LAPACK include dirs")
    set(TPL_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE FILEPATH
      "LAPACK libraries")
    set(TPL_LAPACK_LIBRARY_DIRS "" CACHE PATH
      "LAPACK library dirs")
  endif()

endif()

if (MSVC AND NOT
    (LAPACK_LIBRARY_DIRS  OR
     (NOT "${LAPACK_LIBRARY_NAMES}" STREQUAL "lapack lapack_win32" AND
      NOT "${LAPACK_LIBRARY_NAMES}" STREQUAL "") OR
     LAPACK_INCLUDE_DIRS  OR
     LAPACK_INCLUDE_NAMES OR
     (NOT "${TPL_LAPACK_LIBRARIES}" STREQUAL "lapack" AND
      NOT "${TPL_LAPACK_LIBRARIES}" STREQUAL "") OR
     TPL_LAPACK_INCLUDE_DIRS)
   )
  if(CLAPACK_FOUND)
    advanced_set(TPL_LAPACK_LIBRARIES lapack
        CACHE FILEPATH "Set from MSVC CLAPACK specialization")
  endif()
endif()

#
# Third, call tribits_tpl_find_include_dirs_and_libraries()
#
tribits_tpl_find_include_dirs_and_libraries( LAPACK
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )
# NOTE: If find_package(LAPACK ...) was called and successfully found LAPACK, then
# tribits_tpl_find_include_dirs_and_libraries() will use the already-set
# variables TPL_LAPACK_INCLUDE_DIRS and TPL_LAPACK_LIBRARIES and then print them
# out (and set some other standard variables as well).  This is the final
