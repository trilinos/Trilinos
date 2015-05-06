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
# First search for HDF5 components using standard FIND_PACKAGE(HDF5 ...)
#
TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE(HDF5  HDF5_ALLOW_PREFIND)
IF (HDF5_ALLOW_PREFIND)

  MESSAGE("-- Using FIND_PACKAGE(HDF5 ...) ...") 

  SET(HDF5_COMPNENTS C)
  IF (HDF5_REQUIRE_FORTRAN)
    LIST(APPEND HDF5_COMPNENTS Fortran)
  ENDIF()

  FIND_PACKAGE(HDF5 COMPONENTS ${HDF5_COMPNENTS})

  IF(HDF5_FOUND)
    SET(TPL_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} CACHE PATH
      "HDF5 include dirs")
    SET(TPL_HDF5_LIBRARIES ${HDF5_LIBRARIES} CACHE FILEPATH
      "HDF5 libraries")
    SET(TPL_HDF5_LIBRARY_DIRS ${HDF5_LIBRARY_DIRS} CACHE PATH
      "HDF5 library dirs")
  ENDIF()

ENDIF()

#
# Second, set up default find operation using TriBITS system in case the
# default FIND_PACKAGE(HDF5 command was not used).
#

SET(REQUIRED_HEADERS hdf5.h)
SET(REQUIRED_LIBS_NAMES hdf5)

IF (HDF5_REQUIRE_FORTRAN)
  SET(REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} hdf5_fortran)
ENDIF()

IF (TPL_ENABLE_MPI)
  SET(REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} z)
ENDIF()

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( HDF5
  REQUIRED_HEADERS ${REQUIRED_HEADERS}
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )

# NOTE: If FIND_PACKAGE(HDF5 ...) was called and was successful, then
# TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES() will just use the already set
# variables TPL_HDF5_INCLUDE_DIRS and TPL_HDF5_LIBRARIES and just print them
# out (and set some other standard varaibles).
