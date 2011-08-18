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

INCLUDE(CheckIncludeFileCXX)

IF(MSVC)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE 
    -D_CRT_NONSTDC_NO_DEPRECATE  -D_SCL_SECURE_NO_WARNINGS)
  INCLUDE_DIRECTORIES(${Trilinos_SOURCE_DIR}/commonTools/WinInterface/include)
  # find the CLAPACK built by CMake on the machine for MSVC
  # if found it will set the BLAS and LAPACK libraries
  FIND_PACKAGE(CLAPACK 3.2.1 NO_MODULE)
  IF(CLAPACK_FOUND)
    SET(TPL_BLAS_LIBRARIES blas CACHE INTERNAL "")
    SET(TPL_LAPACK_LIBRARIES lapack CACHE INTERNAL "")
  ENDIF()
ENDIF()

IF (WIN32 AND NOT CYGWIN)
  SET(NATIVE_MS_WINDOWS TRUE)
ELSE()
  SET(NATIVE_MS_WINDOWS FALSE)
ENDIF()

# Probe for non-standard headers

IF (Trilinos_ENABLE_CXX)
  CHECK_INCLUDE_FILE_CXX(sys/time.h HAVE_SYS_TIME_H)
  CHECK_INCLUDE_FILE_CXX(time.h HAVE_TIME_H)
  CHECK_INCLUDE_FILE_CXX(stdint.h HAVE_STDINT_H)
  CHECK_INCLUDE_FILE_CXX(inttypes.h HAVE_INTTYPES_H)
ENDIF()

SET(HAVE_ALGORITHM TRUE)
SET(HAVE_CASSERT TRUE)
SET(HAVE_CCTYPE TRUE)
SET(HAVE_CERRNO TRUE)
SET(HAVE_CLIMITS TRUE)
SET(HAVE_CMATH TRUE)
SET(HAVE_COMPLEX TRUE)
SET(HAVE_CSTDARG TRUE)
SET(HAVE_CSTDIO TRUE)
SET(HAVE_CSTDLIB TRUE)
SET(HAVE_CSTRING TRUE)
SET(HAVE_IOMANIP TRUE)
SET(HAVE_IOSTREAM TRUE)
SET(HAVE_ITERATOR TRUE)
SET(HAVE_LIST TRUE)
SET(HAVE_MAP TRUE)
SET(HAVE_MEMORY TRUE)
SET(HAVE_MUTABLE TRUE)
SET(HAVE_NAMESPACES TRUE)
SET(HAVE_NEW_FOR_SCOPING TRUE)
SET(HAVE_NUMERIC TRUE)
SET(HAVE_NUMERIC_LIMITS TRUE)
SET(HAVE_POW TRUE)
SET(HAVE_SET TRUE)
SET(HAVE_SSTREAM TRUE)
SET(HAVE_FSTREAM TRUE)
SET(HAVE_STDEXCEPT TRUE)
SET(HAVE_STRING TRUE)
SET(HAVE_VECTOR TRUE)

# 2008/12/20: rabartl: Above: All of these defines should be removed
# because we decided that we were going to assume that all compilers
# have these C++98 standard features.  We will deal with cases where
# this is not true but we should not assume the worst right from the
# beginning.

# Find Perl

FIND_PACKAGE(Perl)

# Do Fortran stuff

INCLUDE(TrilinosFortranMangling)

# Get BLAS name mangling
 
INCLUDE(TrilinosBLASMangling)

# Determine C++-0x supported features

IF (Trilinos_ENABLE_CXX11)
  INCLUDE(TrilinosCXX11Support)
  CHECK_CXX11_SUPPORT(Trilinos_ENABLE_CXX11)
  MESSAGE("-- Trilinos_ENABLE_CXX11=${Trilinos_ENABLE_CXX11}")
ENDIF()

# Set up some MPI info

IF (TPL_ENABLE_MPI)
  SET(HAVE_MPI TRUE)
ELSE()
  SET(HAVE_MPI FALSE)
ENDIF()

# OpenMP isn't really a TPL because support is built into the compiler.

IF(Trilinos_ENABLE_OpenMP)
  INCLUDE(FindOpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    # FindOpenMP.cmake doesn't find Fortran flags.  Mike H said this is safe.
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
  ELSE()
    MESSAGE(FATAL_ERROR "Could not find OpenMP, try setting OpenMP_C_FLAGS and OpenMP_CXX_FLAGS directly")
  ENDIF(OPENMP_FOUND)
ENDIF(Trilinos_ENABLE_OpenMP)

# Check if we need the math library or not and find the right one
IF (NOT NATIVE_MS_WINDOWS)
  INCLUDE(MathLibraryNeeded)
ENDIF()

# Check for isnan and isinf support
IF (${PROJECT_NAME}_ENABLE_CXX)
  INCLUDE(FiniteValue)
ENDIF()

# Check for Doxygen/dot - We can use variables set in this check to
# enable/disable the grapical dependency graphs in doxygen Doxyfiles.
INCLUDE(FindDoxygen)

