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


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for IBM AIX machine p90n03
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)
  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET(ENV{OBJECT_MODE} 64)
  SET(COMPILER_VERSION XL-11.01)
  SET(CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  SET(CTEST_BUILD_FLAGS "-j8 -i" )
  SET(Trilinos_REPOSITORY_LOCATION $ENV{HOME}/bench/gitroot/Trilinos.git)
  SET(Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE " ")
  SET(Trilinos_WARNINGS_AS_ERRORS_FLAGS " ")
  #SET(Trilinos_TRACK "Experimental")
  SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES PyTrilinos TriKota Optika)
  IF (COMM_TYPE STREQUAL MPI)
      SET(TPL_ENABLE_MPI ON)
  ENDIF()
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DCMAKE_C_COMPILER:FILEPATH=/usr/vacpp/bin/xlc_r"
    "-DCMAKE_CXX_COMPILER:FILEPATH=/usr/vacpp/bin/xlC_r"
    "-DCMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/xlf_r"
    "-DCMAKE_CXX_FLAGS:STRING=-qrtti=all -qstaticinline"
    "-DCMAKE_Fortran_FLAGS:STRING=-qarch=pwr7 -qtune=pwr7"
    "-DBUILD_SHARED_LIBS=ON"
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    "-DTPL_BLAS_LIBRARIES=$ENV{LAPACK_DIR}/libblas.so"
    "-DTPL_LAPACK_LIBRARIES=$ENV{LAPACK_DIR}/liblapack.so"
    "-DTPL_Pthread_LIBRARIES=-lpthread"
    )

  TRILINOS_CTEST_DRIVER()

ENDMACRO()
