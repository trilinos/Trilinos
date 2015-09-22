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
# Platform/compiler specific options for godel using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET( CTEST_BUILD_FLAGS "-j8 -i" )

  SET( CTEST_COVERAGE_COMMAND /usr/bin/gcov )
  SET( CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

  # Only turn on PyTrilinos for shared libraries
  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES PyTrilinos TriKota Optika)
  
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
    "-DBoost_INCLUDE_DIRS:PATH=/home/trilinos/tpls/boost-1.40.0"
    "-DCOVERAGE_COMMAND:FILEPATH=/usr/bin/gcov"
    "-DMEMORYCHECK_COMMAND:FILEPATH=/usr/bin/valgrind"
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE"
    "-DTPL_ENABLE_BinUtils=ON"
    )

  IF (BUILD_TYPE STREQUAL "DEBUG")
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTrilinos_ENABLE_DEBUG:BOOL=ON"
      )
  ENDIF()

  SET_DEFAULT(COMPILER_VERSION "GCC-4.1.2")
  
  IF (COMM_TYPE STREQUAL MPI)
    SET(TPL_ENABLE_MPI ON)
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      "-DMPI_BASE_DIR:PATH=/usr/local"
      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_godel_openmpi_1.2.7.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )
  
  ELSE()

    # 2009/09/02: rabartl: Let CMake pick its own compilers (fixes Intrepid problem) ...  
    #SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    #  ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
    #  "-DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++"
    #  "-DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc"
    #  "-DCMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/f77"
    #  )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_godel_gcc-4.1.2.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()
