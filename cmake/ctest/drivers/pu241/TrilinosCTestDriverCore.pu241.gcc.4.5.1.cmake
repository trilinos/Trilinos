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
  SET( CTEST_PARALLEL_LEVEL 8 )

  #SET( CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind )
  #SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS )

  # We don't have TPLs for PyTrilinos, Optika, and TriKota
  # 
  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES
     PyTrilinos TriKota Optika  # We don't have TPLs for these
     Sundance Stokhos # Currently have failures and nor currently needed by CASL
     TrilinosFramework # Has 11 failing tests for some reason so disabling for now
     )
  SET(EXTRA_CONFIGURE_OPTIONS
    -DTrilinos_ENABLE_TriKota:BOOL=OFF
    -DTrilinos_ENABLE_Stokhos:BOOL=OFF
    # Allow user to override these by putting theirs at the bottom
    ${EXTRA_CONFIGURE_OPTIONS}
    )
  
  SET_DEFAULT(COMPILER_VERSION "GCC-4.5.1-IFort-12.0.4")

  SET_DEFAULT(Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE Nightly)
  
  IF (COMM_TYPE STREQUAL MPI)
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      "-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=${CTEST_SCRIPT_DIRECTORY}/gcc-4.5.1-mpi-ss-options.cmake"
      "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
      "-DTPL_ENABLE_MPI:BOOL=ON"
      ${EXTRA_CONFIGURE_OPTIONS}
    )
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      "-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=${CTEST_SCRIPT_DIRECTORY}/gcc-4.5.1-serial-ss-options.cmake"
      "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
      "-DTPL_ENABLE_MPI:BOOL=OFF"
      ${EXTRA_CONFIGURE_OPTIONS}
    )

  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()
