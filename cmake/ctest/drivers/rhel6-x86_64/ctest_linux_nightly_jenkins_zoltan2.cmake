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

# This will be modified later
MESSAGE( "Turning CTEST_DO_UPDATES off to prevent attempt to access Sandia libraries." )
SET(CTEST_DO_UPDATES FALSE)

# This will be modified later
#MESSAGE( "Turning CTEST_DO_SUBMIT off to prevent dashboard submissions temporarily." )
#SET(CTEST_DO_SUBMIT FALSE)
#SET(CTEST_SUBMIT_CDASH_SUBPROJECTS_DEPS_FILE FALSE)

MACRO( RunTest testName )

MESSAGE( "Running subtest: ${testName}" )
MESSAGE( "EXTRA_CONFIGURE_OPTIONS: ${EXTRA_CONFIGURE_OPTIONS}" )

# Two hacks that need to be resolved
# I did this so the setups for these tests would be 'natural' config settings - not sure why the orignal format worked with BUILD_TYPE and COMM_TYPE
SET(COMM_TYPE SERIAL )
IF( "${EXTRA_CONFIGURE_OPTIONS}" MATCHES "-DTPL_ENABLE_MPI=ON" OR 
  "${EXTRA_CONFIGURE_OPTIONS}" MATCHES "-DTPL_ENABLE_MPI=TRUE" OR
  "${EXTRA_CONFIGURE_OPTIONS}" MATCHES "-DTPL_ENABLE_MPI:BOOL=ON" OR
  "${EXTRA_CONFIGURE_OPTIONS}" MATCHES "-DTPL_ENABLE_MPI:BOOL=TRUE" )
    SET(COMM_TYPE MPI )
ENDIF()
SET(BUILD_TYPE RELEASE)
IF( "${EXTRA_CONFIGURE_OPTIONS}" MATCHES "-DCMAKE_BUILD_TYPE:STRING=DEBUG" OR
  "${EXTRA_CONFIGURE_OPTIONS}" MATCHES "-DCMAKE_BUILD_TYPE=DEBUG" )
    SET(BUILD_TYPE DEBUG)
ENDIF()

# Currently just for Tech-X setup and not for normal Trilinos CDash testing
IF( DEFINED ENV{OVERRIDE_BINARY_LOCATION_ROOT} )
  SET( CTEST_BINARY_DIRECTORY "$ENV{OVERRIDE_BINARY_LOCATION_ROOT}/${testName}" )
  SET( CTEST_SOURCE_DIRECTORY "${CTEST_SCRIPT_DIRECTORY}/../../../../../Trilinos" )
  MESSAGE( "Custom build location setup set CTEST_BINARY_DIRECTORY to ${CTEST_BINARY_DIRECTORY}")
  MESSAGE( "Custom build location setup set CTEST_SOURCE_DIRECTORY to ${CTEST_SOURCE_DIRECTORY}")
  # Not necessary but convenient - you can comment out the actual run TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER() the end and build all the directories quickly for testing setup
  file(MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY})
ENDIF()


#SET(BUILD_DIR_NAME)
SET(CTEST_PARALLEL_LEVEL 8)
SET(CTEST_TEST_TYPE Experimental)
SET(CTEST_TEST_TIMEOUT 900)
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.generic.jenkins.cmake")
SET(Trilinos_PACKAGES Zoltan2)

MESSAGE( "EXTRA_CONFIGURE_OPTIONS was set to ${EXTRA_CONFIGURE_OPTIONS}" )
TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

ENDMACRO( RunTest )

#
# These are the specific tests we want to run
#

SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=ON"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_DOUBLE:BOOL=TRUE"
  "-DTpetra_INST_INT_INT:BOOL=TRUE"
)
RunTest( par-int-int-double )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=ON"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_FLOAT:BOOL=TRUE"
  "-DTpetra_INST_INT_INT:BOOL=TRUE"
)
RunTest( par-int-int-float )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=ON"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_FLOAT:BOOL=TRUE"
  "-DTpetra_INST_INT_UNSIGNED:BOOL=TRUE"
)
RunTest( par-int-unsigned-float )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=ON"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_FLOAT:BOOL=TRUE"
  "-DTpetra_INST_INT_LONG:BOOL=TRUE"
)
RunTest( par-int-long-float )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=ON"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_FLOAT:BOOL=TRUE"
  "-DTpetra_INST_INT_UNSIGNED_LONG:BOOL=TRUE"
)
RunTest( par-int-unsignedlong-float )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=ON"
  "-DCMAKE_CXX_FLAGS=-DTEST_STK_DATA_TYPES"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
)
RunTest( par-stk )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=ON"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
)
RunTest( par-stk-off )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=ON"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_DOUBLE:BOOL=TRUE"
  "-DTeuchos_ENABLE_LONG_LONG_INT:BOOL=TRUE"
  "-DTpetra_INST_INT_LONG_LONG:BOOL=TRUE"
)
RunTest( par-int-longlong-double )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=ON"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_DOUBLE:BOOL=TRUE"
  "-DTpetra_INST_INT_INT:BOOL=TRUE"
  "-DTPL_ENABLE_Scotch:BOOL=TRUE"
  "-DTPL_ENABLE_ParMETIS:BOOL=TRUE"
)
RunTest( par-int-int-double-scotch-parmetis )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_DOUBLE:BOOL=TRUE"
  "-DTpetra_INST_INT_INT:BOOL=TRUE"
  "-DTrilinos_ENABLE_Galeri:BOOL=TRUE"
)
RunTest( par-int-int-double-galeri )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_DOUBLE:BOOL=TRUE"
  "-DTpetra_INST_INT_INT:BOOL=TRUE"
  "-DTrilinos_ENABLE_Pamgen:BOOL=TRUE"
)
RunTest( par-int-int-double-pamgen )


SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
  "-DTPL_ENABLE_MPI=OFF"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=TRUE"
  "-DTrilinos_ENABLE_Epetra:BOOL=TRUE"
  "-DTpetra_INST_DOUBLE:BOOL=TRUE"
  "-DTpetra_INST_INT_INT:BOOL=TRUE"
)
RunTest( ser-int-int-double )
