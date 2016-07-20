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
# These options are some temporary stuff that will probably be removed later
#

# This will be modified later
MESSAGE( "Turning CTEST_DO_UPDATES off to prevent attempt to access Sandia libraries." )
SET(CTEST_DO_UPDATES FALSE)

# This will be modified later
#MESSAGE( "Turning CTEST_DO_SUBMIT off to prevent dashboard submissions temporarily." )
#SET(CTEST_DO_SUBMIT FALSE)
#SET(CTEST_SUBMIT_CDASH_SUBPROJECTS_DEPS_FILE FALSE)

#
# Main macro for running a test called RunTest
#

MACRO( RunTest parallel epetra galeri pamgen scotchandparmetis loopLOandGOSetting loopScalarSetting )

MESSAGE( "CALLING MACRO WITH ${parallel} ${epetra} ${galeri} ${pamgen} ${scotchandparmetis} ${loopLOandGOSetting} ${loopScalarSetting}" )

# Pattern is set these false and then turn on the ones specific to the loop

SET( Tpetra_INST_INT_INT_SETTING FALSE )
SET( Tpetra_INST_INT_UNSIGNED_SETTING FALSE )
SET( Tpetra_INST_INT_LONG_SETTING FALSE )
SET( Tpetra_INST_INT_LONG_LONG_SETTING FALSE )
SET( Tpetra_INST_INT_UNSIGNED_LONG_SETTING FALSE )

SET( Tpetra_INST_DOUBLE_SETTING FALSE )
SET( Tpetra_INST_FLOAT_SETTING FALSE )

IF( ${parallel} STREQUAL "TRUE" )
  SET( CUSTOM_PREFIX "par" )
  SET(COMM_TYPE MPI )
ELSE()
  SET( CUSTOM_PREFIX "ser" )
  SET(COMM_TYPE SERIAL )
ENDIF()

IF( ${loopLOandGOSetting} STREQUAL "INT_INT" )
  SET( Tpetra_INST_INT_INT_SETTING TRUE )
ELSEIF( ${loopLOandGOSetting} STREQUAL "INT_UNSIGNED" )
  SET( Tpetra_INST_INT_UNSIGNED_SETTING TRUE )
ELSEIF( ${loopLOandGOSetting} STREQUAL "INT_LONG" )
  SET( Tpetra_INST_INT_LONG_SETTING TRUE )
ELSEIF( ${loopLOandGOSetting} STREQUAL "INT_LONG_LONG" )
  SET( Tpetra_INST_INT_LONG_LONG_SETTING TRUE )
ELSEIF( ${loopLOandGOSetting} STREQUAL "INT_UNSIGNED_LONG" )
  SET( Tpetra_INST_INT_UNSIGNED_LONG_SETTING TRUE )
ENDIF()

IF( ${loopScalarSetting} STREQUAL "DOUBLE" )
  SET( Tpetra_INST_DOUBLE_SETTING TRUE )
ELSEIF( loopScalarSetting STREQUAL "FLOAT" )
  SET( Tpetra_INST_FLOAT_SETTING TRUE )
ENDIF()

# Currently just for Tech-X setup and not for normal Trilinos CDash testing
IF( DEFINED ENV{OVERRIDE_BINARY_LOCATION_ROOT} )
  SET( CUSTOM_BUILD_FOLDER_NAME "${CUSTOM_PREFIX}-LOGO-${loopLOandGOSetting}-Scalar-${loopScalarSetting}-Epetra-${epetra}-ScotchParmetis-${scotchandparmetis}-Gal-${galeri}-Pamgen-${pamgen}" )
  SET( CTEST_BINARY_DIRECTORY "$ENV{OVERRIDE_BINARY_LOCATION_ROOT}/${CUSTOM_BUILD_FOLDER_NAME}" )
  SET( CTEST_SOURCE_DIRECTORY "${CTEST_SCRIPT_DIRECTORY}/../../../../../Trilinos" )
  MESSAGE( "Custom build location setup set CTEST_BINARY_DIRECTORY to ${CTEST_BINARY_DIRECTORY}")
  MESSAGE( "Custom build location setup set CTEST_SOURCE_DIRECTORY to ${CTEST_SOURCE_DIRECTORY}")
  # Not necessary but convenient - you can comment out the actual run TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER() the end and build all the directories quickly for testing setup
  file(MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY})
ENDIF()


#
# Set the options specific to this build case - must be BEFORE including the jenkins file
#

SET(BUILD_TYPE DEBUG)
#SET(BUILD_DIR_NAME)
SET(CTEST_PARALLEL_LEVEL 8)
SET(CTEST_TEST_TYPE Experimental)
SET(CTEST_TEST_TIMEOUT 900)
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

# Generic driver script for building & testing Trilinos with Jenkins
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.generic.jenkins.cmake")

SET(Trilinos_PACKAGES Zoltan2)

SET(EXTRA_CONFIGURE_OPTIONS
"-DZoltan2_ENABLE_Experimental:BOOL=$ENV{JENKINS_Zoltan2_ENABLE_Experimental}"
"-DTeuchos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=$ENV{JENKINS_Teuchos_ENABLE_EXPLICIT_INSTANTIATION}"

# Not currently implemented
# "-DTpetra_INST_COMPLEX_DOUBLE:BOOL=$ENV{JENKINS_Tpetra_INST_COMPLEX_DOUBLE}"
# "-DTpetra_INST_COMPLEX_FLOAT:BOOL=$ENV{JENKINS_Tpetra_INST_COMPLEX_FLOAT}"

# Need to decide about this
"-DTeuchos_ENABLE_LONG_LONG_INT:BOOL=TRUE"

"-DTpetra_INST_INT_INT:BOOL=${Tpetra_INST_INT_INT_SETTING}"
"-DTpetra_INST_INT_UNSIGNED:BOOL=${Tpetra_INST_INT_UNSIGNED_SETTING}"
"-DTpetra_INST_INT_LONG:BOOL=${Tpetra_INST_INT_LONG_SETTING}"
"-DTpetra_INST_INT_LONG_LONG:BOOL=${Tpetra_INST_INT_LONG_LONG_SETTING}"
"-DTpetra_INST_INT_UNSIGNED_LONG:BOOL=${Tpetra_INST_INT_UNSIGNED_LONG_SETTING}"

"-DTpetra_INST_DOUBLE:BOOL=${Tpetra_INST_DOUBLE_SETTING}"
"-DTpetra_INST_FLOAT:BOOL=${Tpetra_INST_FLOAT_SETTING}"

"-DTPL_ENABLE_Scotch:BOOL=${scotchandparmetis}"
"-DTPL_ENABLE_ParMETIS:BOOL=${scotchandparmetis}"

"-DTrilinos_ENABLE_Galeri:BOOL=${galeri}"
"-DTrilinos_ENABLE_Pamgen:BOOL=${pamgen}"
"-DTrilinos_ENABLE_Epetra:BOOL=${epetra}"
)

MESSAGE( "EXTRA_CONFIGURE_OPTIONS was set to ${EXTRA_CONFIGURE_OPTIONS}" )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

ENDMACRO( RunTest )


MACRO( RunFullLoop )

  MESSAGE( "Begin the loop of different tests." )

  FOREACH( parallel TRUE )                       # Right now just parallel
  FOREACH( epetra TRUE )                         # Right now just epetra on
  FOREACH( galeri FALSE )                        # Right now just galeri off
  FOREACH( pamgen FALSE )                        # Right now just pamgen off
  FOREACH( scotchandparmetis FALSE TRUE )        # parmetis+scotch off (0) or on (1) - Right now both together
  FOREACH( localglobal RANGE INT_INT INT_UNSIGNED INT_LONG INT_LONG_LONG INT_UNSIGNED_LONG )               # 5 settings for LO and GO
  FOREACH( scalar DOUBLE FLOAT )                 # double or float - need to check on int status

    RunTest( ${parallel} ${epetra} ${galeri} ${pamgen} ${scotchandparmetis} ${localglobal} ${scalar} )

  ENDFOREACH( loopScalarSetting )
  ENDFOREACH( loopLOandGOSetting )
  ENDFOREACH( parmetisScotchSetting )
  ENDFOREACH( pamgenSetting )
  ENDFOREACH( galeriSetting )
  ENDFOREACH( epetraSetting )
  ENDFOREACH( parallelSetting )

ENDMACRO( RunFullLoop )




#        PARALLEL    EPETRA   GALERI      PAMGEN      SCOTCH&PARMETIS       LO-GO       Scalar
#RunTest(  TRUE        TRUE     FALSE      FALSE          FALSE            INT_INT       DOUBLE   )
#RunTest(  TRUE        TRUE     FALSE      FALSE          FALSE            INT_INT       FLOAT    )
#RunTest(  TRUE        TRUE     FALSE      FALSE          FALSE            INT_LONG      FLOAT    )
#RunTest(  TRUE        TRUE     TRUE       FALSE          TRUE             INT_LONG      FLOAT    )


# This will take hours
RunFullLoop()
