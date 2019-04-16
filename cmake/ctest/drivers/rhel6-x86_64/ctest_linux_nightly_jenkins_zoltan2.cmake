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

# For cdash testing setup we disabled updates
IF( DEFINED ENV{CUSTOM_ENV_VARIABLES} )
  SET(CTEST_DO_UPDATES FALSE)
ENDIF()

# Turn these on to disable Cdash testing
# SET(CTEST_DO_SUBMIT FALSE)
# SET(CTEST_SUBMIT_CDASH_SUBPROJECTS_DEPS_FILE FALSE)

# This probably gets changed in the final setup
SET( CTEST_BUILD_NAME $ENV{JENKINS_CTEST_BUILD_NAME} )

# Currently just for Tech-X setup and not for normal Trilinos CDash testing
# We had to change the target build directories
IF( DEFINED ENV{CUSTOM_ENV_VARIABLES} )
  SET( CTEST_BINARY_DIRECTORY "$ENV{OVERRIDE_BINARY_LOCATION_ROOT}/${CTEST_BUILD_NAME}" )
  SET( CTEST_SOURCE_DIRECTORY "${CTEST_SCRIPT_DIRECTORY}/../../../../../Trilinos" )
ENDIF()

SET(CTEST_PARALLEL_LEVEL 8)
SET(CTEST_TEST_TYPE Experimental)
SET(CTEST_TEST_TIMEOUT 900)
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

# Set up the Jenkins parameters
SET(COMM_TYPE $ENV{JENKINS_COMM_TYPE} )
SET(BUILD_TYPE $ENV{JENKINS_BUILD_TYPE} )

# Clear out the base variables - this is probably not necessary anymore
SET( Trilinos_PACKAGES Zoltan2 )
Set( EXTRA_EXCLUDE_PACKAGES "" ) # Start empty
Set( EXTRA_CONFIGURE_OPTIONS "" ) # Start empty

# Macro helper to write a Bool value but skip if string is blank
MACRO( setConfig jenkinsValue trilinosSetting )
  #  MESSAGE( "setConfig called with ${jenkinsValue} ${trilinosSetting}" )
  IF( NOT ${jenkinsValue} STREQUAL "BLANK" OR ${jenkinsValue} STREQUAL "" )
    SET( EXTRA_CONFIGURE_OPTIONS ${EXTRA_CONFIGURE_OPTIONS} "-D${trilinosSetting}=${jenkinsValue}" )
  ENDIF()
ENDMACRO()

# Macro helper to setup a package - currently investigating 4 options
# "BLANK" or "" - Do nothing - no entry added
# "EXCLUDE" - Add to EXTRA_EXCLUDE_PACKAGES and set config setting FALSE
# "INCLUDE" - Add to Trilinos_PACKAGES and set config setting FALSE
# Set based on Bool value - which can be TRUE, FALSE
MACRO( checkPackage jenkinsValue trilinosSetting packageName )
  #  MESSAGE( "checkPackage called with ${jenkinsValue} ${trilinosSetting} ${packageName}" )
  IF( ${jenkinsValue} STREQUAL "BLANK" OR ${jenkinsValue} STREQUAL ""  )
    # Do nothing
  ELSEIF( ${jenkinsValue} STREQUAL "EXCLUDE" )
    # Add to exclude and turn off
    SET(EXTRA_EXCLUDE_PACKAGES ${EXTRA_EXCLUDE_PACKAGES} ${packageName})
    SET( EXTRA_CONFIGURE_OPTIONS ${EXTRA_CONFIGURE_OPTIONS} "-D${trilinosSetting}=FALSE" )
  ELSEIF( ${jenkinsValue} STREQUAL "INCLUDE" )
    # Add to Trilinos_PACKAGES and turn on
    SET(Trilinos_PACKAGES ${Trilinos_PACKAGES} ${packageName})
    SET( EXTRA_CONFIGURE_OPTIONS ${EXTRA_CONFIGURE_OPTIONS} "-D${trilinosSetting}=TRUE" )
  ELSE() # Set TRUE or FALSE based on setting - can be ON or OFF as well
    SET( EXTRA_CONFIGURE_OPTIONS ${EXTRA_CONFIGURE_OPTIONS} "-D${trilinosSetting}=${jenkinsValue}" )
  ENDIF()
ENDMACRO()

# Turning this off has caused difficulties as the code system hard codes it to true
# In our current working system we add this false but also add optionals to
# the EXTRA_EXCLUDE_PACKAGES if we don't want them to build in the first place
setConfig( $ENV{JENKINS_ENABLE_ALL_OPTIONAL_PACKAGES} Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES )

setConfig( $ENV{JENKINS_ENABLE_EXPLICIT_INSTANTIATION}  Trilinos_ENABLE_EXPLICIT_INSTANTIATION )

# setCong skips the addition of the config element if the string is blank: ""
setConfig( $ENV{JENKINS_TPL_ENABLE_PARMETIS}            TPL_ENABLE_ParMETIS )
setConfig( $ENV{JENKINS_TPL_ENABLE_SCOTCH}              TPL_ENABLE_Scotch )

# This exists so we can handle the STK_DATA_TYPES define used in Zoltan2
setConfig( $ENV{JENKINS_CMAKE_CXX_FLAGS} "CMAKE_CXX_FLAGS" )

# These settings write for TRUE or FALSE - we want to turn them off if FALSE
setConfig( $ENV{JENKINS_Tpetra_INST_DOUBLE}             Tpetra_INST_DOUBLE )
setConfig( $ENV{JENKINS_Tpetra_INST_FLOAT}              Tpetra_INST_FLOAT )
setConfig( $ENV{JENKINS_Tpetra_INST_INT_INT}            Tpetra_INST_INT_INT )
setConfig( $ENV{JENKINS_Tpetra_INST_INT_UNSIGNED}       Tpetra_INST_INT_UNSIGNED )
setConfig( $ENV{JENKINS_Tpetra_INST_INT_UNSIGNED_LONG}  Tpetra_INST_INT_UNSIGNED_LONG )
setConfig( $ENV{JENKINS_Tpetra_INST_INT_LONG}           Tpetra_INST_INT_LONG )
setConfig( $ENV{JENKINS_Tpetra_INST_INT_LONG_LONG}      Tpetra_INST_INT_LONG_LONG )

# Special handling of the packages for cdash - handles Blank, TRUE, FALSE, EXCLUDE, or INCLUDE
checkPackage( $ENV{JENKINS_USE_EPETRA}      Trilinos_ENABLE_Epetra    Epetra    )
checkPackage( $ENV{JENKINS_USE_THYRA}       Trilinos_ENABLE_Thyra     Thyra     )
checkPackage( $ENV{JENKINS_USE_XPETRA}      Trilinos_ENABLE_Xpetra    Xpetra    )
checkPackage( $ENV{JENKINS_USE_GALERI}      Trilinos_ENABLE_Galeri    Galeri    )
checkPackage( $ENV{JENKINS_USE_PAMGEN}      Trilinos_ENABLE_Pamgen    Pamgen    )
checkPackage( $ENV{JENKINS_USE_EPETRAEXT}   Trilinos_ENABLE_EpetraExt EpetraExt )
checkPackage( $ENV{JENKINS_USE_TRIUTILS}    Trilinos_ENABLE_TriUtils  TriUtils  )
checkPackage( $ENV{JENKINS_USE_ZOLTAN}      Trilinos_ENABLE_Zoltan    Zoltan    )
checkPackage( $ENV{JENKINS_USE_RTOP}        Trilinos_ENABLE_RTOp      RTOp      )

# Skipping this will cause problems if Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES is OFF
checkPackage( TRUE      Trilinos_ENABLE_Zoltan2 Zoltan2 )

# MESSAGE( "Zoltan2 test: ${CTEST_BUILD_NAME} setup:" )
# MESSAGE( "EXTRA_CONFIGURE_OPTIONS set to: ${EXTRA_CONFIGURE_OPTIONS}" )
# MESSAGE( "EXTRA_EXCLUDE_PACKAGES set to ${EXTRA_EXCLUDE_PACKAGES}" )
# MESSAGE( "Trilinos_PACKAGES set to ${Trilinos_PACKAGES}" )

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.generic.jenkins.cmake")

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

