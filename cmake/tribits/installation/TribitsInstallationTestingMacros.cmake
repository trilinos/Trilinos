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

FUNCTION(FIND_PROJECT_INSTALL)
  IF(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("Searching for ${PROJECT_NAME} installation at ${${PROJECT_NAME}_INSTALLATION_DIR}/include")
    ENDIF()
    FIND_PACKAGE(${PROJECT_NAME} REQUIRED HINTS ${${PROJECT_NAME}_INSTALLATION_DIR})

    IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("Found ${PROJECT_NAME} installation version ${${PROJECT_NAME}_VERSION} at ${${PROJECT_NAME}_DIR}")
    ENDIF()

    #renaming some of the variables so that they do not clash with the variables of the same name.
    SET(${PROJECT_NAME}_INSTALLATION_VERSION           ${${PROJECT_NAME}_VERSION}           PARENT_SCOPE)
    SET(${PROJECT_NAME}_INSTALLATION_INCLUDE_DIRS      ${${PROJECT_NAME}_INCLUDE_DIRS}      PARENT_SCOPE)
    SET(${PROJECT_NAME}_INSTALLATION_LIBRARY_DIRS      ${${PROJECT_NAME}_LIBRARY_DIRS}      PARENT_SCOPE)
    SET(${PROJECT_NAME}_INSTALLATION_LIBRARIES         ${${PROJECT_NAME}_LIBRARIES}         PARENT_SCOPE)
    SET(${PROJECT_NAME}_INSTALLATION_PACKAGE_LIST      ${${PROJECT_NAME}_PACKAGE_LIST}      PARENT_SCOPE)
    SET(${PROJECT_NAME}_INSTALLATION_BUILD_SHARED_LIBS ${${PROJECT_NAME}_BUILD_SHARED_LIBS} PARENT_SCOPE)
    SET(${PROJECT_NAME}_INSTALLATION_TPL_INCLUDE_DIRS  ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}  PARENT_SCOPE)
    SET(${PROJECT_NAME}_INSTALLATION_TPL_LIBRARY_DIRS  ${${PROJECT_NAME}_TPL_LIBRARY_DIRS}  PARENT_SCOPE)
    SET(${PROJECT_NAME}_INSTALLATION_TPL_LIBRARIES     ${${PROJECT_NAME}_TPL_LIBRARIES}     PARENT_SCOPE)

    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGE_LIST})
      SET(${TRIBITS_PACKAGE}_INSTALLATION_INCLUDE_DIRS     ${${TRIBITS_PACKAGE}_INCLUDE_DIRS}     PARENT_SCOPE)
      SET(${TRIBITS_PACKAGE}_INSTALLATION_LIBRARY_DIRS     ${${TRIBITS_PACKAGE}_LIBRARY_DIRS}     PARENT_SCOPE)
      SET(${TRIBITS_PACKAGE}_INSTALLATION_LIBRARIES        ${${TRIBITS_PACKAGE}_LIBRARIES}        PARENT_SCOPE)
      SET(${TRIBITS_PACKAGE}_INSTALLATION_TPL_INCLUDE_DIRS ${${TRIBITS_PACKAGE}_TPL_INCLUDE_DIRS} PARENT_SCOPE)
      SET(${TRIBITS_PACKAGE}_INSTALLATION_TPL_LIBRARY_DIRS ${${TRIBITS_PACKAGE}_TPL_LIBRARY_DIRS} PARENT_SCOPE)
      SET(${TRIBITS_PACKAGE}_INSTALLATION_TPL_LIBRARIES    ${${TRIBITS_PACKAGE}_TPL_LIBRARIES}    PARENT_SCOPE)
    ENDFOREACH()

  ENDIF()
ENDFUNCTION()
