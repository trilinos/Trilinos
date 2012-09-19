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


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test2.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME CONTINUOUS_${COMM_TYPE}_OPT_DEV_SHARED)
#SET(CTEST_TEST_TIMEOUT 900)

#override the default number of processors to run on.
SET( CTEST_BUILD_FLAGS "-j12 -i" )
SET( CTEST_PARALLEL_LEVEL "12" )

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

SET(EXTRA_EXCLUDE_PACKAGES Claps)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_ENABLE_DEBUG:BOOL=ON"
  "-DBUILD_SHARED_LIBS:BOOL=ON"
  "-DMPI_BASE_DIR:PATH=/home/trilinos"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/boost-1.49.0"
  "-DNetcdf_LIBRARY_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/pnetcdf_4.2/lib"
  "-DNetcdf_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/pnetcdf_4.2/include"
  "-DHDF5_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/phdf5-1.8.6/include"
  "-DHDF5_LIBRARY_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/phdf5-1.8.6/lib"
  "-DTPL_ENABLE_Expat:BOOL=ON"
  "-DTPL_ENABLE_LAMMPS:BOOL=ON"
  "-DTPL_ENABLE_couple:BOOL=ON"
  "-DTPL_ENABLE_SPPARKS:BOOL=ON"
  "-DTPL_ENABLE_SuperLU:BOOL=ON"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

SET(CTEST_TEST_TYPE Continuous)

# Set the following variables to reasonable values so the Jenkins job
# provides useful feedback, rather then always returning an error
# because they aren't defined.

# The Jenkins checkout of the scripts repository. 
SET(CTEST_SOURCE_DIRECTORY "$ENV{WORKSPACE}/Trilinos")
# The location where the Jenkins script gets run.
SET(CTEST_BINARY_DIRECTORY "$ENV{WORKSPACE}/continuous-development")
# The CTest command being used by Jenkins.
SET(CTEST_COMMAND "/home/trilinos/bin/ctest")

function(VISIBLE_MESSAGE message)
  message("\n***")
  message("*** ${message}")
  message("***\n")
endfunction()

# Check to see if we need to start with an empty binary directory or enabling
# modified packages only. First, get the current time in seconds past the
# epoch.
execute_process(
  COMMAND date +%s
  OUTPUT_VARIABLE _current_time
)
set(_timestamp_file timestamp.txt)
if(NOT EXISTS ${_timestamp_file})
  VISIBLE_MESSAGE("No timestamp file exists, performing a clean build.")
  SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY ON)
  SET(CTEST_ENABLE_MODIFIED_PACKAGES_ONLY OFF)
  file(WRITE ${_timestamp_file} ${_current_time}) 
else()
  file(READ ${_timestamp_file} _last_time)
  math(EXPR _difference "${_current_time} - ${_last_time}") 
  if(${_difference} GREATER 57600) # 16 hours
    VISIBLE_MESSAGE("Timestamp is more than 16 hours old, performing a clean build.")
    SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY ON)
    SET(CTEST_ENABLE_MODIFIED_PACKAGES_ONLY OFF)
    file(WRITE ${_timestamp_file} ${_current_time}) 
  else()
    VISIBLE_MESSAGE("Timestamp is less than 16 hours old, performing an incremental build.")
    SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY OFF)
    SET(CTEST_ENABLE_MODIFIED_PACKAGES_ONLY ON)
  endif()
endif()

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

