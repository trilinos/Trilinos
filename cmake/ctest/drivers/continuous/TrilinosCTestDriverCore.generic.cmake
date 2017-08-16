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

SET(PythonInterp_FIND_VERSION 2.7)
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for typhon using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME

  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET_DEFAULT( CTEST_PARALLEL_LEVEL "16" )
  SET_DEFAULT( CTEST_BUILD_FLAGS "-j${CTEST_PARALLEL_LEVEL} -i" )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF)

  SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES ${EXTRA_EXCLUDE_PACKAGES} PyTrilinos TriKota )

  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DTrilinos_ENABLE_COMPLEX:BOOL=ON"
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"
    "-DBUILD_SHARED_LIBS:BOOL=ON"
    "-DIntrepid2_ENABLE_DEBUG_INF_CHECK:BOOL=OFF"
    "-DPythonInterp_FIND_VERSION=${PythonInterp_FIND_VERSION}"

    "-DTPL_ENABLE_BLAS:BOOL=ON"
    "-DTPL_ENABLE_LAPACK:BOOL=ON"

    "-DParMETIS_INCLUDE_DIRS:PATH=$ENV{SEMS_PARMETIS_INCLUDE_PATH}"
    "-DParMETIS_LIBRARY_DIRS:PATH=$ENV{SEMS_PARMETIS_LIBRARY_PATH}"

    "-DSuperLU_INCLUDE_DIRS:PATH=$ENV{SEMS_SUPERLU_ROOT}/include"
    "-DSuperLU_LIBRARY_DIRS:PATH=$ENV{SEMS_SUPERLU_ROOT}/lib"
    "-DSuperLU_LIBRARY_NAMES:STRING=superlu"

    "-DBoost_INCLUDE_DIRS:FILEPATH=$ENV{SEMS_BOOST_ROOT}/include"
    "-DBoost_LIBRARY_DIRS:FILEPATH=$ENV{SEMS_BOOST_ROOT}/lib"

    "-DBoostLib_INCLUDE_DIRS:FILEPATH=$ENV{SEMS_BOOST_ROOT}/include"
    "-DBoostLib_LIBRARY_DIRS:FILEPATH=$ENV{SEMS_BOOST_ROOT}/lib"

    "-DHDF5_INCLUDE_DIRS:FILEPATH=$ENV{SEMS_HDF5_ROOT}/include"
    "-DHDF5_LIBRARY_DIRS:FILEPATH=$ENV{SEMS_HDF5_ROOT}/lib"

    "-DNetcdf_INCLUDE_DIRS:FILEPATH=$ENV{SEMS_NETCDF_ROOT}/include"
    "-DNetcdf_LIBRARY_DIRS:FILEPATH=$ENV{SEMS_NETCDF_ROOT}/lib"

    "-DZlib_INCLUDE_DIRS:FILEPATH=$ENV{SEMS_ZLIB_ROOT}/include"
    "-DZlib_LIBRARY_DIRS:FILEPATH=$ENV{SEMS_ZLIB_ROOT}/lib"
 
    )

  IF (BUILD_TYPE STREQUAL "DEBUG")
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTrilinos_ENABLE_DEBUG:BOOL=ON"
      )
  ENDIF()

  SET_DEFAULT(COMPILER_VERSION "$ENV{SEMS_COMPILER_NAME}-$ENV{SEMS_COMPILER_VERSION}")

  #Ensuring that MPI is on for all parallel builds that might be run.
  IF(COMM_TYPE STREQUAL MPI)
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
        ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
        "-DTPL_ENABLE_MPI:BOOL=ON"
        "-DMPI_BASE_DIR=$ENV{SEMS_OPENMPI_ROOT}"
       )
  ELSE()
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
        ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
        "-DCMAKE_CXX_COMPILER:FILEPATH=$ENV{CXX}"
        "-DCMAKE_C_COMPILER:FILEPATH=$ENV{CC}"
        "-DCMAKE_Fortran_COMPILER:FILEPATH=$ENV{FC}"
       )
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()
