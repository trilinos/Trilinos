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
# Platform/compiler specific options for geminga using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  IF(COMM_TYPE STREQUAL MPI)
    string(TOUPPER $ENV{SIERRA_MPI} UC_MPI_NAME)
    SET(BUILD_DIR_NAME ${UC_MPI_NAME}_${BUILD_TYPE}_${BUILD_NAME_DETAILS})
  ELSE()
    SET(BUILD_DIR_NAME ${COMM_TYPE}-${BUILD_TYPE}_${BUILD_NAME_DETAILS})
  ENDIF()

  SET(Trilinos_REPOSITORY_LOCATION_NIGHTLY_DEFAULT "git@gitlab-ex.sandia.gov:trilinos-project/Trilinos.git")
  SET(Trilinos_BRANCH "Tpetra_UVM_Removal" )

  SET(CTEST_DASHBOARD_ROOT  "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )
  SET(CTEST_NOTES_FILES     "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  SET(CTEST_BUILD_FLAGS     "-j45 -i" )

  SET_DEFAULT(CTEST_PARALLEL_LEVEL                  "45" )
  SET_DEFAULT(Trilinos_ENABLE_SECONDARY_TESTED_CODE ON)
  SET(Trilinos_CTEST_DO_ALL_AT_ONCE FALSE)
  SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES             ${EXTRA_EXCLUDE_PACKAGES} TriKota Optika Pamgen)

  # Select package disables
  set (Trilinos_ENABLE_Gtest OFF CACHE BOOL "Gtest just does not build" FORCE)
  set (Trilinos_ENABLE_ShyLU_NodeTacho OFF CACHE BOOL "Can't test Tacho with CUDA without RDC" FORCE)
  set (Trilinos_ENABLE_Shards OFF CACHE BOOL "Shards does not build" FORCE)
  set (Trilinos_ENABLE_Epetra OFF CACHE BOOL "We do not want Epetra" FORCE)

  # Select test disables
  set (KokkosCore_UnitTest_CudaTimingBased_MPI_1_DISABLE ON CACHE BOOL "Not to be run in nightly testing")


  SET(EXTRA_SYSTEM_CONFIGURE_OPTIONS
      "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"

      # Adding the following as a possible fix for github issue #2115.
      #KDD This flag appears to be unnecessary in April 2021, and it
      #KDD breaks building of Zoltan tests
      #KDD "-DCMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON"

      ### ALWAYS AND EVERYWHERE ###
      "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
      "-DBUILD_SHARED_LIBS:BOOL=ON"
      "-DTrilinos_ENABLE_TESTS:BOOL=ON"
      "-DTrilinos_ENABLE_EXAMPLES:BOOL=ON"
      "-DTrilinos_ENABLE_DEPENDENCY_UNIT_TESTS:BOOL=OFF"
      "-DTeuchos_GLOBALLY_REDUCE_UNITTEST_RESULTS:BOOL=ON"
      "-DTrilinos_ENABLE_COMPLEX=ON"
      "-DTeuchos_ENABLE_COMPLEX=ON"
      "-DTpetra_INST_COMPLEX_DOUBLE=ON"

      ### COMPILERS AND FLAGS ###
      "-DCMAKE_CXX_FLAGS:STRING='-Wall -Wno-unknown-pragmas -Wno-unused-but-set-variable -Wno-inline -Wshadow'"
      "-DTrilinos_ENABLE_Fortran:BOOL=OFF"

      ### TPLS ###
      "-DTPL_ENABLE_CUDA:BOOL=ON"
      "-DCMAKE_POLICY_DEFAULT_CMP0074=NEW"
      "-DTPL_ENABLE_CUSPARSE:BOOL=ON"
      "-DTPL_ENABLE_HWLOC:BOOL=OFF"


      # Host Blas is required (https://github.com/kokkos/kokkos-kernels/issues/347) for Kokkos-Kernels to build correctly
      "-DTPL_ENABLE_BLAS:BOOL=ON"
      "-DTPL_ENABLE_LAPACK:BOOL=ON"
      "-DTPL_BLAS_LIBRARIES=/usr/lib64/libblas.so"
      "-DTPL_LAPACK_LIBRARIES=/usr/lib64/liblapack.so"

      ### PACKAGE CONFIGURATION ###
      "-DKokkos_ENABLE_CUDA:BOOL=ON"
      "-DKokkos_ENABLE_CUDA_LAMBDA:BOOL=ON"
      "-DKokkos_ARCH_SKX:BOOL=ON"
      "-DKokkos_ARCH_VOLTA70:BOOL=ON"
      "-DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=OFF"

      "-DTrilinos_ENABLE_Epetra:BOOL=OFF"
      "-DTrilinos_ENABLE_Gtest:BOOL=OFF"
      "-DTrilinos_ENABLE_Pamgen:BOOL=OFF"
      "-DTrilinos_ENABLE_Shards:BOOL=OFF"
      "-DTrilinos_ENABLE_ShyLU_Node:BOOL=OFF"
      "-DTrilinos_ENABLE_ShyLU_NodeTacho:BOOL=OFF"            
      "-DTrilinos_ENABLE_ShyLU:BOOL=OFF"
      "-DTrilinos_ENABLE_ShyLU_DD:BOOL=OFF"
      "-DAmesos2_ENABLE_ShyLU_NodeTacho:BOOL=OFF"
      "-DAmesos2_ENABLE_ShyLU_NodeBasker:BOOL=OFF"

      ### MISC ###
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"
  )

  SET_DEFAULT(COMPILER_VERSION "$ENV{SIERRA_PLATFORM}")

  # Ensure that MPI is on for all parallel builds that might be run.
  IF(COMM_TYPE STREQUAL MPI)

    SET(EXTRA_SYSTEM_CONFIGURE_OPTIONS
        ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
        "-DTPL_ENABLE_MPI:BOOL=ON"
        "-DMPI_BASE_DIR:PATH=$ENV{MPIHOME}"
       )

  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()
