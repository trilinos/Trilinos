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


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.generic.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE $ENV{JENKINS_COMM_TYPE})
SET(BUILD_TYPE $ENV{JENKINS_BUILD_TYPE})
SET(BUILD_DIR_NAME Sierra_${COMM_TYPE}_${BUILD_TYPE}_DEV_ETI_SERIAL-$ENV{JENKINS_DO_SERIAL}_OPENMP-$ENV{JENKINS_DO_OPENMP}_PTHREAD-$ENV{JENKINS_DO_PTHREAD}_CUDA-$ENV{JENKINS_DO_CUDA}_COMPLEX-$ENV{JENKINS_DO_COMPLEX})
SET(CTEST_PARALLEL_LEVEL 16)

# Note: CTEST_TEST_TYPE drives some side-effects in Tribits that should be 
#       taken into account.  If CTEST_TEST_TYPE is Experimental, Tribits will
#       override Trilinos_TRACK and *always* assign results to the Experimental
#       track on CDash.  Also, Tribits does different things based on CTEST_TEST_TYPE
#       being one of the 3 blessed types (Nightly, Continuous, Experimental), so it's
#       best to keep CTEST_TEST_TYPE one of these values.  
SET(CTEST_TEST_TYPE Nightly)
SET(Trilinos_TRACK  $ENV{JENKINS_JOB_TYPE})

SET(CTEST_TEST_TIMEOUT 900)
# SET(CTEST_DO_SUBMIT FALSE)

SET(Trilinos_PACKAGES Amesos Amesos2 Anasazi AztecOO Belos Epetra EpetraExt FEI Ifpack Ifpack2 Intrepid Kokkos ML MueLu NOX Pamgen RTOp Sacado Shards Teuchos Thyra Tpetra TrilinosSS Triutils Xpetra Zoltan Zoltan2)

SET(EXTRA_EXCLUDE_PACKAGES Galeri Intrepid2 Isorropia Stratimikos Teko SEACAS STK)

SET(EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF"
  "-DTPL_ENABLE_SuperLU=OFF"

  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTeuchos_ENABLE_COMPLEX:BOOL=$ENV{JENKINS_DO_COMPLEX}"
  "-DTrilinos_ENABLE_OpenMP:BOOL=$ENV{JENKINS_DO_OPENMP}"
  "-DTPL_ENABLE_HWLOC:STRING=OFF"
  "-DTPL_ENABLE_CUDA:STRING=$ENV{JENKINS_DO_CUDA}"
  
  "-DFEI_AZTECOO:BOOL=ON"
  "-DTrilinos_ENABLE_AztecOO:BOOL=ON"

  "-DTrilinos_ENABLE_EpetraExt:BOOL=ON"

  "-DKokkos_ENABLE_Pthread:BOOL=$ENV{JENKINS_DO_PTHREAD}"
  "-DKokkos_ENABLE_Cuda_UVM:BOOL=$ENV{JENKINS_DO_CUDA}"
  "-DKokkos_ENABLE_OpenMP:BOOL=$ENV{JENKINS_DO_OPENMP}"

  "-DMETIS_INCLUDE_DIRS=$ENV{SEMS_PARMETIS_INCLUDE_PATH}"
  "-DMETIS_LIBRARY_DIRS=$ENV{SEMS_PARMETIS_LIBRARY_PATH}"

  "-DNOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA:BOOL=OFF"
  "-DNOX_ENABLE_EPETRA:BOOL=ON"

  "-DTpetra_INST_SERIAL:BOOL=$ENV{JENKINS_DO_SERIAL}"
  "-DTpetra_INST_OPENMP:BOOL=$ENV{JENKINS_DO_OPENMP}"
  "-DTpetra_INST_PTHREAD:BOOL=$ENV{JENKINS_DO_PTHREAD}"
  "-DTpetra_INST_CUDA:BOOL=$ENV{JENKINS_DO_CUDA}"
  "-DTrilinos_ENABLE_COMPLEX_DOUBLE:BOOL=$ENV{JENKINS_DO_COMPLEX}"
  "-DTpetra_INST_COMPLEX_FLOAT:BOOL=OFF"
  "-DTpetra_INST_INT_INT:BOOL=ON"
  "-DTpetra_INST_INT_LONG:BOOL=ON"
  "-DTpetra_INST_INT_LONG_LONG:BOOL=OFF"
  "-DTpetra_INST_FLOAT:BOOL=OFF"
  "-DTpetra_INST_DOUBLE:BOOL=ON"
  "-DTpetraCore_ENABLE_TSQR:BOOL=ON"
  "-DTpetra_HIDE_DEPRECATED_CODE:BOOL=ON" 
  "-DTrilinos_ENABLE_TpetraTSQR:BOOL=ON" 
  "-DTrilinos_ENABLE_Tpetra:BOOL=ON" 
  -"DKokkos_ENABLE_Cuda_Lambda:BOOL=$ENV{JENKINS_DO_COMPLEX}"

  "-DMueLu_ENABLE_Epetra=OFF"
  "-DMueLu_ENABLE_Tpetra=ON"

  "-DAmesos2_ENABLE_Epetra=OFF"
  "-DAmesos2_CONFIGURE_OPTIONS_FILE:FILEPATH=cmake/ctest/drivers/parameterized/UMFPACK_Settings.cmake"

  "-DCMAKE_CXX_FLAGS:STRING=$ENV{JENKINS_CXX_FLAGS}"
  "-DCMAKE_C_FLAGS:STRING=$ENV{JENKINS_C_FLAGS}"
  "-DCMAKE_Fortran_FLAGS=$ENV{JENKINS_Fortran_FLAGS}"

  "-DTrilinos_C_Standard=gnu11"

)

IF (DEFINED ENV{JENKINS_BLAS_LIBRARY_DIRS})
    SET(EXTRA_CONFIGURE_OPTIONS ${EXTRA_CONFIGURE_OPTIONS}
        "-DBLAS_LIBRARY_DIRS=$ENV{JENKINS_BLAS_LIBRARY_DIRS}")
ENDIF()

IF (DEFINED ENV{JENKINS_LAPACK_LIBRARY_DIRS})
    SET(EXTRA_CONFIGURE_OPTIONS ${EXTRA_CONFIGURE_OPTIONS}
        "-DLAPACK_LIBRARY_DIRS=$ENV{JENKINS_LAPACK_LIBRARY_DIRS}")
ENDIF()

IF (DEFINED ENV{JENKINS_CMAKE_EXE_LINKER_FLAGS})
    SET(EXTRA_CONFIGURE_OPTIONS ${EXTRA_CONFIGURE_OPTIONS}
        "-DTrilinos_EXTRA_LINK_FLAGS=$ENV{JENKINS_CMAKE_EXE_LINKER_FLAGS}")
ENDIF()

IF (DEFINED ENV{JENKINS_CMAKE_DISABLE_MPI_WRAPPER})
    SET(EXTRA_CONFIGURE_OPTIONS ${EXTRA_CONFIGURE_OPTIONS}
        "-DMPI_USE_COMPILER_WRAPPERS=OFF")
ENDIF()

#"-DMPI_EXEC_POST_NUMPROCS_FLAGS:STRING=-bind-to;socket;--map-by;socket"
#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
