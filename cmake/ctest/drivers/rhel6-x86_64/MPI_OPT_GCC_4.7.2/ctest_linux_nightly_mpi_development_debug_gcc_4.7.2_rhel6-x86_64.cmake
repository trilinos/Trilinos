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

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../TrilinosCTestDriverCore.rhel6-x86_64.gcc4.7.2-sems.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME openmpi-1.6.5_Debug_DEV_Werror)
#SET(CTEST_TEST_TIMEOUT 900)

# Note: To send output to the "Specialized" track on CDash, we need to
#       set Trilinos_TRACK to "Specialized" and set CTEST_TEST_TYPE to
#       something other than "Experimental".
#       Since "Nightly" "Continuous" and "Experimental" are special to
#       Tribits, we'll set CTEST_TEST_TYPE to "Nightly" here.
SET(CTEST_TEST_TYPE Nightly)
SET(Trilinos_TRACK  Specialized)

SET(Trilinos_ENABLE_SECONDARY_TESTED_CODE ON)

SET(EXTRA_EXCLUDE_PACKAGES MOOCHO Optika PyTrilinos Didasko)

#
# SET(Trilinos_PACKAGES NOX      EpetraExt   TrilinosCouplings
#                       Ifpack   Isorropia   AztecOO
#                       Belos    Teuchos     Amesos
#                       Sacado   Zoltan      Epetra
#                       Triutils )

SET(EXTRA_CONFIGURE_OPTIONS
  ${EXTRA_CONFIGURE_OPTIONS}
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"

  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DTPL_ENABLE_HDF5:BOOL=ON"
  "-DTPL_ENABLE_Netcdf:BOOL=ON"
  "-DTPL_ENABLE_SuperLU:BOOL=ON"

  "-DTrilinos_SHOW_DEPRECATED_WARNINGS=OFF"

  "-DTrilinos_ENABLE_MOOCHO=OFF"
  "-DTrilinos_ENABLE_Optika=OFF"
  "-DTrilinos_ENABLE_Didasko=OFF"

  "-DZoltan2_ENABLE_Experimental:BOOL=ON"

  "-DCMAKE_CXX_FLAGS:STRING=-Wall -ansi -pedantic -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
