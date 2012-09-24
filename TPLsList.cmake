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
# Define the list of TPLs, their find module names, and their classification
#
# TPL_NAME:
#
#   The name of the TPL used in the CMake cache variables TPL_ENABLE_${TPL_NAME}
#
# TPL_FINDMOD:
#
#   The name of the find module under that is used to get the names of the
#   TPLs.  If ends in '/' then this gives the directory and the standard module
#   name will be used which is FindTPL${TPL_NAME}.cmake.
#
# TPL_CLASSIFICATION:
#
#   PS: Primary Stable TPL
#
#     Primary Stable TPLs are those TPLs that a Trilinos developer must have
#     installed on their machine in order to be able to do Trilinos
#     development.  For example, we require that you have BLAS, LAPACK, and
#     MPI installed in order to do Trilinos development.  These are
#     fundamental dependencies that are needed in order to do precheckin
#     testing.
#
#   SS: Secondary Stable TPL
#
#     Secondary Stable TPLs are those TPLs that are not required in order to
#     be able to develop and test Trilinos before checkins but are none the
#     less offically supported.  Support for SS TPLs is tested as part of the
#     nightly testing process.
#
#   TS: Tertiary Stable TPL
#
#     Tertiary Stable TPLs are those TPLs that are supported TPLs but can not
#     be included in the set of SS TPLs because they may conflicit with other
#     SS Code.  For example, METIS is listed as a TS TPL because it conflicts
#     with ParMETIS which is declared as a SS TPL.
#
#   EX: Experimental TPL
#
#     Experimental TPLs are not offically supported.  They represent
#     experimental capabilities of Trilinos packages.  Support for EX TPLs is
#     never tested as part of the main nightly testing process.  However,
#     package developers are encouraged to set up their own nightly testing
#     for their EX TPLs for their packages.
#
# The default enable for all TPLs is empty "" reguardless of the category.
# The idea is that the enabling of the TPL will be done by the package and
# other enables that the user has to set.
#
# NOTE: The TPLs must be listed in the order of increasing dependencies (if
# such dependencies exist).
#

SET( Trilinos_TPLS_FINDMODS_CLASSIFICATIONS
  MKL             "cmake/TPLs/"    EX
  yaml-cpp        "cmake/TPLs/"    EX
  Peano           "cmake/TPLs/"    EX
  CUDA            "cmake/TPLs/"    SS
  CUSPARSE        "cmake/TPLs/"    SS
  Thrust          "cmake/TPLs/"    SS
  Cusp            "cmake/TPLs/"    SS
  TBB             "cmake/TPLs/"    EX
  Pthread         "cmake/TPLs/"    SS
  HWLOC           "cmake/TPLs/"    SS
  BinUtils        "cmake/TPLs/"    SS
  ARPREC          "cmake/TPLs/"    SS
  QD              "cmake/TPLs/"    SS
  MPI             "${${PROJECT_NAME}_TRIBITS_DIR}/tpls/"    PS
  BLAS            "cmake/TPLs/"    PS
  LAPACK          "cmake/TPLs/"    PS
  Boost           "cmake/TPLs/"    SS
  Scotch          "cmake/TPLs/"    SS
  OVIS            "cmake/TPLs/"    SS
  METIS           "cmake/TPLs/"    TS
  ParMETIS        "cmake/TPLs/"    SS
  PaToH           "cmake/TPLs/"    SS
  CppUnit         "cmake/TPLs/"    SS
  ADOLC           "cmake/TPLs/"    SS
  ADIC            "cmake/TPLs/"    EX
  TVMET           "cmake/TPLs/"    SS
  MF              "cmake/TPLs/"    SS
  ExodusII        "cmake/TPLs/"    SS
  Nemesis         "cmake/TPLs/"    SS
  XDMF            "cmake/TPLs/"    TS
  Netcdf          "cmake/TPLs/"    SS
  y12m            "cmake/TPLs/"    SS
  SuperLUDist     "cmake/TPLs/"    SS
  SuperLUMT	  "cmake/TPLs/"	   SS
  SuperLU         "cmake/TPLs/"    SS
  Zlib            "cmake/TPLs/"    SS
  UMFPACK         "cmake/TPLs/"    SS
  MA28            "cmake/TPLs/"    TS
  AMD             "cmake/TPLs/"    TS
  CSparse         "cmake/TPLs/"    EX
  PETSC           "cmake/TPLs/"    SS
  HYPRE           "cmake/TPLs/"    EX
  BLACS           "cmake/TPLs/"    SS
  SCALAPACK       "cmake/TPLs/"    SS
  MUMPS           "cmake/TPLs/"    SS
  PARDISO_MKL     "cmake/TPLs/"    EX
  Oski            "cmake/TPLs/"    SS
  TAUCS           "cmake/TPLs/"    SS
  ForUQTK         "cmake/TPLs/"    EX
  Dakota          "cmake/TPLs/"    EX
  HIPS            "cmake/TPLs/"    EX
  HDF5            "cmake/TPLs/"    EX
  MATLAB          "cmake/TPLs/"    EX
  CASK            "cmake/TPLs/"    EX
  SPARSKIT        "cmake/TPLs/"    SS
  QT              "cmake/TPLs/"    SS
  gtest           "cmake/TPLs/"    EX
  BoostLib        "cmake/TPLs/"    SS
  OpenNURBS       "cmake/TPLs/"    EX
  Portals         "cmake/TPLs/"    SS
  CrayPortals     "cmake/TPLs/"    SS
  Gemini          "cmake/TPLs/"    SS
  InfiniBand      "cmake/TPLs/"    SS
  Pablo           "cmake/TPLs/"    SS
  HPCToolkit      "cmake/TPLs/"    SS
  Pnetcdf         "cmake/TPLs/"    SS
  Clp             "cmake/TPLs/"    EX
  GLPK            "cmake/TPLs/"    EX
  qpOASES         "cmake/TPLs/"    EX
  Matio           "cmake/TPLs/"    SS
  PAPI            "cmake/TPLs/"    SS
  )

# NOTES:
#
# (*) ParMETIS must be listed after Scotch because the
#     ParMETIS include directories must come before the
#     Scotch include directories.
#
