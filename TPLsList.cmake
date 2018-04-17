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

TRIBITS_REPOSITORY_DEFINE_TPLS(
  MKL             "cmake/TPLs/"    EX
  yaml-cpp        "cmake/TPLs/"    EX
  Peano           "cmake/TPLs/"    EX
  CUDA            "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"    ST
  CUSPARSE        "cmake/TPLs/"    ST
  Thrust          "cmake/TPLs/"    ST
  Cusp            "cmake/TPLs/"    ST
  TBB             "cmake/TPLs/"    EX
  Pthread         "cmake/TPLs/"    PT
  HWLOC           "cmake/TPLs/"    ST
  QTHREAD         "cmake/TPLs/"    ST
  BinUtils        "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    ST
  ARPREC          "packages/teuchos/cmake/tpls/"    ST
  QD              "packages/teuchos/cmake/tpls/"    ST
  MPI             "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/" PT
  BLAS            "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"   PT
  LAPACK          "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"   PT
  Boost           "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"   PT
  Scotch          "cmake/TPLs/"    ST
  OVIS            "cmake/TPLs/"    ST
  gpcd            "cmake/TPLs/"    ST
  METIS           "cmake/TPLs/"    TS
  MTMETIS         "cmake/TPLs/"    EX
  ParMETIS        "cmake/TPLs/"    PT
  PuLP            "cmake/TPLs/"    EX
  TopoManager     "cmake/TPLs/"    EX
  LibTopoMap      "cmake/TPLs/"    ST
  PaToH           "cmake/TPLs/"    ST
  CppUnit         "cmake/TPLs/"    ST
  ADOLC           "cmake/TPLs/"    ST
  ADIC            "cmake/TPLs/"    EX
  TVMET           "cmake/TPLs/"    ST
  MF              "cmake/TPLs/"    ST
  ExodusII        "cmake/TPLs/"    ST
  Nemesis         "cmake/TPLs/"    ST
  XDMF            "cmake/TPLs/"    TS
  Zlib            "cmake/TPLs/"    PT
  HDF5            "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"  PT
  CGNS            "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"  PT
  Pnetcdf         "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"  PT
  Netcdf          "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"  PT
  y12m            "cmake/TPLs/"    ST
  SuperLUDist     "cmake/TPLs/"    ST
  SuperLUMT	  "cmake/TPLs/"	   ST
  SuperLU         "cmake/TPLs/"    PT
  Cholmod	  "cmake/TPLs/"	   EX
  UMFPACK         "cmake/TPLs/"    ST
  MA28            "cmake/TPLs/"    TS
  AMD             "cmake/TPLs/"    TS
  CSparse         "cmake/TPLs/"    EX
  HYPRE           "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"   EX
  PETSC           "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"   ST
  BLACS           "cmake/TPLs/"    ST
  SCALAPACK       "cmake/TPLs/"    ST
  MUMPS           "cmake/TPLs/"    ST
  PARDISO_MKL     "cmake/TPLs/"    EX
  PARDISO         "cmake/TPLs/"    EX
  Oski            "cmake/TPLs/"    ST
  TAUCS           "cmake/TPLs/"    ST
  ForUQTK         "cmake/TPLs/"    EX
  Dakota          "cmake/TPLs/"    EX
  HIPS            "cmake/TPLs/"    EX
  MATLAB          "cmake/TPLs/"    EX
  CASK            "cmake/TPLs/"    EX
  SPARSKIT        "cmake/TPLs/"    ST
  QT              "packages/teuchos/cmake/tpls/"    ST
  gtest           "cmake/TPLs/"    EX
  BoostLib        "cmake/TPLs/"    PT
  BoostAlbLib     "cmake/TPLs/"    ST
  OpenNURBS       "cmake/TPLs/"    EX
  Portals         "cmake/TPLs/"    ST
  CrayPortals     "cmake/TPLs/"    ST
  Gemini          "cmake/TPLs/"    ST
  InfiniBand      "cmake/TPLs/"    ST
  BGPDCMF         "cmake/TPLs/"    ST
  BGQPAMI         "cmake/TPLs/"    ST
  Pablo           "cmake/TPLs/"    ST
  HPCToolkit      "cmake/TPLs/"    ST
  Clp             "cmake/TPLs/"    EX
  GLPK            "cmake/TPLs/"    EX
  qpOASES         "cmake/TPLs/"    EX
  Matio           "cmake/TPLs/"    ST
  PAPI            "cmake/TPLs/"    ST
  MATLABLib       "cmake/TPLs/"    EX
  Eigen           "packages/teuchos/cmake/tpls/"    EX
  X11             "cmake/TPLs/"    ST
  Lemon           "cmake/TPLs/"    EX
  GLM             "cmake/TPLs/"    EX
  quadmath        "cmake/TPLs/"    EX
  CAMAL           "cmake/TPLs/"    ST
  RTlib           "cmake/TPLs/"    ST
  DLlib           "cmake/TPLs/"    ST
  AmgX            "cmake/TPLs/"    EX
  CGAL            "cmake/TPLs/"    EX
  CGALCore        "cmake/TPLs/"    EX
  VTune           "cmake/TPLs/"    ST
  TASMANIAN       "cmake/TPLs/"    EX
  ArrayFireCPU    "cmake/TPLs/"    EX
  SimMesh         "SCOREC/cmake/TPLs/"    EX
  SimModel        "SCOREC/cmake/TPLs/"    EX
  SimParasolid    "SCOREC/cmake/TPLs/"    EX
  SimAcis         "SCOREC/cmake/TPLs/"    EX
  SimField        "SCOREC/cmake/TPLs/"    EX
  Valgrind        "cmake/TPLs/"    EX
  QUO             "cmake/TPLs/"    EX
  ViennaCL        "cmake/TPLs/"    EX
  Avatar          "cmake/TPLs/"    EX
  )

# NOTES:
#
# (*) ParMETIS must be listed after Scotch because the
#     ParMETIS include directories must come before the
#     Scotch include directories.
#
