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
# Define the TriBITS package or not and decide what version of TriBITS to
# point to.
#
IF (TriBITS_SOURCE_DIR)
  # The TriBITS repo is defined, so let that repo define the TriBITS package.
  SET(TRIBITS_PACKAGE_LINE)
  # NOTE: This is the case when Trilinos is just a repo in a larger TriBITS
  # meta-project that will be pulling in TriBITS as its own repo and that repo
  # will define the TriBITS package and will ignore the version of TriBITS
  # under Trilinos/cmake/tribits/.
ELSEIF(${PROJECT_NAME}_TRIBITS_PACKAGE_USE_TRIBITS_DIR)
  # There is no TriBITS repo defined so let the Trilinos repo define the
  # TriBITS package and we want to use the version of TriBITS in
  # ${PROJECT_NAME}_TRIBITS_DIR.  WARNING: Only use this when
  # ${PROJECT_NAME}_TRIBITS_DIR is a subdir of ${PROJECT_SOURCE_DIR}!
  SET(TRIBITS_PACKAGE_LINE
    TriBITS   "${${PROJECT_NAME}_TRIBITS_DIR}"   PS
    )
  # NOTE: We use ${${PROJECT_NAME}_TRIBITS_DIR}, *not* hard-coded
  # cmake/tribits/ in case we are pointing to a different TriBITS
  # implementation when doing TriBITS/Trilinos co-development.
ELSE()
  # There is no TriBITS repo defined so let the Trilinos repo define the
  # TriBITS package and use the version of TriBITS in Trilinos for the TriBITS
  # package, not the one in ${PROJECT_NAME}_TRIBITS_DIR!
  SET(TRIBITS_PACKAGE_LINE
    TriBITS   "cmake/tribits"   PS
    )
  # NOTE: It is important *not* to use ${PROJECT_NAME}_TRIBITS_DIR when doing
  # automated CTest/CDash testing where ${PROJECT_NAME}_TRIBITS_DIR is defined
  # outside of PROJECT_SOURCE_DIR!  We want to be running the inner TriBITS
  # source build, not the outer TriBITS in this case.  This also maintains
  # prefect backward compatibility w.r.t. the definition of a TriBITS package
  # under Trilinos.
ENDIF()


#
# Define the Trilinos packages
#
TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  ${TRIBITS_PACKAGE_LINE}
  Teuchos               packages/teuchos                  PS
  ThreadPool            packages/ThreadPool               PS # Depends on ptheads system library
  RTOp                  packages/rtop                     PS
  Gtest                 commonTools/gtest                 SS
  Kokkos                packages/kokkos                   PS
  Sacado                packages/sacado                   PS
  Epetra                packages/epetra                   PS
  Zoltan                packages/zoltan                   PS
  Shards                packages/shards                   PS
  GlobiPack             packages/globipack                PS
  Triutils              packages/triutils                 PS
  Tpetra                packages/tpetra                   PS
  EpetraExt             packages/epetraext                PS
  Xpetra                packages/xpetra                   PS
  Thyra                 packages/thyra                    PS
  OptiPack              packages/optipack                 PS
  Isorropia             packages/isorropia                PS
  Pliris                packages/pliris                   PS
  Claps                 packages/claps                    EX
  AztecOO               packages/aztecoo                  PS
  Galeri                packages/galeri                   PS
  Amesos                packages/amesos                   PS
  Amesos2               packages/amesos2                  SS
  Pamgen                packages/pamgen                   PS
  SEACAS                packages/seacas                   SS # Depends on netcdf, optionally hdf5, xdmf, pamgen
  Trios                 packages/trios                    EX #temporary
  Ifpack                packages/ifpack                   PS
  Komplex               packages/komplex                  PS
  ML                    packages/ml                       PS
  Belos                 packages/belos                    PS
  Anasazi               packages/anasazi                  PS
  Zoltan2               packages/zoltan2                  SS
  Ifpack2               packages/ifpack2                  PS
  ShyLU                 packages/shylu                    EX
  Stratimikos           packages/stratimikos              PS
  FEI                   packages/fei                      PS
  Teko                  packages/teko                     SS
  RBGen                 packages/rbgen                    PS
  TriKota               packages/TriKota                  SS
  Intrepid              packages/intrepid                 PS
  STK                   packages/stk                      SS # Depends on boost
  Phalanx               packages/phalanx                  SS
  Phdmesh               packages/phdmesh                  EX # to be replaced by STK
  NOX                   packages/nox                      PS
  Moertel               packages/moertel                  PS
  MueLu                 packages/muelu                    SS
  Rythmos               packages/rythmos                  PS
  MOOCHO                packages/moocho                   PS
  Aristos               packages/aristos                  EX
  Stokhos               packages/stokhos                  SS
  Piro                  packages/piro                     SS
  Panzer                packages/panzer                   SS
  Sundance              packages/Sundance                 SS # Could be PS based on deps (BUG: 4669)
  CTrilinos             packages/CTrilinos                SS # Switched to SS to speed up checkin testing
  ForTrilinos           packages/ForTrilinos              EX
  PyTrilinos            packages/PyTrilinos               SS
  WebTrilinos           packages/WebTrilinos              EX # Should be SS
  Didasko               packages/didasko                  SS
  NewPackage            packages/new_package              EX # Should be SS
  Optika		packages/optika		          SS
  Mesquite              packages/mesquite                 PS
  MeshingGenie          packages/meshinggenie             EX
  TrilinosCouplings     packages/trilinoscouplings        SS
  FEApp                 demos/FEApp                       SS # Capability demonstration package
  )


#
# Disable certain packages on certain platforms.
#
# NOTE: This just makes the packages experimental 'EX' and therefore still
# allows the user to enable the package explicitly but the package will not
# get enabled implicitly.
#

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(MOOCHO Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Phalanx Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Phdmesh Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(PyTrilinos Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Sundance Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Tpetra Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Ifpack2 Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(TriKota Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Pamgen Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(STK Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(SEACAS Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Anasazi Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Zoltan Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Isorropia Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Teko Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Mesquite AIX)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Trios Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Xpetra Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Panzer Windows)
