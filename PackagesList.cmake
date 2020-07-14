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
# Define the Trilinos packages
#
TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  TrilinosFrameworkTests  commonTools/framework           PT
  TrilinosATDMConfigTests cmake/std/atdm                  PT
  Gtest                 commonTools/gtest                 PT
  Kokkos                packages/kokkos                   PT
  Teuchos               packages/teuchos                  PT
  KokkosKernels         packages/kokkos-kernels           PT
  RTOp                  packages/rtop                     PT
  Sacado                packages/sacado                   PT
  MiniTensor            packages/minitensor               PT
  Epetra                packages/epetra                   PT
  SCOREClion            SCOREC/lion                       ST
  SCORECpcu             SCOREC/pcu                        ST
  SCORECgmi             SCOREC/gmi                        ST
  SCORECgmi_sim         SCOREC/gmi_sim                    ST
  SCORECapf             SCOREC/apf                        ST
  SCORECapf_sim         SCOREC/apf_sim                    ST
  SCORECmds             SCOREC/mds                        ST
  SCORECparma           SCOREC/parma                      ST
  SCORECspr             SCOREC/spr                        ST
  AvatarT               packages/avatart                  EX
  Zoltan                packages/zoltan                   PT
  Shards                packages/shards                   PT
  Triutils              packages/triutils                 PT
  EpetraExt             packages/epetraext                PT	
  Tpetra                packages/tpetra                   PT
  TrilinosSS            packages/common/auxiliarySoftware/SuiteSparse PT # Auxiliary software.
  Domi                  packages/domi                     PT
  Thyra                 packages/thyra                    PT
  Xpetra                packages/xpetra                   PT
  Isorropia             packages/isorropia                PT
  Pliris                packages/pliris                   ST
  AztecOO               packages/aztecoo                  PT
  Galeri                packages/galeri                   PT
  Amesos                packages/amesos                   PT
  Pamgen                packages/pamgen                   PT
  Zoltan2Core           packages/zoltan2/core             PT
  Ifpack                packages/ifpack                   PT
  ML                    packages/ml                       PT
  Belos                 packages/belos                    PT
  ShyLU_Node            packages/shylu/shylu_node         PT
  Amesos2               packages/amesos2                  PT
  SEACAS                packages/seacas                   PT # Depends on netcdf, optionally hdf5, xdmf, pamgen
  Komplex               packages/komplex                  ST
  Anasazi               packages/anasazi                  PT
  Ifpack2               packages/ifpack2                  PT
  Stratimikos           packages/stratimikos              PT
  FEI                   packages/fei                      PT
  Teko                  packages/teko                     PT
  TriKota               packages/TriKota                  ST
  Intrepid              packages/intrepid                 PT
  Intrepid2             packages/intrepid2                PT
  Compadre              packages/compadre                 ST
  STK                   packages/stk                      PT # Depends on boost
  Percept               packages/percept                  PT # Depends on boost
  SCORECapf_zoltan      SCOREC/zoltan                     ST
  SCORECapf_stk         SCOREC/stk                        ST
  SCORECma              SCOREC/ma                         ST
  SCORECpumi            SCOREC/pumi                       ST
  SCOREC                SCOREC                            ST
  Phalanx               packages/phalanx                  PT
  NOX                   packages/nox                      PT
  Moertel               packages/moertel                  ST
  MueLu                 packages/muelu                    PT
  Zoltan2Sphynx         packages/zoltan2/sphynx           PT
  Zoltan2               packages/zoltan2                  PT
  ShyLU_DD              packages/shylu/shylu_dd           PT
  ShyLU                 packages/shylu                    PT
  Rythmos               packages/rythmos                  PT
  Tempus                packages/tempus                   PT
  MOOCHO                packages/moocho                   ST
  Stokhos               packages/stokhos                  PT
  ROL                   packages/rol                      PT
  Piro                  packages/piro                     PT
  Panzer                packages/panzer                   PT
  CTrilinos             packages/CTrilinos                ST # Switched to ST to speed up checkin testing
#  ForTrilinos           packages/ForTrilinos              EX
  PyTrilinos            packages/PyTrilinos               ST
  WebTrilinos           packages/WebTrilinos              EX # Should be ST
  NewPackage            packages/new_package              EX # Should be ST
  Optika		packages/optika		          EX
  Adelus                packages/adelus                   PT
  TrilinosCouplings     packages/trilinoscouplings        PT
  Pike                  packages/pike                     PT
  xSDKTrilinos          packages/xSDKTrilinos             ST
  )

# Allow builds even if some packages are missing
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(AvatarT)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCOREC)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCOREClion)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECgmi)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECgmi_sim)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECpcu)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf_sim)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECmds)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECparma)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECspr)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf_stk)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf_zoltan)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECma)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECpumi)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Avatar)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(MOOCHO)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Sundance)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(CTrilinos)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(ForTrilinos)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Optika)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Mesquite)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(WebTrilinos)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(xSDKTrilinos)

#
# Disable certain packages on certain platforms.
#
# NOTE: This just makes the packages experimental 'EX' and therefore still
# allows the user to enable the package explicitly but the package will not
# get enabled implicitly.
#

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(MOOCHO Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Phalanx Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(PyTrilinos Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Sundance Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Tpetra Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Ifpack2 Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(TriKota Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Pamgen Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(STK Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(SEACAS Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Anasazi Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Isorropia Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Zoltan Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Teko Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Panzer Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Compadre Windows)
