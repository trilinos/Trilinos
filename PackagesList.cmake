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


INCLUDE(TribitsListHelpers)

#
# Define the Trilinos package names, directories, and classification.
#
# Package classifications are:
#
#   PS: Primary Stable Package
#
#     Primary Stable Packages have at least some Primary Stable Code which is
#     expected to be fully tested before every push to the global repo.  The
#     default enable for PS packages is empty "" which allows the PS package
#     to be enabled implicitly based on other criteria.  The option
#     Trilinos_ENABLE_ALL_PACKAGES=ON will cause all PS packages to be enabled
#     unless they are explicitly disabled.
#
#   SS: Secondary Stable Package
#
#     Secondary Stable Packages have no PS code or they would be classified as
#     PS packages.  A package must be classified as SS if it has a required
#     dependency on another SS package or SS TPL.  A package may also be
#     declared SS to avoid requiring it to be tested before every push to the
#     global repo.  For example, a package that does not provide any
#     significant functionally like Didasko is classified as a SS package even
#     through it could be classified as PS just based on its required package
#     and TPL dependencies.  SS packages will have their default enables set
#     to empty "".  This allows them to be enabled implicilty.  When
#     Trilinos_ENABLE_ALL_PACKAGES=ON but
#     Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF, the SS packages will not be
#     enabled.  However, when Trilinos_ENABLE_ALL_PACKAGES=ON and
#     Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON, then SS packages will be
#     enabled if they are not explicitly disabled.  Packages that are SS but
#     not PS must be disabled in pre-push testing.  However, SS packages are
#     tested by the post-push CI and nightly testing processes.
#
#   EX: Experimental Package
#
#     Experimental packages are those packages that contain no PS or SS
#     code. The default enable for EX packages is always OFF which requires
#     that they be explicitly enabled in order to be turned on. EX packages
#     must be disabled in pre-push testring and are not tested as part of the
#     post-push CI or nightly testing processes.  However, package developers
#     of EX pacakges are encouraged to set up their own nightly testing for
#     thier EX packages.
#
# NOTE: These packages must be listed in strictly assending order in terms of
# package dependencies.  If you get the order wrong, then an error message
# will be printed during configuration with CMake.
#

SET( Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  TrilinosFramework     cmake                             PS # Only tests, no libraries/capabilities!
  Teuchos               packages/teuchos                  PS
  ThreadPool            packages/ThreadPool               PS # Depends on ptheads system library
  Sacado                packages/sacado                   PS
  RTOp                  packages/rtop                     PS
  Kokkos                packages/kokkos                   PS
  Epetra                packages/epetra                   PS
  Zoltan                packages/zoltan                   PS
  Shards                packages/shards                   PS
  GlobiPack             packages/globipack                PS
  Triutils              packages/triutils                 PS
  Tpetra                packages/tpetra                   PS
  EpetraExt             packages/epetraext                PS
  Thyra                 packages/thyra                    PS
  OptiPack              packages/optipack                 PS
  Isorropia             packages/isorropia                PS
  Pliris                packages/pliris                   PS
  Claps                 packages/claps                    EX
  AztecOO               packages/aztecoo                  PS
  Galeri                packages/galeri                   PS
  Amesos2               packages/amesos2                  EX
  Amesos                packages/amesos                   PS
  Pamgen                packages/pamgen                   PS
  SEACAS                packages/seacas                   SS # Depends on netcdf, optionally hdf5, xdmf, pamgen
  Trios                 packages/trios                    EX #temporary
  Ifpack                packages/ifpack                   PS
  Komplex               packages/komplex                  PS
  ML                    packages/ml                       PS
  Belos                 packages/belos                    PS
  Ifpack2               packages/ifpack2                  PS
  ShyLU                 packages/shylu                    EX
  Stratimikos           packages/stratimikos              PS
  FEI                   packages/fei                      PS
  Anasazi               packages/anasazi                  PS
  Teko                  packages/teko                     SS
  RBGen                 packages/rbgen                    PS
  TriKota               packages/TriKota                  SS
  Intrepid              packages/intrepid                 PS
  STK                   packages/stk                      SS # Depends on boost
  Phalanx               packages/phalanx                  SS
  Phdmesh               packages/phdmesh                  EX # to be replaced by STK
  NOX                   packages/nox                      PS
  Moertel               packages/moertel                  PS
  TrilinosCouplings     packages/trilinoscouplings        SS
  Rythmos               packages/rythmos                  PS
  MOOCHO                packages/moocho                   PS
  Aristos               packages/aristos                  EX
  Xpetra                packages/xpetra                   SS
  MueLu                 packages/muelu                    SS
  Stokhos               packages/stokhos                  SS
  Piro                  packages/piro                     SS
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
  FEApp                 demos/FEApp                       SS # Capability demonstration package
  )


#
# Disable certain packages on certain platforms.
#
# NOTE: This just makes the packages experimental 'EX' and therefore still
# allows the user to enable the package explicitly but the package will not
# get enabled implicitly.
#

PACKAGE_DISABLE_ON_PLATFORMS(MOOCHO Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Phalanx Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Phdmesh Windows)
PACKAGE_DISABLE_ON_PLATFORMS(PyTrilinos Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Sundance Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Tpetra Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Ifpack2 Windows)
PACKAGE_DISABLE_ON_PLATFORMS(TriKota Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Pamgen Windows)
PACKAGE_DISABLE_ON_PLATFORMS(STK Windows)
PACKAGE_DISABLE_ON_PLATFORMS(SEACAS Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Anasazi Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Zoltan Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Isorropia Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Teko Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Mesquite AIX)
PACKAGE_DISABLE_ON_PLATFORMS(Trios Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Xpetra Windows)
