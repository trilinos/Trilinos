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


# This list is just used for unit testing the dependency handling
# CMake code.  The reason that we have a separate list is so that we
# can keep very stable unit tests.


SET( Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  TrilinosFramework   cmake                           PS
  Teuchos             packages/teuchos                PS
  RTOp                packages/rtop                   PS
  Epetra              packages/epetra                 PS 
  Zoltan              packages/zoltan                 PS
  Shards              packages/shards                 PS
  Triutils            packages/triutils               PS
  Tpetra              packages/tpetra                 PS
  EpetraExt           packages/epetraext              PS
  Stokhos             packages/stokhos                EX
  Sacado              packages/sacado                 SS
  Thyra               packages/thyra                  PS
  Isorropia           packages/isorropia              PS
  AztecOO             packages/aztecoo                PS
  Galeri              packages/galeri                 PS
  Amesos              packages/amesos                 PS
  Intrepid            packages/intrepid               PS
  Ifpack              packages/ifpack                 PS
  ML                  packages/ml                     PS
  Belos               packages/belos                  SS
  Stratimikos         packages/stratimikos            PS
  RBGen               packages/rbgen                  PS
  Phalanx             packages/phalanx                SS
  Panzer              packages/panzer                 SS
  )

# NOTE: Sacado is really PS but for testing purpose it is made SS
# NOTE: Belos is really PS but for testing purpose it is made SS

PACKAGE_DISABLE_ON_PLATFORMS(ML BadSystem1)
PACKAGE_DISABLE_ON_PLATFORMS(Ifpack BadSystem1 BadSystem2)
