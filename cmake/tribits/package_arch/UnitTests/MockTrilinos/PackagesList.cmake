# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
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
