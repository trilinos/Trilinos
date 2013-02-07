// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Basic testing of Zoltan2::BasicIdentifierAdapter 

#include <Zoltan2_BasicIdentifierInput.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0, gfail=0;

  // Create global identifiers with weights

  lno_t numLocalIds = 10;
  int weightDim = 2;

  gno_t *myIds = new gno_t [numLocalIds];
  scalar_t *weights = new scalar_t [numLocalIds*weightDim];
  gno_t base = rank * numLocalIds * numLocalIds;

  for (lno_t i=0; i < numLocalIds; i++){
    myIds[i] = base+i;
    weights[i*weightDim] = 1.0;
    weights[i*weightDim + 1] = (nprocs-rank) / (i+1);
  }

  // Create a Zoltan2::BasicIdentifierAdapter object
  // and verify that it is correct

  typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> userTypes_t;
  std::vector<const scalar_t *> weightValues;
  std::vector<int> strides;

  weightValues.push_back(weights);
  weightValues.push_back(weights + 1);
  strides.push_back(2);
  strides.push_back(2);

  Zoltan2::BasicIdentifierAdapter<userTypes_t> ia( numLocalIds, myIds,
    weightValues, strides);

  if (!fail && ia.getLocalNumberOfIdentifiers() != size_t(numLocalIds)){
    fail = 4;
  }

  if (!fail && ia.getNumberOfWeights() != weightDim)
    fail = 5;

  const gno_t *globalIdsIn;
  scalar_t const *weightsIn[2];
  int weightStridesIn[2];

  if (!fail && ia.getIdentifierList(globalIdsIn) != size_t(numLocalIds))
    fail = 6;

  for (int w=0; !fail && w < weightDim; w++){
    if (ia.getIdentifierWeights(w, weightsIn[w], weightStridesIn[w]) <
        size_t(numLocalIds * weightStridesIn[w]))
      fail = 20;
  }

  const scalar_t *w1 = weightsIn[0];
  const scalar_t *w2 = weightsIn[1];
  int incr1 = weightStridesIn[0];
  int incr2 = weightStridesIn[1];

  for (lno_t i=0; !fail && i < numLocalIds; i++){

    if (globalIdsIn[i] != base+i)
      fail = 8;
    
    if (!fail && w1[i*incr1] != 1.0)
      fail = 9;
    
    if (!fail && w2[i*incr2] != weights[i*weightDim+1])
      fail = 10;
  }

  delete [] myIds;
  delete [] weights;

  gfail = globalFail(comm, fail);
  if (gfail)
    printFailureCode(comm, fail);   // will exit(1)

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

