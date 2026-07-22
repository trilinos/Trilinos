// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Basic testing of Zoltan2::BasicIdentifierAdapter 

#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>

int main(int narg, char *arg[]) {
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0, gfail = 0;

  // Create global identifiers with weights
  zlno_t numLocalIds = 10;
  int nWeights = 2;

  zgno_t *myIds = new zgno_t[numLocalIds];
  zscalar_t *weights = new zscalar_t [numLocalIds*nWeights];
  zgno_t myFirstId = rank * numLocalIds * numLocalIds;

  for (zlno_t i=0; i < numLocalIds; i++){
    myIds[i] = zgno_t(myFirstId+i);
    weights[i*nWeights] = 1.0;
    weights[i*nWeights + 1] = (nprocs-rank) / (i+1);
  }

  typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> userTypes_t;
  std::vector<const zscalar_t *> weightValues;
  std::vector<int> strides;

  weightValues.push_back(weights);
  weightValues.push_back(weights + 1);
  strides.push_back(2);
  strides.push_back(2);

  Zoltan2::BasicIdentifierAdapter<userTypes_t> ia(numLocalIds, myIds,
    weightValues, strides);

  if (!fail && ia.getLocalNumIDs() != size_t(numLocalIds)){
    fail = 4;
  }
  if (!fail && ia.getNumWeightsPerID() != nWeights) {
    fail = 5;
  }

  const zgno_t *globalIdsIn;
  zscalar_t const *weightsIn[2];
  int weightStridesIn[2];

  ia.getIDsView(globalIdsIn);

  for (int w=0; !fail && w < nWeights; w++) {
    ia.getWeightsView(weightsIn[w], weightStridesIn[w], w);
  }

  const zscalar_t *w1 = weightsIn[0];
  const zscalar_t *w2 = weightsIn[1];
  int incr1 = weightStridesIn[0];
  int incr2 = weightStridesIn[1];

  for (zlno_t i=0; !fail && i < numLocalIds; i++){
    if (globalIdsIn[i] != zgno_t(myFirstId+i)) {
      fail = 8;
    }
    if (!fail && w1[i*incr1] != 1.0) {
      fail = 9;
    }
    if (!fail && w2[i*incr2] != weights[i*nWeights+1]) {
      fail = 10;
    }
  }

  delete [] myIds;
  delete [] weights;

  gfail = globalFail(*comm, fail);
  if (gfail) {
    printFailureCode(*comm, fail);   // will exit(1)
  }
  if (rank == 0) {
    std::cout << "PASS" << std::endl;
  }
}

