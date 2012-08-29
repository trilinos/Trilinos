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
// Test the following:
//         PartitioningSolutionQuality class
//         MetricValues class
//         Metric related namespace methods


#include <Zoltan2_PartitioningSolutionQuality.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicIdentifierInput.hpp>
#include <stdlib.h>
#include <vector>

typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> user_t;
typedef Zoltan2::BasicIdentifierInput<user_t> idInput_t;
typedef Zoltan2::PartitioningSolutionQuality<idInput_t> quality_t;

using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::arcp;

using namespace std;
using std::endl;
using std::cout;

typedef zoltan2_partId_t partId_t;

void doTest(RCP<const Comm<int> > comm, int numLocalObj,
  int weightDim, int numLocalParts, bool givePartSizes);

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();

  doTest(comm, 10, 0, -1, false);
  doTest(comm, 10, 0, 1, false);
  doTest(comm, 10, 0, 1, true);
  doTest(comm, 10, 1, 1, false);
  doTest(comm, 10, 1, 1, true);
  doTest(comm, 10, 2, 1, false);
  doTest(comm, 10, 2, 1, true);
  doTest(comm, 10, 1, 2, true);
  doTest(comm, 10, 1, 2, false);
  doTest(comm, 10, 1, -1, false);
  doTest(comm, 10, 1, -1, true);
  doTest(comm, 10, 2, -1, false);

  if (rank==0)
    cout << "PASS" << endl;
}

// Assumes numLocalObj is the same on every process.

#include <Zoltan2_GetParameter.hpp>


void doTest(RCP<const Comm<int> > comm, int numLocalObj,
  int weightDim, int numLocalParts, bool givePartSizes)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail=0;
  srand(rank+1);
  bool testEmptyParts = (numLocalParts < 1);
  int numGlobalParts = 0;

  if (testEmptyParts){
    numGlobalParts = nprocs / 2;
    if (numGlobalParts >= 1)
      numLocalParts = (rank < numGlobalParts ? 1 : 0);
    else{
      numLocalParts = 1;
      testEmptyParts = false;
    }
  }
  else{
    numGlobalParts = nprocs * numLocalParts;
  }

  if (rank == 0){
    cout << endl;
    cout << "Test: weight dimension " << weightDim;
    cout << ", desired number of parts " << numGlobalParts;
    if (givePartSizes)
      cout << ", with differing part sizes." << endl;
    else
      cout << ", with uniform part sizes." << endl;
    cout << "Number of procs " << nprocs;
    cout << ", each with " << numLocalObj << " objects, part = rank." << endl;
  }

  // An environment.  This is usually created by the problem.

  Teuchos::ParameterList pl("test list");
  pl.set("num_local_parts", double(numLocalParts));
  
  RCP<const Zoltan2::Environment> env = 
    rcp(new Zoltan2::Environment(pl, comm));

  // A simple identifier map.  Usually created by the model.

  gno_t *myGids = new gno_t [numLocalObj];
  for (int i=0, x=rank*numLocalObj; i < numLocalObj; i++, x++){
    myGids[i] = x;
  }

  ArrayRCP<const gno_t> gidArray(myGids, 0, numLocalObj, true);

  RCP<const Zoltan2::IdentifierMap<user_t> > idMap = 
    rcp(new Zoltan2::IdentifierMap<user_t>(env, comm, gidArray)); 

  // Part sizes.  Usually supplied by the user to the Problem.
  // Then the problem supplies them to the Solution.

  int partSizeDim = (givePartSizes ? (weightDim ? weightDim : 1) : 0);
  ArrayRCP<ArrayRCP<partId_t> > ids(partSizeDim);
  ArrayRCP<ArrayRCP<scalar_t> > sizes(partSizeDim);

  if (givePartSizes && numLocalParts > 0){
    partId_t *myParts = new partId_t [numLocalParts];
    myParts[0] = rank * numLocalParts;
    for (int i=1; i < numLocalParts; i++)
      myParts[i] = myParts[i-1] + 1;
    ArrayRCP<partId_t> partNums(myParts, 0, numLocalParts, true);

    scalar_t sizeFactor = nprocs/2 - rank;
    if (sizeFactor < 0) sizeFactor *= -1;
    sizeFactor += 1;

    for (int dim=0; dim < partSizeDim; dim++){
      scalar_t *psizes = new scalar_t [numLocalParts];
      for (int i=0; i < numLocalParts; i++)
        psizes[i] = sizeFactor;
      sizes[dim] = arcp(psizes, 0, numLocalParts, true);
      
      ids[dim] = partNums;
    }
  }

  // An input adapter with random weights.  Created by the user.

  std::vector<const scalar_t *> weights;
  std::vector<int> strides;   // default to 1

  int len = numLocalObj*weightDim;
  ArrayRCP<scalar_t> wgtBuf;
  scalar_t *wgts = NULL;

  if (len > 0){
    wgts = new scalar_t [len];
    wgtBuf = arcp(wgts, 0, len, true);
    for (int i=0; i < len; i++)
      wgts[i] = (scalar_t(rand()) / scalar_t(RAND_MAX)) + 1.0;
  }

  for (int i=0; i < weightDim; i++, wgts+=numLocalObj)
    weights.push_back(wgts);

  RCP<const idInput_t> ia;

  try{
    ia = rcp(new idInput_t(numLocalObj, myGids, weights, strides));
  }
  catch (std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "create adapter", 1);

  // A solution (usually created by a problem)

  RCP<Zoltan2::PartitioningSolution<idInput_t> > solution;

  try{
    if (givePartSizes)
      solution = rcp(new Zoltan2::PartitioningSolution<idInput_t>(
        env, comm, idMap, weightDim,
        ids.view(0,partSizeDim), sizes.view(0,partSizeDim)));
    else
      solution = rcp(new Zoltan2::PartitioningSolution<idInput_t>(
        env, comm, idMap, weightDim));
  }
  catch (std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "create solution", 1);

  // Part assignment for my objects: The algorithm usually calls this. 

  partId_t *partNum = new partId_t [numLocalObj];
  ArrayRCP<partId_t> partAssignment(partNum, 0, numLocalObj, true);
  for (int i=0; i < numLocalObj; i++)
    partNum[i] = rank;

  solution->setParts(gidArray, partAssignment);
  RCP<const Zoltan2::PartitioningSolution<idInput_t> > solutionConst =
    rcp_const_cast<const Zoltan2::PartitioningSolution<idInput_t> >(solution);

  // create metric object (also usually created by a problem)

  RCP<quality_t> metricObject;

  try{
    metricObject = rcp(new quality_t(env, comm, ia, solutionConst));
  }
  catch (std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "compute metrics", 1);


  if (rank==0){
    scalar_t imb;
    try{
      metricObject->getObjectCountImbalance(imb); 
      cout << "Object imbalance: " << imb << endl;
    }
    catch (std::exception &e){
      fail=1;
    }
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "getObjectCountImbalance", 1);

  if (rank==0 && weightDim > 0){
    scalar_t imb;
    try{
      for (int i=0; i < weightDim; i++){
        metricObject->getWeightImbalance(imb, i);
        cout << "Weight " << i << " imbalance: " << imb << endl;
      }
    }
    catch (std::exception &e){
      fail=10;
    }
    if (!fail && weightDim > 1){
      try{
        metricObject->getNormedImbalance(imb);
        cout << "Normed weight imbalance: " << imb << endl;
      }
      catch (std::exception &e){
        fail=11;
      }
    }
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "get imbalances", 1);
  
  if (rank==0){
    try{
      metricObject->printMetrics(cout);
    }
    catch (std::exception &e){
      fail=1;
    }
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "print metrics", 1);
}
