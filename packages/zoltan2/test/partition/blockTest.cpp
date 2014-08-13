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
#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>

using namespace std;
using Teuchos::Comm;
using Teuchos::RCP;

int main(int narg, char **arg)
{
  int fail=0, gfail=0;
  Teuchos::GlobalMPISession session(&narg, &arg);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  // Construct the user input
  int numGlobalIdentifiers = 100;
  int numMyIdentifiers = numGlobalIdentifiers / nprocs;
  if (rank < numGlobalIdentifiers % nprocs)
    numMyIdentifiers += 1;

  gno_t myBaseId = gno_t(numGlobalIdentifiers * rank);

  gid_t *myIds = new gid_t[numMyIdentifiers];
  scalar_t *myWeights = new scalar_t[numMyIdentifiers];

  if (!myIds || !myWeights){
    fail = 1;
  }

  gfail = globalFail(comm, fail);

  if (gfail){
    if (rank==0){
      std::cout << "Memory allocation failure" << std::endl;
      std::cout << "FAIL" << std:: endl;
    }
    return 1;
  }

  scalar_t origsumwgts = 0;
  for (int i=0; i < numMyIdentifiers; i++){
    myIds[i] = myBaseId+i;
    myWeights[i] = rank%3 + 1;
    origsumwgts += myWeights[i];
  }

  // Some output
  int *origcnt = new int[nprocs];
  scalar_t *origwgts = new scalar_t[nprocs];
  Teuchos::gather<int, int>(&numMyIdentifiers, 1, origcnt, 1, 0, *comm);
  Teuchos::gather<int, scalar_t>(&origsumwgts, 1, origwgts, 1, 0, *comm);
  if (rank == 0) {
    cout << "BEFORE PART CNTS: ";
    for (int i = 0; i < nprocs; i++)
      cout << origcnt[i] << " ";
    cout << endl;
    cout << "BEFORE PART WGTS: ";
    for (int i = 0; i < nprocs; i++)
      cout << origwgts[i] << " ";
    cout << endl;
  }
  delete [] origcnt;
  delete [] origwgts;

  // Building Zoltan2 adapters
  std::vector<const scalar_t *> weightValues;
  std::vector<int> weightStrides;   // default is one
  weightValues.push_back(const_cast<const scalar_t *>(myWeights));

  typedef Zoltan2::BasicUserTypes<scalar_t, gid_t, lno_t, gno_t> mydata_t;
  typedef Zoltan2::BasicIdentifierAdapter<mydata_t> adapter_t;
  typedef adapter_t::part_t part_t;

  adapter_t adapter(lno_t(numMyIdentifiers),myIds,weightValues,weightStrides);

  // Set up the parameters and problem
  bool useWeights = true;
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("weights", "no-weights", &useWeights,
                "Indicated whether to use identifier weights in partitioning");
  cmdp.parse(narg, arg);

  Teuchos::ParameterList params("test parameters");
  params.set("compute_metrics", "true");
  params.set("num_global_parts", nprocs);
  params.set("algorithm", "block");
  params.set("partitioning_approach", "partition");
  if (!useWeights) params.set("partitioning_objective", "balance_object_count");
  
  Zoltan2::PartitioningProblem<adapter_t> problem(&adapter, &params);

  problem.solve();

  Zoltan2::PartitioningSolution<adapter_t> solution = problem.getSolution();

  // Some output 
  scalar_t *totalWeight = new scalar_t [nprocs];
  scalar_t *sumWeight = new scalar_t [nprocs];
  memset(totalWeight, 0, nprocs * sizeof(scalar_t));
  int *totalCnt = new int [nprocs];
  int *sumCnt = new int [nprocs];
  memset(totalCnt, 0, nprocs * sizeof(int));

  const gid_t *idList = solution.getIdList();
  const part_t *partList = solution.getPartList();
  const scalar_t libImbalance = problem.getWeightImbalance();

  for (int i=0; !fail && i < numMyIdentifiers; i++){
    if (idList[i] != myIds[i])
      fail = 1;

    if (!fail) {
      totalCnt[partList[i]]++;
      totalWeight[partList[i]] += myWeights[i];
    }
  }

  gfail = globalFail(comm, fail);

  if (gfail){
    if (rank==0){
      std::cout << "failure in solution data" << std::endl;
      std::cout << "FAIL" << std:: endl;
    }
    return 1;
  }

  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, nprocs, 
                               totalCnt, sumCnt);
  Teuchos::reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM, nprocs, 
                                    totalWeight, sumWeight);

  double epsilon = 10e-6;

  if (rank == 0){
    std::cout << "AFTER PART CNTS: ";
    for (int i=0; i < nprocs; i++)
      std::cout << sumCnt[i] << " ";
    std::cout << std::endl;

    scalar_t total = 0;
    std::cout << "AFTER PART WGTS: ";
    for (int i=0; i < nprocs; i++){
      std::cout << sumWeight[i] << " ";
      total += sumWeight[i];
    }
    std::cout << std::endl;

    scalar_t avg = total / scalar_t(nprocs);

    scalar_t imbalance = -1.0;

    for (int i=0; i < nprocs; i++){
      scalar_t imb = 0;
      if (sumWeight[i] > avg)
        imb = (sumWeight[i] - avg) / avg;
      else
        imb = (avg - sumWeight[i]) / avg;

      if (imb > imbalance)
        imbalance = imb;
    }
    imbalance += 1.0;

    std::cout << "Computed imbalance: " << imbalance << std::endl;
    std::cout << "Library's imbalance: " << libImbalance << std::endl;

    double err;
    if (imbalance > libImbalance)
      err = imbalance - libImbalance;
    else
      err = libImbalance - imbalance;

    if (err > epsilon)
      fail = 1;
  }
  else{
      fail = 0;
  }

  gfail = globalFail(comm, fail);

  if (gfail){
    if (rank==0){
      std::cout << "failure in solution's imbalance data" << std::endl;
      std::cout << "FAIL" << std:: endl;
    }
    return 1;
  }

  if (rank==0)
    std::cout << "PASS" << std:: endl;

  delete [] myWeights;
  delete [] myIds;
  delete [] sumCnt;
  delete [] totalCnt;
  delete [] sumWeight;
  delete [] totalWeight;

  return 0;
}
