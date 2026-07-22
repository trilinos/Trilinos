// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>


int main(int narg, char **arg)
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int fail=0, gfail=0;

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  // Construct the user input
  int numGlobalIdentifiers = 100;
  int numMyIdentifiers = numGlobalIdentifiers / nprocs;
  if (rank < numGlobalIdentifiers % nprocs)
    numMyIdentifiers += 1;

  zgno_t myBaseId = zgno_t(numGlobalIdentifiers * rank);

  zgno_t *myIds = new zgno_t[numMyIdentifiers];
  zscalar_t *myWeights = new zscalar_t[numMyIdentifiers];

  if (!myIds || !myWeights){
    fail = 1;
  }

  gfail = globalFail(*comm, fail);

  if (gfail){
    if (rank==0){
      std::cout << "Memory allocation failure" << std::endl;
      std::cout << "FAIL" << std::endl;
    }
    return 1;
  }

  zscalar_t origsumwgts = 0;
  for (int i=0; i < numMyIdentifiers; i++){
    myIds[i] = myBaseId+i;
    myWeights[i] = rank%3 + 1;
    origsumwgts += myWeights[i];
  }

  // Some output
  int *origcnt = new int[nprocs];
  zscalar_t *origwgts = new zscalar_t[nprocs];
  Teuchos::gather<int, int>(&numMyIdentifiers, 1, origcnt, 1, 0, *comm);
  Teuchos::gather<int, zscalar_t>(&origsumwgts, 1, origwgts, 1, 0, *comm);
  if (rank == 0) {
    std::cout << "BEFORE PART CNTS: ";
    for (int i = 0; i < nprocs; i++)
      std::cout << origcnt[i] << " ";
    std::cout << std::endl;
    std::cout << "BEFORE PART WGTS: ";
    for (int i = 0; i < nprocs; i++)
      std::cout << origwgts[i] << " ";
    std::cout << std::endl;
  }
  delete [] origcnt;
  delete [] origwgts;

  // Building Zoltan2 adapters
  std::vector<const zscalar_t *> weightValues;
  std::vector<int> weightStrides;   // default is one
  weightValues.push_back(const_cast<const zscalar_t *>(myWeights));

  typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> mydata_t;
  typedef Zoltan2::BasicIdentifierAdapter<mydata_t> adapter_t;
  typedef Zoltan2::EvaluatePartition<adapter_t> quality_t;
  typedef adapter_t::part_t part_t;

  adapter_t *adapter = 
    new adapter_t(zlno_t(numMyIdentifiers),myIds,weightValues,weightStrides);

  // Set up the parameters and problem
  bool useWeights = true;
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("weights", "no-weights", &useWeights,
                "Indicated whether to use identifier weights in partitioning");
  cmdp.parse(narg, arg);

  Teuchos::ParameterList params("test parameters");
  //params.set("compute_metrics", true); // bool parameter
  params.set("num_global_parts", nprocs);
  params.set("algorithm", "block");
  params.set("partitioning_approach", "partition");
  if (!useWeights) params.set("partitioning_objective", "balance_object_count");
  
  Zoltan2::PartitioningProblem<adapter_t> problem(adapter, &params);

  problem.solve();

  Zoltan2::PartitioningSolution<adapter_t> solution = problem.getSolution();

  // create metric object

  quality_t *metricObject = new quality_t(adapter, &params, comm, &solution);

  // Some output 
  zscalar_t *totalWeight = new zscalar_t [nprocs];
  zscalar_t *sumWeight = new zscalar_t [nprocs];
  memset(totalWeight, 0, nprocs * sizeof(zscalar_t));
  int *totalCnt = new int [nprocs];
  int *sumCnt = new int [nprocs];
  memset(totalCnt, 0, nprocs * sizeof(int));

  const part_t *partList = solution.getPartListView();
  zscalar_t libImbalance = metricObject->getWeightImbalance(0);
  delete metricObject;
  delete adapter;

  for (int i=0; i < numMyIdentifiers; i++){
    totalCnt[partList[i]]++;
    totalWeight[partList[i]] += myWeights[i];
  }

  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, nprocs, 
                               totalCnt, sumCnt);
  Teuchos::reduceAll<int, zscalar_t>(*comm, Teuchos::REDUCE_SUM, nprocs, 
                                    totalWeight, sumWeight);

  double epsilon = 10e-6;

  if (rank == 0){
    std::cout << "AFTER PART CNTS: ";
    for (int i=0; i < nprocs; i++)
      std::cout << sumCnt[i] << " ";
    std::cout << std::endl;

    zscalar_t total = 0;
    std::cout << "AFTER PART WGTS: ";
    for (int i=0; i < nprocs; i++){
      std::cout << sumWeight[i] << " ";
      total += sumWeight[i];
    }
    std::cout << std::endl;

    zscalar_t avg = total / zscalar_t(nprocs);

    zscalar_t imbalance = -1.0;

    for (int i=0; i < nprocs; i++){
      zscalar_t imb = 0;
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

  gfail = globalFail(*comm, fail);

  if (gfail){
    if (rank==0){
      std::cout << "failure in solution's imbalance data" << std::endl;
      std::cout << "FAIL" << std::endl;
    }
    return 1;
  }

  if (rank==0)
    std::cout << "PASS" << std::endl;

  delete [] myWeights;
  delete [] myIds;
  delete [] sumCnt;
  delete [] totalCnt;
  delete [] sumWeight;
  delete [] totalWeight;

  return 0;
}
