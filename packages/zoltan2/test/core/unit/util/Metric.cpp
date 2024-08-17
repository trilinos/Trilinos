// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Test the following:
//         EvaluatePartition class
//         MetricValues class
//         Metric related namespace methods


#include <Zoltan2_EvaluatePartition.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <stdlib.h>
#include <vector>


using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::arcp;

using namespace std;
using std::endl;
using std::cout;

template<class idInput_t>
void doTest(RCP<const Comm<int> > comm, int numLocalObj,
  int nWeights, int numLocalParts, bool givePartSizes, bool useDegreeAsWeight=false);

typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> user_t;

// for testing basic input adapter
typedef Zoltan2::BasicIdentifierAdapter<user_t> basic_idInput_t;

// for testing graph adapter
typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tcrsGraph_t;
typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t, user_t> graph_idInput_t;

// creates this so we can run the test suite over BasicIdentifierAdapter
// and XpetraCrsGraphAdapter
template<class idInput_t> void runTestSuite(RCP<const Comm<int> > comm, bool bCanTestDegreeAsWeights) {
  doTest<idInput_t>(comm, 10, 0, -1, false);
  doTest<idInput_t>(comm, 10, 0,  1, false);
  doTest<idInput_t>(comm, 10, 0,  1, true);
  doTest<idInput_t>(comm, 10, 1,  1, false);
  doTest<idInput_t>(comm, 10, 1,  1, true);
  doTest<idInput_t>(comm, 10, 2,  1, false);
  doTest<idInput_t>(comm, 10, 2,  1, true);
  doTest<idInput_t>(comm, 10, 1,  2, true);
  doTest<idInput_t>(comm, 10, 1,  2, false);
  doTest<idInput_t>(comm, 10, 1, -1, false);
  doTest<idInput_t>(comm, 10, 1, -1, true);
  doTest<idInput_t>(comm, 10, 2, -1, false);

  if(bCanTestDegreeAsWeights) {
    doTest<idInput_t>(comm, 10, 1, 1, true, true); // with degreeAsWeights
  }
}

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();

  // do some tests with BasicIdentifierAdapter
  runTestSuite<basic_idInput_t>(comm, false);

  // now some tests with XpetraCrsGraphAdapter
  // Note that right now these are all going to produce the same graph
  // metrics but could be developed further
  runTestSuite<graph_idInput_t>(comm, true);

  comm->barrier();
  if (rank==0)
    std::cout << "PASS" << std::endl;
}

// to validate the results, we call evaluate_adapter_results which is
// templated so it can, for example, check graph metrics only for the graph
// adapter. Currently both basic and graph adapter setup imbalance metrics so
// we also do a check on that with a universal call to
// evaluate_imbalance_results. Currently this needs no specialization.
// If we add more adapters later this could be flexible to accomodate that.
template<class idInput_t>
void evaluate_imbalance_results(RCP<const Comm<int> > comm,
  RCP<Zoltan2::EvaluatePartition<idInput_t>> metricObject, int numLocalObj,
  int nWeights, int original_numLocalParts, bool givePartSizes) {
  int fail = 0;
  int rank = comm->getRank();

  zscalar_t object_count_imbalance;
  try{
    object_count_imbalance = metricObject->getObjectCountImbalance();
    if(rank == 0) {
      cout << "Object imbalance: " << object_count_imbalance << endl;
    }
  }
  catch (std::exception &e){
    fail=1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getObjectCountImbalance", 1);

  if (nWeights > 0){
    try{
      for (int i=0; i < nWeights; i++){
        zscalar_t imb = metricObject->getWeightImbalance(i);
        if(rank == 0){
          cout << "Weight " << i << " imbalance: " << imb << endl;
        }
      }
    }
    catch (std::exception &e){
      fail=10;
    }
    if (!fail && nWeights > 1){
      try{
        zscalar_t imb = metricObject->getNormedImbalance();
        if(rank == 0){
          cout << "Normed weight imbalance: " << imb << endl;
        }
      }
      catch (std::exception &e){
        fail=11;
      }
    }
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "get imbalances", 1);
}

template<class idInput_t>
void evaluate_adapter_results(RCP<const Comm<int> > comm,
  RCP<Zoltan2::EvaluatePartition<idInput_t>> metricObject, int numLocalObj,
  int nWeights, int original_numLocalParts, bool givePartSizes) {
    throw std::logic_error("evaluate_result not implemented.");
}

template<>
void evaluate_adapter_results<graph_idInput_t>(RCP<const Comm<int> > comm,
  RCP<Zoltan2::EvaluatePartition<graph_idInput_t>> metricObject, int numLocalObj,
  int nWeights, int original_numLocalParts, bool givePartSizes) {
  int fail = 0;
  int rank = comm->getRank();

  int total_edge_cut = -1;
  try{
    // TODO: the unweighted getTotalEdgeCut is an integer
    // maybe the API should be changed for this and other similar cases
    total_edge_cut = static_cast<int>(metricObject->getTotalEdgeCut());
    if(rank == 0){
      cout << "Total Edge Cut: " << total_edge_cut << endl;
    }
  }
  catch (std::exception &e){
    fail=1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getTotalEdgeCut", 1);

  int max_edge_cut = -1;
  try{
    max_edge_cut = static_cast<int>(metricObject->getMaxEdgeCut());
    if(rank == 0){
      cout << "Max Edge Cut: " << max_edge_cut << endl;
    }
  }
  catch (std::exception &e){
    fail=1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getMaxEdgeCut", 1);

  int total_messages = -1;
  try{
    total_messages = static_cast<int>(metricObject->getTotalMessages());
    if(rank == 0){
      cout << "Total Messages: " << total_messages << endl;
    }
  }
  catch (std::exception &e){
    fail=1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getTotalMessages", 1);

  int max_messages = -1;
  try{
    max_messages = static_cast<int>(metricObject->getMaxMessages());
    if(rank == 0){
      cout << "Max Messages: " << max_messages << endl;
    }
  }
  catch (std::exception &e){
    fail=1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getMaxMessages", 1);

  // Now let's check our numbers.
  // Here we do a calculation of what getTotalEdgeCut should return based on
  // how we set things up in create_adapter.
  // Currently the algorithm simply has every object create two links, one
  // to the first global id and one to the last.
  // Two of the procs will contain one of those global ids so they only have
  // edge cuts equal to numLocalObjs to send to the other.
  // So that is the (2 * numLocalObjs) term.
  // All other procs will contain neither of those global ids so they have
  // to send their objects to two procs.
  // So that is the ((num_procs-2) * numLocalObjs * 2 term.
  int num_procs = comm->getSize();
  int expected_total_edge_cuts = (num_procs == 1) ? 0 :
    (2 * numLocalObj) + ((num_procs-2) * numLocalObj * 2);
  TEST_FAIL_AND_EXIT(*comm, total_edge_cut == expected_total_edge_cuts,
    "getTotalEdgeCut is not the expected ", 1);

  // we can also calculate max edge cuts
  // if num_procs 1, then it's 0
  // if num_procs 2, then it's numLocalObjs
  // otherwise it's 2 * numLocalObjs because at least one proc is sending
  // to two other procs
  int expected_max_edge_cuts = (num_procs == 1) ? 0 :
    (num_procs == 2) ? numLocalObj : numLocalObj * 2;
  TEST_FAIL_AND_EXIT(*comm, max_edge_cut == expected_max_edge_cuts,
      "getMaxEdgeCut is not the expected value", 1);

  // now check total messages - in present form we can simply divide but in
  // future things could be generalized
  int expected_total_messages = expected_total_edge_cuts / numLocalObj;
  TEST_FAIL_AND_EXIT(*comm, total_messages == expected_total_messages,
      "getTotalMessages is not the expected value", 1);

  // now check max messages - in present form we can simply divide but in
  // future things could be more generalized
  int expected_max_messages = expected_max_edge_cuts / numLocalObj;
  TEST_FAIL_AND_EXIT(*comm, max_messages == expected_max_messages,
      "getMaxMessages is not the expected value", 1);

  evaluate_imbalance_results(comm, metricObject,
    numLocalObj, nWeights, original_numLocalParts, givePartSizes);
}

// for basic_idInput_t we just call the common evaluate_imbalance_results
// no other specialized data to consider
template<>
void evaluate_adapter_results<basic_idInput_t>(RCP<const Comm<int> > comm,
  RCP<Zoltan2::EvaluatePartition<basic_idInput_t>> metricObject, int numLocalObj,
  int nWeights, int original_numLocalParts, bool givePartSizes) {
  evaluate_imbalance_results(comm, metricObject,
    numLocalObj, nWeights, original_numLocalParts, givePartSizes);
}

template<class idInput_t>
idInput_t * create_adapter(RCP<const Comm<int> > comm,
  int numLocalObj, zgno_t *myGids,
  std::vector<const zscalar_t *> & weights,
  std::vector<int> & strides,
  bool useDegreeAsWeight) {
    throw std::logic_error("create_adapter not implemented.");
}

template<>
graph_idInput_t * create_adapter<graph_idInput_t>(RCP<const Comm<int> > comm,
  int numLocalObj, zgno_t *myGids,
  std::vector<const zscalar_t *> & weights,
  std::vector<int> & strides,
  bool useDegreeAsWeight) {

  typedef Tpetra::Map<zlno_t, zgno_t> map_t;
  typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t> matrix_t;

  const zgno_t gNvtx = numLocalObj * comm->getSize();
  const Teuchos::ArrayView<const zgno_t> indexList(myGids, numLocalObj);
  Teuchos::RCP<const map_t> map = rcp(new map_t(gNvtx, indexList, 0, comm));

  // Make some stuff in the graph
  size_t maxRowLen = 2;
  Teuchos::RCP<matrix_t> matrix = rcp(new matrix_t(map, maxRowLen));

  // I picked this graph as a simple test case.
  // Something we can easily calculate the final result for as we'd like to
  // validate this but not end up rewriting the algorithm we are testing.
  // I have each graph element create two links to the
  // first global index and last. That means two procs will have edge cuts
  // equal to their numLocalObj while the rest will have 2 * numLocalObj
  //
  // Two of the procs will have only 1 message to send
  // The other procs will have 2 messages to send
  // Message max is 2
  // Message total is going to be (2)*2 + (numProcs-2)*2
  Teuchos::Array<zgno_t> col(2);
  Teuchos::Array<zscalar_t> val(2); val[0] = 1.; val[1] = 1.;
  zgno_t first_id = map->getMinAllGlobalIndex();
  zgno_t last_id = map->getMaxAllGlobalIndex();
  for (zlno_t i = 0; i < numLocalObj; i++) {
    zgno_t id = map->getGlobalElement(i);
    col[0] = first_id;
    col[1] = last_id;
    matrix->insertGlobalValues(id, col(), val());
  }

  matrix->fillComplete(map, map);

  size_t nVwgts = weights.size();
  graph_idInput_t * adapter = new graph_idInput_t(matrix->getCrsGraph(), nVwgts);

  // Set the weights
  for (size_t j = 0; j < nVwgts; j++) {
    adapter->setWeights(weights[j], 1, j);
  }

  // set degreeAsWeight if enabled
  if(useDegreeAsWeight) {
    for (size_t j = 0; j < nVwgts; j++) {
      adapter->setWeightIsDegree(j);
    }
  }

  return adapter;
}

template<>
basic_idInput_t * create_adapter<basic_idInput_t>(RCP<const Comm<int> > comm,
  int numLocalObj, zgno_t *myGids,
  std::vector<const zscalar_t *> & weights,
  std::vector<int> & strides,
  bool useDegreeAsWeight) {
  // useDegreeAsWeight is ignored
  return new basic_idInput_t(numLocalObj, myGids, weights, strides);
}

// Assumes numLocalObj is the same on every process.
template<class idInput_t>
void doTest(RCP<const Comm<int> > comm, int numLocalObj,
  int nWeights, int numLocalParts, bool givePartSizes, bool useDegreeAsWeight)
{
  typedef Zoltan2::EvaluatePartition<idInput_t> quality_t;

  typedef typename idInput_t::part_t part_t;

  int rank = comm->getRank();

  int original_numLocalParts = numLocalParts; // save for log and error checking

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

  if(rank == 0) {
    cout << endl
      << "Test: number of weights " << nWeights
      << ", desired number of parts " << numGlobalParts
      << ", calculated num local parts " << numLocalParts
      << ", original num local parts " << original_numLocalParts
      << (givePartSizes ? ", with differing part sizes." :
        ", with uniform part sizes.")
      << ", Number of procs " << nprocs
      << ", each with " << numLocalObj << " objects, part = rank."
      << (useDegreeAsWeight ? ", use degree as weights" : "")
      << endl;
  }

  // An environment.  This is usually created by the problem.

  Teuchos::ParameterList pl("test list");
  pl.set("num_local_parts", numLocalParts);
  
  RCP<const Zoltan2::Environment> env = 
    rcp(new Zoltan2::Environment(pl, comm));

  // A simple identifier map.  Usually created by the model.

  zgno_t *myGids = new zgno_t [numLocalObj];
  for (int i=0, x=rank*numLocalObj; i < numLocalObj; i++, x++){
    myGids[i] = x;
  }

  // Part sizes.  Usually supplied by the user to the Problem.
  // Then the problem supplies them to the Solution.

  int partSizeDim = (givePartSizes ? (nWeights ? nWeights : 1) : 0);
  ArrayRCP<ArrayRCP<part_t> > ids(partSizeDim);
  ArrayRCP<ArrayRCP<zscalar_t> > sizes(partSizeDim);

  if (givePartSizes && numLocalParts > 0){
    part_t *myParts = new part_t [numLocalParts];
    myParts[0] = rank * numLocalParts;
    for (int i=1; i < numLocalParts; i++)
      myParts[i] = myParts[i-1] + 1;
    ArrayRCP<part_t> partNums(myParts, 0, numLocalParts, true);

    zscalar_t sizeFactor = nprocs/2 - rank;
    if (sizeFactor < 0) sizeFactor *= -1;
    sizeFactor += 1;

    for (int dim=0; dim < partSizeDim; dim++){
      zscalar_t *psizes = new zscalar_t [numLocalParts];
      for (int i=0; i < numLocalParts; i++)
        psizes[i] = sizeFactor;
      sizes[dim] = arcp(psizes, 0, numLocalParts, true);
      ids[dim] = partNums;
    }
  }

  // An input adapter with random weights.  Created by the user.

  std::vector<const zscalar_t *> weights;
  std::vector<int> strides;   // default to 1

  int len = numLocalObj*nWeights;
  ArrayRCP<zscalar_t> wgtBuf;
  zscalar_t *wgts = NULL;

  if (len > 0){
    wgts = new zscalar_t [len];
    wgtBuf = arcp(wgts, 0, len, true);
    for (int i=0; i < len; i++)
      wgts[i] = (zscalar_t(rand()) / zscalar_t(RAND_MAX)) + 1.0;
  }

  for (int i=0; i < nWeights; i++, wgts+=numLocalObj)
    weights.push_back(wgts);

  idInput_t *ia = NULL;

  try {
    ia = create_adapter<idInput_t>(comm, numLocalObj, myGids, weights, strides, useDegreeAsWeight);
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
        env, comm, nWeights,
        ids.view(0,partSizeDim), sizes.view(0,partSizeDim)));
    else
      solution = rcp(new Zoltan2::PartitioningSolution<idInput_t>(
        env, comm, nWeights));
  }
  catch (std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "create solution", 1);

  // Part assignment for my objects: The algorithm usually calls this. 

  part_t *partNum = new part_t [numLocalObj];
  ArrayRCP<part_t> partAssignment(partNum, 0, numLocalObj, true);
  for (int i=0; i < numLocalObj; i++)
    partNum[i] = rank;

  solution->setParts(partAssignment);

  // create metric object (also usually created by a problem)

  RCP<quality_t> metricObject;

  try{
    metricObject = rcp(new quality_t(ia, &pl, comm, solution.getRawPtr()));
  }
  catch (std::exception &e){
    fail=1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "compute metrics", 1);
  
  try{
    if(rank == 0){
      metricObject->printMetrics(cout);
    }
  }
  catch (std::exception &e){
    fail=1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "print metrics", 1);

  // will call TEST_FAIL_AND_EXIT at each internal step
  evaluate_adapter_results<idInput_t>(comm, metricObject,
    numLocalObj, nWeights, original_numLocalParts, givePartSizes);

  delete ia;
}
