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
  int nWeights, int numLocalParts, bool givePartSizes);

typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> user_t;

// for testing basic input adapter
typedef Zoltan2::BasicIdentifierAdapter<user_t> basic_idInput_t;

// for testing graph adapter
typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tcrsGraph_t;
typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t, user_t> graph_idInput_t;

// creates this so we can run the test suite over BasicIdentifierAdapter
// and XpetraCrsGraphAdapter
template<class idInput_t> void runTestSuite(RCP<const Comm<int> > comm) {
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
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();

  // do some tests with BasicIdentifierAdapter
  runTestSuite<basic_idInput_t>(comm);

  // now some tests with XpetraCrsGraphAdapter
  // Note that right now these are all going to produce the same graph
  // metrics but could be developed further
  runTestSuite<graph_idInput_t>(comm);
  
  if (rank==0)
    cout << "PASS" << endl;
}

template<class idInput_t>
int evaluate_result(RCP<const Comm<int> > comm,
  RCP<Zoltan2::EvaluatePartition<idInput_t>> metric,
  int numLocalObjs, int numGlobalParts) {
    throw std::logic_error("evaluate_result not implemented.");
}

template<>
int evaluate_result<graph_idInput_t>(RCP<const Comm<int> > comm,
  RCP<Zoltan2::EvaluatePartition<graph_idInput_t>> metric,
  int numLocalObjs, int numGlobalParts) {

  int fail = 0;

  // first just confirm we can read getTotalEdgeCut()
  if(comm->getRank() == 0) { // getTotalEdgeCut() is global, just check rank 0
    int total_edge_cut;
    try{
      // TODO: the unweighted getTotalEdgeCut is an integer
      // maybe the API should be changed for this and other similar cases
      total_edge_cut = static_cast<int>(metric->getTotalEdgeCut());
      cout << "Total Edge Cut: " << total_edge_cut << endl;
    }
    catch (std::exception &e){
      fail=1;
    }

    int max_edge_cut;
    try{
      max_edge_cut = static_cast<int>(metric->getMaxEdgeCut());
      cout << "Max Edge Cut: " << max_edge_cut << endl;
    }
    catch (std::exception &e){
      fail=1;
    }

    int total_messages;
    try{
      total_messages = static_cast<int>(metric->getTotalMessages());
      cout << "Total Messages: " << total_messages << endl;
    }
    catch (std::exception &e){
      fail=1;
    }

    int max_messages;
    try{
      max_messages = static_cast<int>(metric->getMaxMessages());
      cout << "Max Messages: " << max_messages << endl;
    }
    catch (std::exception &e){
      fail=1;
    }

    // Now let's check our numbers
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
    // also handle edge case of num_procs less than 2 - no edge cuts then
    int expected_total_edge_cuts = (num_procs == 1) ? 0 :
      (2 * numLocalObjs) + ((num_procs-2) * numLocalObjs * 2);

    // now we validate
    if(total_edge_cut != expected_total_edge_cuts) {
      printf("We expected total edge cuts to be %d but got %d. ",
        expected_total_edge_cuts, total_edge_cut);
      fail = 1;
    }

    // we can also calculate max edge cuts
    // if num_procs 1, then it's 0
    // if num_procs 2, then it's numLocalObjs
    // otherwise it's 2 * numLocalObjs because at least one proc is sending
    // to two other procs
    int expected_max_edge_cuts = (num_procs == 1) ? 0 :
      (num_procs == 2) ? numLocalObjs : numLocalObjs * 2;

    // now we validate
    if(max_edge_cut != expected_max_edge_cuts) {
      printf("We expected max edge cuts to be %d but got %d. ",
        expected_total_edge_cuts, total_edge_cut);
      fail = 1;
    }

    // now check total messages - in present form we can simply divide but in
    // future things could be generalized
    int expected_total_messages = expected_total_edge_cuts / numLocalObjs;
    // now we validate
    if(total_messages != expected_total_messages) {
      printf("We expected total messages to be %d but got %d. ",
        expected_total_messages, total_messages);
      fail = 1;
    }

    // now check max messages - in present form we can simply divide but in
    // future things could be more generalized
    int expected_max_messages = expected_max_edge_cuts / numLocalObjs;
    // now we validate
    if(max_messages != expected_max_messages) {
      printf("We expected max messages to be %d but got %d. ",
        expected_max_messages, max_messages);
      fail = 1;
    }
  }
  
  return fail;
}

template<>
int evaluate_result<basic_idInput_t>(RCP<const Comm<int> > comm,
  RCP<Zoltan2::EvaluatePartition<basic_idInput_t>> metric,
  int numLocalObjs, int numGlobalParts) {

  // TODO - I set evaluate_result up to implement some validation of the
  // graph metrics utility methods. However we could extend this to imbalance
  // metrics. Just want to review the setup before going to far.
  return 0;
}

template<class idInput_t>
idInput_t * create_adapter(RCP<const Comm<int> > comm,
  int numLocalObj, zgno_t *myGids,
  std::vector<const zscalar_t *> & weights,
  std::vector<int> & strides) {
    throw std::logic_error("create_adapter not implemented.");
}

template<>
graph_idInput_t * create_adapter<graph_idInput_t>(RCP<const Comm<int> > comm,
  int numLocalObj, zgno_t *myGids,
  std::vector<const zscalar_t *> & weights,
  std::vector<int> & strides) {

  typedef Tpetra::Map<zlno_t, zgno_t> map_t;
  typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t> matrix_t;

  const zgno_t gNvtx = numLocalObj * comm->getSize();
  const Teuchos::ArrayView<const zgno_t> indexList(myGids, numLocalObj);
  Teuchos::RCP<const map_t> map = rcp(new map_t(gNvtx, indexList, 0, comm));

  // Make some stuff in the graph
  size_t maxRowLen = 1;
  Teuchos::RCP<matrix_t> matrix = rcp(new matrix_t(map, maxRowLen));
  
  // TODO: Need to decide what a good test graph is
  // Something we can easily calculate the final result for as we'd like to
  // validate this but not end up rewriting the algorithm we are testing.
  // For starters I simply have each graph element create two links to the
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

  return adapter;
}

template<>
basic_idInput_t * create_adapter<basic_idInput_t>(RCP<const Comm<int> > comm,
  int numLocalObj, zgno_t *myGids,
  std::vector<const zscalar_t *> & weights,
  std::vector<int> & strides) {
  return new basic_idInput_t(numLocalObj, myGids, weights, strides);
}

// Assumes numLocalObj is the same on every process.
template<class idInput_t>
void doTest(RCP<const Comm<int> > comm, int numLocalObj,
  int nWeights, int numLocalParts, bool givePartSizes)
{
  typedef Zoltan2::EvaluatePartition<idInput_t> quality_t;

  typedef typename idInput_t::part_t part_t;

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
    cout << "Test: number of weights " << nWeights;
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
    ia = create_adapter<idInput_t>(comm, numLocalObj, myGids, weights, strides);
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


  if (rank==0){
    ;
    try{
      zscalar_t imb = metricObject->getObjectCountImbalance();
      cout << "Object imbalance: " << imb << endl;
    }
    catch (std::exception &e){
      fail=1;
    }
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "getObjectCountImbalance", 1);

  if (rank==0 && nWeights > 0){
    try{
      for (int i=0; i < nWeights; i++){
    	zscalar_t imb = metricObject->getWeightImbalance(i);
        cout << "Weight " << i << " imbalance: " << imb << endl;
      }
    }
    catch (std::exception &e){
      fail=10;
    }
    if (!fail && nWeights > 1){
      try{
    	zscalar_t imb = metricObject->getNormedImbalance();
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
  
  // Now evaluate the result for precision of the results
  // TODO: I'd like to discuss these ideas further before continuing
  // To prototype this, I'll just validate the graph adapter getMaxEdgeCut()
  // method.
  fail = evaluate_result<idInput_t>(comm, metricObject,
    numLocalObj, numGlobalParts);
  TEST_FAIL_AND_EXIT(*comm, fail==0, "evaluate_result", 1);

  delete ia;
}
