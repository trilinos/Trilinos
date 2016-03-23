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

/*! \file rcb_C.cpp
    \brief An example of partitioning coordinates with RCB.
*/

#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Tpetra_Map.hpp>
#include <vector>
#include <cstdlib>

using namespace std;
using std::vector;
using Teuchos::RCP;
using Zoltan2::Environment;

/*! \example rcb_C.cpp
    An example of the use of the RCB algorithm to partition coordinate data.
*/

int main(int argc, char *argv[])
{
#ifdef HAVE_ZOLTAN2_MPI                   
  MPI_Init(&argc, &argv);
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int rank=0, nprocs=1;
#endif

  // For convenience, we'll use the Tpetra defaults for local/global ID types
  // Users can substitute their preferred local/global ID types
  typedef Tpetra::Map<> Map_t;
  typedef Map_t::local_ordinal_type localId_t;
  typedef Map_t::global_ordinal_type globalId_t;

  typedef double scalar_t;
  typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;

  // TODO explain
  typedef Zoltan2::BasicVectorAdapter<myTypes> inputAdapter_t;
  typedef Zoltan2::EvaluatePartition<inputAdapter_t, myTypes> quality_t;
  typedef inputAdapter_t::part_t part_t;
  typedef inputAdapter_t::base_adapter_t base_adapter_t;

  ///////////////////////////////////////////////////////////////////////
  // Create input data.

  size_t localCount = 40;
  int dim = 3;

  scalar_t *coords = new scalar_t [dim * localCount];

  scalar_t *x = coords; 
  scalar_t *y = x + localCount; 
  scalar_t *z = y + localCount; 

  // Create coordinates that range from 0 to 10.0

  srand(rank);
  scalar_t scalingFactor = 10.0 / RAND_MAX;

  for (size_t i=0; i < localCount*dim; i++){
    coords[i] = scalar_t(rand()) * scalingFactor;
  }

  // Create global ids for the coordinates.

  globalId_t *globalIds = new globalId_t [localCount];
  globalId_t offset = rank * localCount;

  for (size_t i=0; i < localCount; i++)
    globalIds[i] = offset++;
   
  ///////////////////////////////////////////////////////////////////////
  // Create parameters for an RCB problem

  double tolerance = 1.1;

  if (rank == 0)
    std::cout << "Imbalance tolerance is " << tolerance << std::endl;

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("debug_procs", "0");
  params.set("error_check_level", "debug_mode_assertions");

  params.set("algorithm", "rcb");
  params.set("imbalance_tolerance", tolerance );
  params.set("num_global_parts", nprocs);

  params.set("bisection_num_test_cuts", 1);
   
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // A simple problem with no weights.
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // Create a Zoltan2 input adapter for this geometry. TODO explain

  inputAdapter_t *ia1 = new inputAdapter_t(localCount,globalIds,x,y,z,1,1,1);

  // Create a Zoltan2 partitioning problem

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem1 =
           new Zoltan2::PartitioningProblem<inputAdapter_t>(ia1, &params);
   
  // Solve the problem

  problem1->solve();

  // An environment.  This is usually created by the problem.
  // Note:  These RCPs will go away in Spring 2016 when we finish simplication
  // of the EvaluatePartition interface.

  RCP<const Environment> env1 = problem1->getEnvironment();

  // create metric object

  quality_t *metricObject1 = new quality_t(env1, problem1->getComm(), ia1,
					   &problem1->getSolution());
  // Check the solution.

  if (rank == 0) {
    metricObject1->printMetrics(cout);
  }

  if (rank == 0){
    scalar_t imb;
    metricObject1->getWeightImbalance(imb, 0);
    if (imb <= tolerance)
      std::cout << "pass: " << imb << std::endl;
    else
      std::cout << "fail: " << imb << std::endl;
    std::cout << std::endl;
  }
  delete metricObject1;
   
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // Try a problem with weights 
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  scalar_t *weights = new scalar_t [localCount];
  for (size_t i=0; i < localCount; i++){
    weights[i] = 1.0 + scalar_t(rank) / scalar_t(nprocs);
  }

  // Create a Zoltan2 input adapter that includes weights.

  vector<const scalar_t *>coordVec(2);
  vector<int> coordStrides(2);

  coordVec[0] = x; coordStrides[0] = 1;
  coordVec[1] = y; coordStrides[1] = 1;

  vector<const scalar_t *>weightVec(1);
  vector<int> weightStrides(1);

  weightVec[0] = weights; weightStrides[0] = 1;

  inputAdapter_t *ia2=new inputAdapter_t(localCount, globalIds, coordVec, 
                                         coordStrides,weightVec,weightStrides);

  // Create a Zoltan2 partitioning problem

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem2 =
           new Zoltan2::PartitioningProblem<inputAdapter_t>(ia2, &params);

  // Solve the problem

  problem2->solve();

  // An environment.  This is usually created by the problem.
  // Note:  These RCPs will go away in Spring 2016 when we finish simplication
  // of the EvaluatePartition interface.

  RCP<const Environment> env2 = problem2->getEnvironment();

  // create metric object

  quality_t *metricObject2 = new quality_t(env2, problem2->getComm(), ia2,
					   &problem2->getSolution());
  // Check the solution.

  if (rank == 0) {
    metricObject2->printMetrics(cout);
  }

  if (rank == 0){
    scalar_t imb;
    metricObject2->getWeightImbalance(imb, 0);
    if (imb <= tolerance)
      std::cout << "pass: " << imb << std::endl;
    else
      std::cout << "fail: " << imb << std::endl;
    std::cout << std::endl;
  }
  delete metricObject2;

  if (localCount > 0){
    delete [] weights;
    weights = NULL;
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // Try a problem with multiple weights.
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // Add to the parameters the multicriteria objective.

  params.set("partitioning_objective", "multicriteria_minimize_total_weight");

  // Create the new weights.

  weights = new scalar_t [localCount*3];
  srand(rank);

  for (size_t i=0; i < localCount*3; i+=3){
    weights[i] = 1.0 + rank / nprocs;      // weight idx 1
    weights[i+1] = rank<nprocs/2 ? 1 : 2;  // weight idx 2
    weights[i+2] = rand()/RAND_MAX +.5;    // weight idx 3
  }

  // Create a Zoltan2 input adapter with these weights.

  weightVec.resize(3);
  weightStrides.resize(3);

  weightVec[0] = weights;   weightStrides[0] = 3;
  weightVec[1] = weights+1; weightStrides[1] = 3;
  weightVec[2] = weights+2; weightStrides[2] = 3;

  inputAdapter_t *ia3=new inputAdapter_t(localCount, globalIds, coordVec,
                                         coordStrides,weightVec,weightStrides);

  // Create a Zoltan2 partitioning problem.

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem3 =
           new Zoltan2::PartitioningProblem<inputAdapter_t>(ia3, &params);

  // Solve the problem

  problem3->solve();

  // An environment.  This is usually created by the problem.
  // Note:  These RCPs will go away in Spring 2016 when we finish simplication
  // of the EvaluatePartition interface.

  RCP<const Environment> env3 = problem3->getEnvironment();

  // create metric object

  quality_t *metricObject3 = new quality_t(env3, problem3->getComm(), ia3,
					   &problem3->getSolution());
  // Check the solution.

  if (rank == 0) {
    metricObject3->printMetrics(cout);
  }

  if (rank == 0){
    scalar_t imb;
    metricObject3->getWeightImbalance(imb, 0);
    if (imb <= tolerance)
      std::cout << "pass: " << imb << std::endl;
    else
      std::cout << "fail: " << imb << std::endl;
    std::cout << std::endl;
  }
  delete metricObject3;

  ///////////////////////////////////////////////////////////////////////
  // Try the other multicriteria objectives.

  bool dataHasChanged = false;    // default is true

  params.set("partitioning_objective", "multicriteria_minimize_maximum_weight");
  problem3->resetParameters(&params);
  problem3->solve(dataHasChanged);    

  // Objective changed!

  env3 = problem3->getEnvironment();

  // Solution changed!

  metricObject3 = new quality_t(env3, problem3->getComm(), ia3,
                                &problem3->getSolution());
  if (rank == 0){
    metricObject3->printMetrics(cout);
    scalar_t imb;
    metricObject3->getWeightImbalance(imb, 0);
    if (imb <= tolerance)
      std::cout << "pass: " << imb << std::endl;
    else
      std::cout << "fail: " << imb << std::endl;
    std::cout << std::endl;
  }
  delete metricObject3;

  params.set("partitioning_objective", "multicriteria_balance_total_maximum");
  problem3->resetParameters(&params);
  problem3->solve(dataHasChanged);    

  // Objective changed!

  env3 = problem3->getEnvironment();

  // Solution changed!

  metricObject3 = new quality_t(env3, problem3->getComm(), ia3,
                                &problem3->getSolution());
  if (rank == 0){
    metricObject3->printMetrics(cout);
    scalar_t imb;
    metricObject3->getWeightImbalance(imb, 0);
    if (imb <= tolerance)
      std::cout << "pass: " << imb << std::endl;
    else
      std::cout << "fail: " << imb << std::endl;
    std::cout << std::endl;
  }
  delete metricObject3;

  if (localCount > 0){
    delete [] weights;
    weights = NULL;
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // Using part sizes, ask for some parts to be empty.
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // Change the number of parts to twice the number of processes to
  // ensure that we have more than one global part.

  params.set("num_global_parts", nprocs*2);

  // Using the initial problem that did not have any weights, reset
  // parameter list, and give it some part sizes.

  problem1->resetParameters(&params);

  part_t partIds[2];
  scalar_t partSizes[2];

  partIds[0] = rank*2;    partSizes[0] = 0;
  partIds[1] = rank*2+1;  partSizes[1] = 1;

  problem1->setPartSizes(2, partIds, partSizes);

  // Solve the problem.  The argument "dataHasChanged" indicates 
  // that we have not changed the input data, which allows the problem
  // so skip some work when re-solving. 

  dataHasChanged = false;

  problem1->solve(dataHasChanged);

  // Obtain the solution

  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution4 =
    problem1->getSolution();

  // Check it.  Part sizes should all be odd.

  const part_t *partAssignments = solution4.getPartListView();

  int numInEmptyParts = 0;
  for (size_t i=0; i < localCount; i++){
    if (partAssignments[i] % 2 == 0)
      numInEmptyParts++;
  }

  if (rank == 0)
    std::cout << "Request that " << nprocs << " parts be empty." <<std::endl;

  // Solution changed!

  metricObject1 = new quality_t(env1, problem1->getComm(), ia1,
                                &problem1->getSolution());
  // Check the solution.

  if (rank == 0) {
    metricObject1->printMetrics(cout);
  }

  if (rank == 0){
    scalar_t imb;
    metricObject1->getWeightImbalance(imb, 0);
    if (imb <= tolerance)
      std::cout << "pass: " << imb << std::endl;
    else
      std::cout << "fail: " << imb << std::endl;
    std::cout << std::endl;
  }
  delete metricObject1;

  if (coords)
    delete [] coords;

  if (globalIds)
    delete [] globalIds;

  delete problem1;
  delete ia1;
  delete problem2;
  delete ia2;
  delete problem3;
  delete ia3;

#ifdef HAVE_ZOLTAN2_MPI
  MPI_Finalize();
#endif

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

