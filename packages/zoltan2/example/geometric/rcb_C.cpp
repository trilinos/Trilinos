// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

/*! \example rcb_C.cpp
    An example of the use of the RCB algorithm to partition coordinate data.
*/

int main(int argc, char *argv[])
{
  Tpetra::ScopeGuard tscope(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  // For convenience, we'll use the Tpetra defaults for local/global ID types
  // Users can substitute their preferred local/global ID types
  typedef Tpetra::Map<> Map_t;
  typedef Map_t::local_ordinal_type localId_t;
  typedef Map_t::global_ordinal_type globalId_t;

  typedef Tpetra::Details::DefaultTypes::scalar_type scalar_t;
  typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;

  // TODO explain
  typedef Zoltan2::BasicVectorAdapter<myTypes> inputAdapter_t;
  typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;
  typedef inputAdapter_t::part_t part_t;

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

  // create metric object where communicator is Teuchos default

  quality_t *metricObject1 = new quality_t(ia1, &params, //problem1->getComm(),
					   &problem1->getSolution());
  // Check the solution.

  if (rank == 0) {
    metricObject1->printMetrics(std::cout);
  }

  if (rank == 0){
    scalar_t imb = metricObject1->getObjectCountImbalance();
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

  std::vector<const scalar_t *>coordVec(2);
  std::vector<int> coordStrides(2);

  coordVec[0] = x; coordStrides[0] = 1;
  coordVec[1] = y; coordStrides[1] = 1;

  std::vector<const scalar_t *>weightVec(1);
  std::vector<int> weightStrides(1);

  weightVec[0] = weights; weightStrides[0] = 1;

  inputAdapter_t *ia2=new inputAdapter_t(localCount, globalIds, coordVec, 
                                         coordStrides,weightVec,weightStrides);

  // Create a Zoltan2 partitioning problem

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem2 =
           new Zoltan2::PartitioningProblem<inputAdapter_t>(ia2, &params);

  // Solve the problem

  problem2->solve();

  // create metric object for MPI builds

#ifdef HAVE_ZOLTAN2_MPI
  quality_t *metricObject2 = new quality_t(ia2, &params, //problem2->getComm()
					   MPI_COMM_WORLD,
					   &problem2->getSolution());
#else
  quality_t *metricObject2 = new quality_t(ia2, &params, problem2->getComm(),
					   &problem2->getSolution());
#endif
  // Check the solution.

  if (rank == 0) {
    metricObject2->printMetrics(std::cout);
  }

  if (rank == 0){
    scalar_t imb = metricObject2->getWeightImbalance(0);
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

  // create metric object where Teuchos communicator is specified

  quality_t *metricObject3 = new quality_t(ia3, &params, problem3->getComm(),
					   &problem3->getSolution());
  // Check the solution.

  if (rank == 0) {
    metricObject3->printMetrics(std::cout);
  }

  if (rank == 0){
    scalar_t imb = metricObject3->getWeightImbalance(0);
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

  // Solution changed!

  metricObject3 = new quality_t(ia3, &params, problem3->getComm(),
                                &problem3->getSolution());
  if (rank == 0){
    metricObject3->printMetrics(std::cout);
    scalar_t imb = metricObject3->getWeightImbalance(0);
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

  // Solution changed!

  metricObject3 = new quality_t(ia3, &params, problem3->getComm(),
                                &problem3->getSolution());
  if (rank == 0){
    metricObject3->printMetrics(std::cout);
    scalar_t imb = metricObject3->getWeightImbalance(0);
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

  metricObject1 = new quality_t(ia1, &params, //problem1->getComm(),
                                &problem1->getSolution());
  // Check the solution.

  if (rank == 0) {
    metricObject1->printMetrics(std::cout);
  }

  if (rank == 0){
    scalar_t imb = metricObject1->getObjectCountImbalance();
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

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

