// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file mj_int_coordinates.cpp
    \brief Generate a test to partition integer coordinates
    See definition of int_scalar_t
*/

#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Tpetra_Map.hpp>
#include <vector>
#include <cstdlib>

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int rank = comm->getRank(); 
  int nprocs = comm->getSize();
  int nFail = 0;

  typedef Tpetra::Map<> Map_t;
  typedef Map_t::local_ordinal_type localId_t;
  typedef Map_t::global_ordinal_type globalId_t;
  typedef int int_scalar_t;  // This is the case we are testing here.

  typedef Zoltan2::BasicUserTypes<int_scalar_t, localId_t, globalId_t> myTypes;
  typedef Zoltan2::BasicVectorAdapter<myTypes> inputAdapter_t;

  typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;

  ///////////////////////////////////////////////////////////////////////
  // Create input data.

  size_t localCount = 40;
  int dim = 3;

  // Create coordinates that range from 0 to 9
  int_scalar_t *coords = new int_scalar_t [dim * localCount];
  int_scalar_t *x = coords; 
  int_scalar_t *y = x + localCount; 
  int_scalar_t *z = y + localCount; 

  srand(rank);
  for (size_t i=0; i < localCount*dim; i++)
    coords[i] = int_scalar_t(rand() % 10);

  // Create global ids for the coordinates.
  globalId_t *globalIds = new globalId_t [localCount];
  globalId_t offset = rank * localCount;
  for (size_t i=0; i < localCount; i++) globalIds[i] = offset++;
   
  ///////////////////////////////////////////////////////////////////////
  // Create parameters for an MJ problem

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("error_check_level", "debug_mode_assertions");

  params.set("algorithm", "multijagged");
  params.set("num_global_parts", nprocs+1);

  ///////////////////////////////////////////////////////////////////////
  // Test one:  No weights

  inputAdapter_t *ia1 = new inputAdapter_t(localCount,globalIds,x,y,z,1,1,1);

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem1 =
           new Zoltan2::PartitioningProblem<inputAdapter_t>(ia1, &params);
   
  problem1->solve();

  quality_t *metricObject1 = new quality_t(ia1, &params, comm,
					   &problem1->getSolution());
  if (rank == 0){

    metricObject1->printMetrics(std::cout);

    double imb = metricObject1->getObjectCountImbalance();
    if (imb <= 1.01)  // Should get perfect balance
      std::cout << "no weights -- balance satisfied: " << imb << std::endl;
    else {
      std::cout << "no weights -- balance failure: " << imb << std::endl;
      nFail++;
    }
    std::cout << std::endl;
  }
  delete metricObject1;
  delete problem1;
  delete ia1;
   
  ///////////////////////////////////////////////////////////////////////
  // Test two:  weighted
  // Create a Zoltan2 input adapter that includes weights.

  int_scalar_t *weights = new int_scalar_t [localCount];
  for (size_t i=0; i < localCount; i++) weights[i] = 1 + int_scalar_t(rank);

  std::vector<const int_scalar_t *>coordVec(3);
  std::vector<int> coordStrides(3);

  coordVec[0] = x; coordStrides[0] = 1;
  coordVec[1] = y; coordStrides[1] = 1;
  coordVec[2] = z; coordStrides[2] = 1;

  std::vector<const int_scalar_t *>weightVec(1);
  std::vector<int> weightStrides(1);

  weightVec[0] = weights; weightStrides[0] = 1;

  inputAdapter_t *ia2=new inputAdapter_t(localCount, globalIds, coordVec, 
                                         coordStrides,weightVec,weightStrides);

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem2 =
           new Zoltan2::PartitioningProblem<inputAdapter_t>(ia2, &params);

  problem2->solve();

  quality_t *metricObject2 = new quality_t(ia2, &params, comm,
					   &problem2->getSolution());
  if (rank == 0){

    metricObject2->printMetrics(std::cout);

    double imb = metricObject2->getWeightImbalance(0);
    if (imb <= 1.01)
      std::cout << "weighted -- balance satisfied " << imb << std::endl;
    else {
      std::cout << "weighted -- balance failed " << imb << std::endl;
      nFail++;
    }
    std::cout << std::endl;
  }
  delete metricObject2;

  if (weights) delete [] weights;
  if (coords) delete [] coords;
  if (globalIds) delete [] globalIds;
  delete problem2;
  delete ia2;

  if (rank == 0) { 
    if (nFail == 0) std::cout << "PASS" << std::endl;
    else  std::cout << "FAIL:  " << nFail << " tests failed" << std::endl;
  }

  return 0;
}

