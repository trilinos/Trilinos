// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file mj_imbalanced.cpp
    \brief Generate a test to partition lots of points, all at the same coord
    Exposed an overflow in Multijagged partitioning, determining weights to
    the left of each cut.  #4785
*/

#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Tpetra_Map.hpp>
#include <vector>
#include <cstdlib>

typedef Tpetra::Map<> myMap_t;
typedef myMap_t::local_ordinal_type myLocalId_t;
typedef myMap_t::global_ordinal_type myGlobalId_t;
typedef double myScalar_t;

///////////////////////////////////////////////////////////////////////////////
template <typename A, typename B>
void copyArrays(size_t n, A *&a, B *b)
{
  a = new A[n];
  for (size_t i = 0; i < n; i++) a[i] = b[i];
}

template <typename A, typename B>
void copyArrays(size_t n, A *&a, A *b)
{
  a = b;
}

///////////////////////////////////////////////////////////////////////////////
template <typename A, typename B>
void freeArrays(A *&a, B *b)
{
  delete [] a;
}

template <typename A, typename B>
void freeArrays(A *&a, A *b)
{
  // no delete needed since only copied pointer
}

///////////////////////////////////////////////////////////////////////////////
void test_no_weights(
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
  Teuchos::ParameterList &params,
  size_t localCount, 
  myGlobalId_t *globalIds, 
  myScalar_t *coords, 
  int &nFail
)
{
  typedef Tpetra::Map<>::node_type myNode_t;

  typedef Zoltan2::BasicUserTypes<myScalar_t, myLocalId_t, myGlobalId_t,
                                  myNode_t> myTypes;
  typedef Zoltan2::BasicVectorAdapter<myTypes> inputAdapter_t;
  typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;

  myScalar_t *sx = coords; 
  myScalar_t *sy = sx + localCount;
  myScalar_t *sz = sy + localCount;

  inputAdapter_t *ia = new inputAdapter_t(localCount,globalIds,sx,sy,sz,1,1,1);

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem =
           new Zoltan2::PartitioningProblem<inputAdapter_t>(ia, &params);
   
  problem->solve();

  quality_t *metricObject = new quality_t(ia, &params, comm,
					   &problem->getSolution());
  if (comm->getRank() == 0){

    metricObject->printMetrics(std::cout);

    double imb = metricObject->getObjectCountImbalance();
    if (imb <= 1.01)  // Should get perfect balance
      std::cout << "no weights -- balance satisfied: " << imb << std::endl;
    else {
      std::cout << "no weights -- balance failure: " << imb << std::endl;
      nFail++;
    }
    std::cout << std::endl;
  }

  delete metricObject;
  delete problem;
  delete ia;
}

///////////////////////////////////////////////////////////////////////////////
void test_weights(
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
  Teuchos::ParameterList &params,
  size_t localCount, 
  myGlobalId_t *globalIds, 
  myScalar_t *coords, 
  myScalar_t *weights, 
  int &nFail
)
{
  typedef Tpetra::Map<>::node_type myNode_t;

  typedef Zoltan2::BasicUserTypes<myScalar_t, myLocalId_t, myGlobalId_t,
    myNode_t> myTypes;
  typedef Zoltan2::BasicVectorAdapter<myTypes> inputAdapter_t;
  typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;

  std::vector<const myScalar_t *> coordVec(3);
  std::vector<int> coordStrides(3);

  coordVec[0] = coords; coordStrides[0] = 1;
  coordVec[1] = coords + localCount; coordStrides[1] = 1;
  coordVec[2] = coords + localCount + localCount; coordStrides[2] = 1;

  std::vector<const myScalar_t *> weightVec(1);
  std::vector<int> weightStrides(1);

  weightVec[0] = weights; weightStrides[0] = 1;

  inputAdapter_t *ia=new inputAdapter_t(localCount, globalIds, coordVec, 
                                         coordStrides,weightVec,weightStrides);

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem =
           new Zoltan2::PartitioningProblem<inputAdapter_t>(ia, &params);

  problem->solve();

  quality_t *metricObject = new quality_t(ia, &params, comm,
                                          &problem->getSolution());
  if (comm->getRank() == 0){

    metricObject->printMetrics(std::cout);

    double imb = metricObject->getWeightImbalance(0);
    if (imb <= 1.01)
      std::cout << "weighted -- balance satisfied " << imb << std::endl;
    else {
      std::cout << "weighted -- balance failed " << imb << std::endl;
      nFail++;
    }
    std::cout << std::endl;
  }

  delete metricObject;
  delete problem;
  delete ia;
}

///////////////////////////////////////////////////////////////////////////////
int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank(); 
  int np = comm->getSize();
  int nFail = 0;

  ///////////////////////////////////////////////////////////////////////
  // Create parameters for an MJ problem

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("error_check_level", "debug_mode_assertions");
 
  params.set("algorithm", "multijagged");
  params.set("num_global_parts", 64);
  
  ///////////////////////////////////////////////////////////////////////
  // Create input data.

  const size_t N = 9000000;
  size_t maxLocalCount = N / np; // biggest test we'll run

  // Create coordinates that range from 0 to 999
  const int dim = 3;
  myScalar_t *coords = new myScalar_t[dim * maxLocalCount];

  srand(me);
  for (size_t i=0; i < maxLocalCount*dim; i++) 
    coords[i] = myScalar_t(0);

  // Create weights for the coordinates
  myScalar_t *weights = new myScalar_t[maxLocalCount];
  for (size_t i=0; i < maxLocalCount; i++) 
    weights[i] = myScalar_t(1+me);

  // Allocate space for global IDs; they will be generated below
  myGlobalId_t *globalIds = new myGlobalId_t[maxLocalCount];

  size_t localCount = N / np;

  // Create consecutive global ids for the coordinates
  myGlobalId_t offset = me * localCount;
  for (size_t i=0; i < localCount; i++)
    globalIds[i] = offset++;

  if (me == 0) {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "myGlobalId_t = " << typeid(offset).name()
              << " " << sizeof(offset)
              << "; localCount = " << localCount
              << "; globalCount = " << np * localCount << std::endl;
  }

  ///////////////////////////////////////////////////////////////////////
  // Test one:  No weights
  if (me == 0) std::cout << "Test:  no weights, scalar = double" << std::endl;
  test_no_weights(comm, params, localCount, globalIds, coords, nFail);

  ///////////////////////////////////////////////////////////////////////
  // Test two:  weighted
  if (me == 0) std::cout << "Test:  weights, scalar = double" << std::endl;
  test_weights(comm, params, localCount, globalIds, coords, weights, nFail);

  // Early exit when failure is detected
  int gnFail;
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MAX, 1, &nFail, &gnFail);

  delete [] weights;
  delete [] coords;
  delete [] globalIds;

  if (me == 0) { 
    if (nFail == 0) std::cout << "PASS" << std::endl;
    else  std::cout << "FAIL:  " << nFail << " tests failed" << std::endl;
  }

  return 0;
}

