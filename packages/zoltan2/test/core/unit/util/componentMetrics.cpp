// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Zoltan2_componentMetrics.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <iostream>
#include <Teuchos_RCP.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>


/////////////////////////////////////////////////////////////////////////////
// Test program for componentMetrics code.
// Usage:
//     a.out
// Karen Devine, 2017
/////////////////////////////////////////////////////////////////////////////

typedef Tpetra::Map<zlno_t, zgno_t> zmap_t;
typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t> zmatrix_t;
typedef Tpetra::CrsGraph<zlno_t, zgno_t> zgraph_t;

typedef Zoltan2::XpetraCrsMatrixAdapter<zmatrix_t> matrixAdapter_t;
typedef Zoltan2::XpetraCrsGraphAdapter<zgraph_t> graphAdapter_t;

/////////////////////////////////////////////////////////////////////////////
template <typename CC_T>
int checkResult(
  std::string &name,   // test name (including rank)
  CC_T cc,             // the perProcessorComponentMetrics object for the test
  size_t nccAnswer,    // Expected answers for the test
  size_t maxAnswer,
  size_t minAnswer,
  double avgAnswer
)
{
  int ierr = 0;

  if (cc.getNumComponents() != nccAnswer) {
    std::cout << name << "Invalid number of components "
              << cc.getNumComponents() << " should be " << nccAnswer
              << std::endl;
    ierr++;
  }

  if (cc.getMaxComponentSize() != maxAnswer) {
    std::cout << name << "Maximum component size "
              << cc.getMaxComponentSize() << " should be " << maxAnswer
              << std::endl;
    ierr++;
  }

  if (cc.getMinComponentSize() != minAnswer) {
    std::cout << name << "Minimum component size "
              << cc.getMinComponentSize() << " should be " << minAnswer
              << std::endl;
    ierr++;
  }

  if (cc.getAvgComponentSize() != avgAnswer) {
    std::cout << name << "Average component size "
              << cc.getAvgComponentSize() << " should be " << avgAnswer
              << std::endl;
    ierr++;
  }

  return ierr;
}

/////////////////////////////////////////////////////////////////////////////
int test_every_third(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Testing a graph with every third vertex in the same component
  // Test using a GraphAdapter.

  std::ostringstream namestream;
  namestream << comm->getRank() << " every_third ";
  std::string name = namestream.str();

  int ierr = 0;

  // Create a default map
  const size_t gNvtx = 27;

  Teuchos::RCP<const zmap_t> map = rcp(new zmap_t(gNvtx, 0, comm));
  size_t nVtx = map->getLocalNumElements();

  // Create a Tpetra::Matrix with every third local row in the same component
  size_t maxRowLen = 3;
  Teuchos::RCP<zmatrix_t> matrix = rcp(new zmatrix_t(map, maxRowLen));

  Teuchos::Array<zgno_t> col(3);
  Teuchos::Array<zscalar_t> val(3); val[0] = 1.; val[1] = 1.; val[2] = 1.;

  for (size_t i = 0; i < nVtx; i++) {
    zgno_t id = map->getGlobalElement(i);
    col[0] = (id+3)%gNvtx;
    col[1] = (id+6)%gNvtx;
    col[2] = (id+9)%gNvtx;
    matrix->insertGlobalValues(id, col(), val());
  }

  matrix->fillComplete(map, map);

  // Create a Zoltan2 XpetraGraphAdapter
  graphAdapter_t ia(matrix->getCrsGraph(), 0);

  // Call connected components utility
  typedef Zoltan2::perProcessorComponentMetrics<graphAdapter_t> cc_t;
  cc_t cc(ia, *comm);

  // Check result:
  size_t nccAnswer = std::min<size_t>(3, nVtx);
  size_t maxAnswer = nVtx/3 + ((nVtx%3)!=0);
  size_t minAnswer = (nVtx ? std::max<size_t>(nVtx/3, 1) : 0);
  double avgAnswer = (nVtx ? double(nVtx) / double(std::min<size_t>(nVtx,3))
                           : 0.);

  ierr = checkResult<cc_t>(name, cc,
                           nccAnswer, maxAnswer, minAnswer, avgAnswer);

  return ierr;
}

/////////////////////////////////////////////////////////////////////////////
int test_dist_component(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Testing a graph with single component, distributed among procs
  // Test using a GraphAdapter.

  std::ostringstream namestream;
  namestream << comm->getRank() << " dist_component ";
  std::string name = namestream.str();

  int ierr = 0;

  // Create a default map
  const size_t gNvtx = 25;

  Teuchos::RCP<const zmap_t> map = rcp(new zmap_t(gNvtx, 0, comm));
  size_t nVtx = map->getLocalNumElements();

  // Create a Tpetra::Matrix with a single component
  size_t maxRowLen = 3;
  Teuchos::RCP<zmatrix_t> matrix = rcp(new zmatrix_t(map, maxRowLen));

  Teuchos::Array<zgno_t> col(3);
  Teuchos::Array<zscalar_t> val(3); val[0] = 1.; val[1] = 1.; val[2] = 1.;

  for (size_t i = 0; i < nVtx; i++) {
    zgno_t id = map->getGlobalElement(i);
    col[0] = (id+4)%gNvtx;
    col[1] = (id+1)%gNvtx;
    col[2] = (id+3)%gNvtx;
    matrix->insertGlobalValues(id, col(), val());
  }

  matrix->fillComplete(map, map);

  // Create a Zoltan2 XpetraGraphAdapter
  graphAdapter_t ia(matrix->getCrsGraph(), 0);

  // Call connected components utility
  typedef Zoltan2::perProcessorComponentMetrics<graphAdapter_t> cc_t;
  cc_t cc(ia, *comm);

  // Check result:
  // one component with all local vertices
  size_t nccAnswer = size_t(nVtx > 0);
  size_t maxAnswer = nVtx;
  size_t minAnswer = nVtx;
  double avgAnswer = nVtx;

  ierr = checkResult<cc_t>(name, cc,
                           nccAnswer, maxAnswer, minAnswer, avgAnswer);

  return ierr;
}

/////////////////////////////////////////////////////////////////////////////
int test_one_proc(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Testing a graph with all vertices and edges on a single proc
  // Test using a MatrixAdapter.

  std::ostringstream namestream;
  namestream << comm->getRank() << " one_proc ";
  std::string name = namestream.str();

  int ierr = 0;

  int me = comm->getRank();
  int np = comm->getSize();
  bool allOnThisProc = (me == np-1);

  // Create a map with all data on a single proc
  const size_t gNvtx = 25;
  size_t nVtx;

  if (allOnThisProc) nVtx = gNvtx;
  else               nVtx = 0;

  Teuchos::RCP<const zmap_t> map = rcp(new zmap_t(gNvtx, nVtx, 0, comm));

  // Create a Tpetra::Matrix with one component
  size_t maxRowLen = 2;
  Teuchos::RCP<zmatrix_t> matrix = rcp(new zmatrix_t(map, maxRowLen));
  if (allOnThisProc) {
    Teuchos::Array<zgno_t> col(2);
    Teuchos::Array<zscalar_t> val(2); val[0] = 1.; val[1] = 1.;
    for (size_t i = 0; i < nVtx; i++) {
      zgno_t id = map->getGlobalElement(i);
      col[0] = id;
      col[1] = (id+1)%gNvtx;
      matrix->insertGlobalValues(id, col(), val());
    }
  }
  matrix->fillComplete(map, map);

  // Create a Zoltan2 XpetraMatrixAdapter
  matrixAdapter_t ia(matrix, 0);

  // Call connected components utility
  typedef Zoltan2::perProcessorComponentMetrics<matrixAdapter_t> cc_t;
  cc_t cc(ia, *comm);

  // Check result:
  if (allOnThisProc) {
    // one component with all global vertices
    size_t nccAnswer = 1;
    size_t maxAnswer = gNvtx;
    size_t minAnswer = gNvtx;
    double avgAnswer = double(gNvtx);

    ierr = checkResult<cc_t>(name, cc,
                             nccAnswer, maxAnswer, minAnswer, avgAnswer);
  }
  else {
    // one component with all global vertices
    size_t nccAnswer = 0;
    size_t maxAnswer = 0;
    size_t minAnswer = 0;
    double avgAnswer = 0.;

    ierr = checkResult<cc_t>(name, cc,
                             nccAnswer, maxAnswer, minAnswer, avgAnswer);
  }

  return ierr;
}

/////////////////////////////////////////////////////////////////////////////
int test_no_graph(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Testing a graph with no edges
  // Test using a MatrixAdapter.
  // Test using an MPI communicator rather than Teuchos::Comm, if appropriate.

  std::ostringstream namestream;
  namestream << comm->getRank() << " no_graph ";
  std::string name = namestream.str();

  int ierr = 0;

  // Create a default map
  const size_t gNvtx = 25;

  Teuchos::RCP<const zmap_t> map = rcp(new zmap_t(gNvtx, 0, comm));
  size_t nVtx = map->getLocalNumElements();

  // Create a Tpetra::Matrix with no edges
  size_t maxRowLen = 1;
  Teuchos::RCP<zmatrix_t> matrix = rcp(new zmatrix_t(map, maxRowLen));
  matrix->fillComplete(map, map);

  // Create a Zoltan2 XpetraMatrixAdapter
  matrixAdapter_t ia(matrix, 0);

#ifdef HAVE_ZOLTAN2_MPI
  // For testing only, extract MPI communicator from Teuchos::Comm
  // There's no need to do this in real life; use Teuchos::Comm if you have it.
  // We just want to exercise the MPI_Comm interface here.
  if (comm->getRank() == 0) std::cout << "  using MPI_Comm " << std::endl;
  MPI_Comm useThisComm =  Teuchos::getRawMpiComm(*comm);
#else
  const Teuchos::Comm<int> &useThisComm = *comm;
#endif

  // Call connected components utility
  typedef Zoltan2::perProcessorComponentMetrics<matrixAdapter_t> cc_t;
  cc_t cc(ia, useThisComm);

  // Check result:
  // With no edges, every vertex should be a component
  size_t nccAnswer = nVtx;
  size_t maxAnswer = size_t(nVtx > 0);
  size_t minAnswer = size_t(nVtx > 0);
  double avgAnswer = double(nVtx > 0);

  ierr = checkResult<cc_t>(name, cc,
                           nccAnswer, maxAnswer, minAnswer, avgAnswer);

  return ierr;
}

/////////////////////////////////////////////////////////////////////////////
int main(int narg, char** arg)
{
  // Establish session.
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  int testReturn = 0;
  if (me == 0) std::cout << "test_one_proc..." << std::endl;
  testReturn += test_one_proc(comm);         // all data on a single proc
  if (me == 0) std::cout << "test_no_graph..." << std::endl;
  testReturn += test_no_graph(comm);         // no edges in graph
  if (me == 0) std::cout << "test_dist_component..." << std::endl;
  testReturn += test_dist_component(comm);   // one component per rank
  if (me == 0) std::cout << "test_every_third..." << std::endl;
  testReturn += test_every_third(comm);      // every 3rd vtx in same component

  int gtestReturn = 0;
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MAX, 1,
                               &testReturn, &gtestReturn);
  if (me == 0) {
    if (gtestReturn) std::cout << "FAIL" << std::endl;
    else             std::cout << "PASS" << std::endl;
  }

  return gtestReturn;
}
