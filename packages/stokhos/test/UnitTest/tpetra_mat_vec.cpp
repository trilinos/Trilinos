// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Testing utilties
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_UnitTestHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Tpetra
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Belos solver
#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_BELOS
#include "BelosTpetraAdapter.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Tpetra_CrsMatrix, MatVec, Scalar, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using LocalOrdinal = Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = Tpetra::Map<>::global_ordinal_type;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2)));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    size_t ncol = 1;

    // if (row != nrow-1) {
    //   columnIndices[1] = row+1;
    //   ncol = 2;
    // }
    graph->insertGlobalIndices(row, columnIndices(0,ncol));
  }
  graph->fillComplete();
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph));

  // Set values in matrix
  Array<Scalar> vals(2);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    vals[0] = Scalar(2.0);
    size_t ncol = 1;

    // if (row != nrow-1) {
    //   columnIndices[1] = row+1;
    //   vals[1] = Scalar(2.0);
    //   ncol = 2;
    // }
    matrix->replaceGlobalValues(row, columnIndices(0,ncol), vals(0,ncol));
  }
  matrix->fillComplete();

  // Fill vector
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  x->putScalar(1.0);

  matrix->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
                   Teuchos::VERB_EXTREME);

  x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
              Teuchos::VERB_EXTREME);

  // Multiply
  RCP<Tpetra_Vector> y = Tpetra::createVector<Scalar>(map);
  matrix->apply(*x, *y);

  y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
              Teuchos::VERB_EXTREME);

  // Check
  auto y_view = y->getLocalViewHost(Tpetra::Access::ReadOnly);
  for (size_t i=0; i<num_my_row; ++i) {
    //const GlobalOrdinal row = myGIDs[i];
    Scalar val = 2.0;

    // if (row != nrow-1) {
    //   val += 2.0;
    // }
    TEST_EQUALITY( y_view(i,0), val );
  }
}

#if defined(HAVE_STOKHOS_BELOS)

//
// Test Belos GMRES solve for a simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Tpetra_CrsMatrix, BelosGMRES, Scalar, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using LocalOrdinal = Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = Tpetra::Map<>::global_ordinal_type;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2)));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    size_t ncol = 1;
    // if (row != nrow-1) {
    //   columnIndices[1] = row+1;
    //   ncol = 2;
    // }
    graph->insertGlobalIndices(row, columnIndices(0,ncol));
  }
  graph->fillComplete();
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph));

  // Set values in matrix
  Array<Scalar> vals(2);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    vals[0] = Scalar(2.0);
    size_t ncol = 1;

    // if (row != nrow-1) {
    //   columnIndices[1] = row+1;
    //   vals[1] = Scalar(2.0);
    //   ncol = 2;
    // }
    matrix->replaceGlobalValues(row, columnIndices(0,ncol), vals(0,ncol));
  }
  matrix->fillComplete();

  // Fill RHS vector
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  {
      auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll);
      for (size_t i=0; i<num_my_row; ++i) {
          b_view(i, 0) = Scalar(1.0);
      }
  }

  // Solve
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  typedef Belos::LinearProblem<Scalar,MV,OP> BLinProb;
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  RCP< BLinProb > problem = rcp(new BLinProb(matrix, x, b));
  RCP<ParameterList> belosParams = rcp(new ParameterList);
  typename ST::magnitudeType tol = 1e-9;
  belosParams->set("Flexible Gmres", false);
  belosParams->set("Num Blocks", 100);
  belosParams->set("Convergence Tolerance", Scalar(tol));
  belosParams->set("Maximum Iterations", 100);
  belosParams->set("Verbosity", 33);
  belosParams->set("Output Style", 1);
  belosParams->set("Output Frequency", 1);
  belosParams->set("Output Stream", out.getOStream());
  RCP<Belos::SolverManager<Scalar,MV,OP> > solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,MV,OP>(problem, belosParams));
  problem->setProblem();
  Belos::ReturnType ret = solver->solve();
  TEST_EQUALITY_CONST( ret, Belos::Converged );

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check -- Correct answer is:
  //     [ 0, 0,   ..., 0            ]
  //     [ 1, 1/2, ..., 1/VectorSize ]
  //     [ 0, 0,   ..., 0            ]
  //     [ 1, 1/2, ..., 1/VectorSize ]
  //     ....
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly);
  Scalar val = Scalar(0.5);
  for (size_t i=0; i<num_my_row; ++i) {
    // const GlobalOrdinal row = myGIDs[i];
    // if (row % 2) {
    //   val = Scalar(0.5);
    // }
    // else
    //   val = Scalar(0.0);

    // Set small values to zero
    Scalar v = x_view(i,0);
    if (ST::magnitude(v) < tol)
      v = Scalar(0.0);

    TEST_FLOATING_EQUALITY(v, val, tol);
  }
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Tpetra_CrsMatrix, BelosGMRES, Scalar, Node )
{}

#endif

using Node = Tpetra::Map<>::node_type;

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(
  Tpetra_CrsMatrix, MatVec, double, Node )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(
  Tpetra_CrsMatrix, BelosGMRES, double, Node )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  return ret;
}
