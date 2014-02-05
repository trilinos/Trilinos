// Testing utilties
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_UnitTestHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Tpetra
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix, MatVec, Scalar, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<Node> node = rcp(new Node);
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal>(
      nrow, comm, node);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2), Tpetra::StaticProfile));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
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
  ArrayRCP<Scalar> y_view = y->get1dViewNonConst();
  ArrayRCP<Scalar> x_view = y->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    Scalar val = 2.0;

    // if (row != nrow-1) {
    //   val += 2.0;
    // }
    TEST_EQUALITY( y_view[i], val );
  }
}


//typedef KokkosClassic::DefaultNode::DefaultNodeType Node;
typedef Kokkos::Compat::KokkosThreadsWrapperNode Node;
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(
  Tpetra_CrsMatrix, MatVec, double, int, int, Node )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  return ret;
}
