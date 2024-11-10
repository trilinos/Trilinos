// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsGraphTransposer.hpp>

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ScalarTraits;

  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::createUniformContigMapWithNode;
  using Tpetra::MatrixMarket::Reader;
  using Tpetra::CrsGraphTransposer;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraphTransposer, RectangularTranspose, LO, GO, Node )
  {
    typedef CrsMatrix<>::scalar_type Scalar;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef CrsGraph<LO,GO,Node> GRAPH;
    auto comm = Tpetra::getDefaultComm();

    RCP<MAT> matrix  = Reader<MAT>::readSparseFile("../RowMatrixTransposer/a.mtx", comm);
    RCP<MAT> matrixT = Reader<MAT>::readSparseFile("../RowMatrixTransposer/atrans.mtx", comm);

    RCP<const GRAPH> graph  = matrix->getCrsGraph();
    RCP<const GRAPH> graphT = matrixT->getCrsGraph();

    CrsGraphTransposer<LO, GO, Node> gt (graph);
    RCP<const GRAPH> calculatedT = gt.createTranspose();

    TEUCHOS_ASSERT(calculatedT->isIdenticalTo(*graphT));
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraphTransposer, Symmetrization, LO, GO, Node )
  {
    typedef CrsMatrix<>::scalar_type Scalar;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef CrsGraph<LO,GO,Node> GRAPH;
    auto comm = Tpetra::getDefaultComm();

    RCP<MAT> matrix    = Reader<MAT>::readSparseFile("L.mtx", comm);
    RCP<MAT> matrixT   = Reader<MAT>::readSparseFile("LT.mtx", comm);
    RCP<MAT> matrixSym = Reader<MAT>::readSparseFile("L+LT.mtx", comm);

    RCP<const GRAPH> graph    = matrix->getCrsGraph();
    RCP<const GRAPH> graphT   = matrixT->getCrsGraph();
    RCP<const GRAPH> graphSym = matrixSym->getCrsGraph();

    CrsGraphTransposer<LO, GO, Node> gt (graph);
    RCP<const GRAPH> calculatedT = gt.createTranspose();

    TEUCHOS_ASSERT(calculatedT->isIdenticalTo(*graphT));

    RCP<const GRAPH> calculatedSym = gt.symmetrize();

    Tpetra::MatrixMarket::Writer<MAT>::writeSparseGraphFile("test",*calculatedSym);

    TEUCHOS_ASSERT(calculatedSym->isIdenticalTo(*graphSym));
  }


#define UNIT_TEST_GROUP( LO, GO, NODE ) \
            TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraphTransposer, RectangularTranspose, LO, GO, NODE ) \
            TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraphTransposer, Symmetrization, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

}
