/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

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
