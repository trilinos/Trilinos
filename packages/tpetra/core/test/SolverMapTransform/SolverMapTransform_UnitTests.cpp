// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_SolverMap_LinearProblem.hpp>
#include "TestCase1_decl.hpp"

namespace {

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( SolverMapTransform, TestCase_1, LO, GO, Node )
  {
    using Scalar_t  = Tpetra::CrsMatrix<>::scalar_type;
    using Problem_t = Tpetra::LinearProblem<Scalar_t, LO, GO, Node>;
    auto comm = Tpetra::getDefaultComm();

    TestCase1<Scalar_t, LO, GO, Node> testCase( comm );
    Teuchos::RCP<Problem_t> linearProblem = Teuchos::rcp<Problem_t>( testCase.linearProblem(), false );

    Tpetra::SolverMap_LinearProblem<Scalar_t, LO, GO, Node> solverMapTransform;
    Teuchos::RCP<Problem_t> transformedLinearProblem = solverMapTransform( linearProblem );
    solverMapTransform.fwd();

    bool checkResult = testCase.checkTransformedProblem( transformedLinearProblem.get() );
    TEUCHOS_ASSERT(checkResult);
  }

#define UNIT_TEST_GROUP( LO, GO, NODE ) \
            TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( SolverMapTransform, TestCase_1, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

}
