// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace {
using std::endl;

using Teuchos::arcp;
using Teuchos::arcpClone;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::arrayView;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::broadcast;
using Teuchos::Comm;
using Teuchos::FancyOStream;
using Teuchos::null;
using Teuchos::OrdinalTraits;
using Teuchos::outArg;
using Teuchos::ParameterList;
using Teuchos::parameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ScalarTraits;
using Teuchos::tuple;

using Tpetra::createCrsMatrix;
using Tpetra::createUniformContigMapWithNode;
using Tpetra::createVector;
using Tpetra::CrsGraph;
using Tpetra::CrsMatrix;
using Tpetra::DoNotOptimizeStorage;
using Tpetra::DoOptimizeStorage;
using Tpetra::global_size_t;
using Tpetra::GloballyDistributed;
using Tpetra::Import;
using Tpetra::INSERT;
using Tpetra::Map;
using Tpetra::OptimizeOption;
using Tpetra::RowMatrix;
using Tpetra::Vector;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
}

//
// UNIT TESTS
//

////
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(CrsMatrix, NodeConversion, N2) {
}

//
// INSTANTIATIONS
//

#define NC_TESTS(N2) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(CrsMatrix, NodeConversion, N2)

TPETRA_ETI_MANGLING_TYPEDEFS()
}  // namespace
