// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "SingletonFiltering_TestUtils.hpp"

#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHelpers.hpp>
#include <Tpetra_CrsSingletonFilter_LinearProblem.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <TpetraExt_MatrixMatrix.hpp>

namespace {

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Tpetra::TestingUtilities::getDefaultComm;

// Unit Tests
// --------------------------------------------------------------------------

// Singleton Filtering Test
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SingletonFiltering_P1, fwd, LO, GO, Scalar, Node) {
  auto Comm = Tpetra::getDefaultComm();

  test_Singleton_fwd<Scalar, LO, GO, Node>(
      "SF1_Matrix_Original.mm", "SF1_LHS_Original.mm", "SF1_RHS_Original.mm",
      "SF1_Matrix_Reduced.mm", "SF1_LHS_Reduced.mm", "SF1_RHS_Reduced.mm",
      Comm, out, success);
}

#define UNIT_TEST_GROUP(SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SingletonFiltering_P1, fwd, LO, GO, SCALAR, NODE)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN(UNIT_TEST_GROUP)

}  // namespace
