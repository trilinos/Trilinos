// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace {  // (anonymous)

template <typename Scalar, typename LO, typename GO, typename Node, int Tag>
void getLocalDiagCopyTest(Teuchos::FancyOStream& out, bool& success) {
  using Teuchos::RCP;

  using map_type        = Tpetra::Map<LO, GO, Node>;
  using vec_type        = Tpetra::Vector<Scalar, LO, GO, Node>;
  using crs_matrix_type = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;

  using STS           = Teuchos::ScalarTraits<Scalar>;
  using MT            = typename STS::magnitudeType;
  const Scalar SC_ONE = STS::one();

  using LOT           = Teuchos::OrdinalTraits<LO>;
  const LO LO_INVALID = LOT::invalid();
  const LO LO_ONE     = LOT::one();
  const GO GO_ONE     = Teuchos::OrdinalTraits<GO>::one();

  int lclSuccess = 0;
  int gblSuccess = 0;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const size_t numProc                = comm->getSize();
  const size_t myProc                 = comm->getRank();

  // create a Map
  RCP<const map_type> map = Tpetra::createContigMapWithNode<LO, GO, Node>(LO_INVALID,
                                                                          LO_ONE + LO_ONE,
                                                                          comm);

  // Create a matrix with at most 3 entries per row
  RCP<crs_matrix_type> matrix = Teuchos::rcp(new crs_matrix_type(map, 3));
  const Scalar rankAsScalar   = static_cast<Scalar>(static_cast<MT>(comm->getRank()));

  Teuchos::Array<Scalar> vals = {{SC_ONE, rankAsScalar + SC_ONE, SC_ONE}};
  for (size_t lclRowIdx = 0; lclRowIdx < 2; ++lclRowIdx) {
    const GO gblRowIdx      = Teuchos::as<GO>(2 * myProc + lclRowIdx);
    Teuchos::Array<GO> cols = {{gblRowIdx - GO_ONE, gblRowIdx, gblRowIdx + GO_ONE}};

    if ((myProc == 0) && (lclRowIdx == 0)) {  // First row of the matrix
      matrix->insertGlobalValues(gblRowIdx, cols(1, 2), vals(1, 2));
    } else if ((myProc == numProc - 1) && (lclRowIdx == 1)) {  // Last row of the matrix
      matrix->insertGlobalValues(gblRowIdx, cols(0, 2), vals(0, 2));
    } else {
      matrix->insertGlobalValues(gblRowIdx, cols(), vals());
    }
  }

  matrix->fillComplete();

  // Make sure that all processes got this far.
  {
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MIN, lclSuccess, Teuchos::outArg(gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }

  RCP<vec_type> diag = Teuchos::rcp(new vec_type(map));
  diag->putScalar(-1);

  if constexpr (Tag == 0) {
    matrix->resumeFill();
  }
  matrix->getLocalDiagCopy(*diag);
}

// Unit test of getLocalDiagCopy
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, getLocalDiagCopy, Scalar, LO, GO, Node) {
  getLocalDiagCopyTest<Scalar, LO, GO, Node, 1>(out, success);
}

// Unit test of getLocalDiagCopy
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, getLocalDiagCopyFillActive, Scalar, LO, GO, Node) {
  getLocalDiagCopyTest<Scalar, LO, GO, Node, 0>(out, success);
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP(SCALAR, LO, GO, NODE)                                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, getLocalDiagCopy, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, getLocalDiagCopyFillActive, SCALAR, LO, GO, NODE)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN(UNIT_TEST_GROUP)

}  // namespace
