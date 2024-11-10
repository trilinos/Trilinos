// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file 762.cpp
/// \brief Regression test for Github Issue #762

#include "Ifpack2_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Ifpack2_UnitTestHelpers.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include <sstream>

namespace { // (anonymous)

using Tpetra::global_size_t;
using Teuchos::FancyOStream;
using Teuchos::getFancyOStream;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::cerr;
using std::endl;

template<class Scalar, class LO, class GO, class Node>
void
test762 (FancyOStream& out,
         bool& success,
         const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& A,
         const Scalar& alpha,
         const Scalar& beta,
         const size_t numVecs)
{
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  int lclSuccess = 1; // to be updated below
  int gblSuccess = 0; // output argument
  std::ostringstream errStrm;

  Teuchos::OSTab tab0 (out);
  RCP<const Teuchos::Comm<int> > comm = A.getMap ()->getComm ();

  out << "Creating test problem" << endl;
  MV x (A.getRowMap (), numVecs);
  MV y (A.getRowMap (), numVecs);
  x.putScalar (STS::one ());

  out << "Creating copies of x and y" << endl;
  MV x_copy = Tpetra::createCopy (x);
  MV y_copy = Tpetra::createCopy (y);

  out << "Testing apply() for alpha = " << alpha
      << " and beta = " << beta << endl;
  try {
    A.apply (x_copy, y_copy, Teuchos::NO_TRANS, alpha, beta);
    out << "apply (alpha = " << alpha << ", beta = " << beta << ") returned!"
        << endl;
  }
  catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << comm->getRank () << ": CrsMatrix::apply threw "
      "an exception: " << e.what () << endl;
  }
  catch (...) {
    lclSuccess = 0;
    errStrm << "Process " << comm->getRank () << ": CrsMatrix::apply threw "
      "an exception not a subclass of std::exception." << endl;
  }
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *comm);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2, Issue762, Scalar, LO, GO)
{
  typedef tif_utest::Node Node;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // Whether to test the known bad case (that triggered #762).
  const bool testKnownBad = true;

  // Whether to print to cerr, for more immediate output ('out'
  // doesn't actually print until the test returns, so it's not so
  // helpful if the test segfaults).
  const bool printToCerr = true;

  int lclSuccess = 1; // to be updated below
  int gblSuccess = 0; // output argument
  std::ostringstream errStrm; // for collecting error output on each MPI process

  RCP<FancyOStream> foutPtr;
  if (printToCerr) {
    foutPtr = getFancyOStream (rcpFromRef (cerr));
  }
  else {
    foutPtr = rcpFromRef (out);
  }
  FancyOStream& fout = *foutPtr;

  fout << "Ifpack2: Test #762" << endl;
  Teuchos::OSTab tab1 (fout);

  const global_size_t num_rows_per_proc = 5;
  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO, GO, Node> (num_rows_per_proc);
  TEST_ASSERT( ! rowmap.is_null () );
  if (rowmap.is_null ()) {
    return; // that's the best we can do
  }
  RCP<const Teuchos::Comm<int> > comm = rowmap->getComm ();
  TEST_ASSERT( ! comm.is_null () );
  if (comm.is_null ()) {
    return; // that's the best we can do
  }
  if (comm->getSize () > 1) {
    fout << "This test may only be run in serial "
      "or with a single MPI process." << endl;
    return;
  }

  fout << "Creating matrix" << endl;
  RCP<const crs_matrix_type> crsmatrix;
  try {
    crsmatrix = tif_utest::create_test_matrix2<Scalar,LO,GO,Node>(rowmap);
    if (printToCerr) {
      cerr << "create_test_matrix2 returned!" << endl;
    }
    // Don't print to 'out' here, since we only care about per-process
    // printing in this case.
  }
  catch (std::exception& e) {
    success = false;
    errStrm << "Process " << comm->getRank () << ": create_test_matrix2 threw "
      "an exception: " << e.what () << endl;
  }
  catch (...) {
    success = false;
    errStrm << "Process " << comm->getRank () << ": create_test_matrix2 threw "
      "an exception not a subclass of std::exception." << endl;
  }
  TEST_ASSERT( ! crsmatrix.is_null () );
  lclSuccess = (lclSuccess == 0 || ! success) ? 0 : 1;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *comm);
    return;
  }

  const Scalar ZERO = STS::zero ();
  const Scalar ONE = STS::one ();
  const Scalar TWO = ONE + ONE;

  // Test the following values of alpha in Y := alpha*(A*X) + beta*Y.
  const int numAlphaValues = 4;
  const Scalar alphaValues[] = { ZERO, ONE, -ONE, TWO };

  const Scalar beta = ZERO; // we only test one value of beta

  // Test the following numbers of vectors (columns) in X and Y.
  const int numNumVecsValues = 7;
  const size_t numVecsValues[] = {1, 2, 3, 4, 5, 6, 7};

  for (int whichNumVecs = 0; whichNumVecs < numNumVecsValues; ++whichNumVecs) {
    const size_t numVecs = numVecsValues[whichNumVecs];

    for (int whichAlpha = 0; whichAlpha < numAlphaValues; ++whichAlpha) {
      const Scalar alpha = alphaValues[whichAlpha];

      const bool theBadCase =
        (numVecs == 2 && alpha == TWO && beta == ZERO); // the known bad case
      const bool runTheTest =
        ! theBadCase || (theBadCase && testKnownBad);
      if (runTheTest) {
        fout << "Running case alpha = " << alpha << ", beta = " << beta
             << ", numVecs = " << numVecs << endl;
        test762<Scalar, LO, GO, Node> (fout, success, *crsmatrix,
                                       alpha, beta, numVecs);
      }
      else {
        fout << "Skipping case alpha = " << alpha << ", beta = " << beta
             << ", numVecs = " << numVecs << endl;
      }
    }
  }
}

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2, Issue762, SC, LO, GO )

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types.

IFPACK2_INSTANTIATE_SLG( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

