// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <type_traits>

// mfh 14 Apr 2019: As part of #4827 debugging, I found that the only
// test in Trilinos that exercises Tpetra::deep_copy with MultiVectors
// that alias each other is an Anasazi::TraceMinDavdison test.  The
// issue happens after Anasazi prints "Copying X".  I want Tpetra's
// own tests to exercise this use case, so I'm adding this test.

namespace { // (anonymous)

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, AliasedDeepCopy, SC, LO, GO, NT )
  {
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    using map_type = Tpetra::Map<LO, GO, NT>;
    using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
    using mag_type = typename MV::mag_type;
    using STS = Teuchos::ScalarTraits<SC>;

    out << "Test Tpetra::deep_copy with aliased Tpetra::MultiVectors" << endl;
    Teuchos::OSTab tab1 (out);

    out << "Create Map" << endl;
    auto comm = Tpetra::getDefaultComm ();
    const LO lclNumRows = 31;
    const GO gblNumRows =
      static_cast<GO> (comm->getSize ()) * static_cast<GO> (lclNumRows);
    const size_t numCols = 11;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

    // Create the MultiVector to view.
    out << "Create \"original\" MultiVector X" << endl;
    MV X (map, numCols);

    // Make a MultiVector, X_noncontig, which views a subset of columns of X.
    out << "Create X_noncontig, which views a noncontigous subset "
      "of the columns of X" << endl;
    const size_t numColsSubset = 4;
    Teuchos::Array<size_t> colsSubset (4);
    // Some columns are adjacent; some aren't.
    colsSubset[0] = 2;
    colsSubset[1] = 5;
    colsSubset[2] = 6;
    colsSubset[3] = 9;
    RCP<MV> X_noncontig = X.subViewNonConst (colsSubset ());

    // Fill each column of X_noncontig with a different number,
    // corresponding to its column index in X.
    out << "Fill X_noncontig" << endl;
    {
      SC curVal;
      for (size_t j = 0; j < numColsSubset; ++j) {
        curVal = static_cast<mag_type> (colsSubset[j]) * STS::one ();
        X_noncontig->getVectorNonConst (j)->putScalar (curVal);
      }
    }

    // 0th column of X_noncontig aliases 2nd column of X_sub.
    // Teuchos::Range1D expresses an inclusive range [0, 3].
    auto X_sub = X.subViewNonConst (Teuchos::Range1D (0, 3));

    out << "Attempt Tpetra::deep_copy" << endl;

    Tpetra::deep_copy (*X_sub, *X_noncontig);

    out << "Check values in target of Tpetra::deep_copy" << endl;

    // Columns 0,1,2,3 of X got overwritten by columns 2,5,6,9 of X.
    // The overlap means that the final values in column 0 of X are
    // not well defined.  We actually care about column 2 of X.

    auto X_col = X.getVector (2);
    const mag_type expectedNormInf = mag_type (6);
    const mag_type actualNormInf = X_col->normInf ();
    TEST_EQUALITY( expectedNormInf, actualNormInf );

    out << "Make sure all processes had correct results" << endl;
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));

    if (gblSuccess != 1) {
      out << "Test FAILED on some process in the communicator!" << endl;
      success = false;
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, AliasedDeepCopy, SC, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

