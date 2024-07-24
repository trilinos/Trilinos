// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test Tpetra::CrsMatrix::get{Global,Local}RowCopy for the case where
// the matrix has never yet been made fill
// complete.

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace {

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, GetGlobalRowCopy, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::tuple;
    using std::endl;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;
    typedef typename crs_matrix_type::nonconst_global_inds_host_view_type g_indices_type;
    typedef typename crs_matrix_type::nonconst_values_host_view_type values_type;

    RCP<const Comm<int> > comm = getDefaultComm ();
    const int numProcs = comm->getSize ();

    out << "Test Tpetra::CrsMatrix::getGlobalRowCopy test" << endl;
    Teuchos::OSTab tab1 (out);

    out << "Create a row Map" << endl;

    const LO lclNumRows = 10;
    const GO gblNumRows = lclNumRows * numProcs;
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp (new map_type (static_cast<GST> (gblNumRows),
                         static_cast<size_t> (lclNumRows),
                         indexBase, comm));
    // Leave room for three locally owned entries per row.
    // We will only fill in two entries per row.
    crs_matrix_type A (rowMap, 3);
    const GO gblNumCols = gblNumRows;

    out << "Fill the matrix" << endl;

    Teuchos::Array<GO> gblColInds (2);
    Teuchos::Array<Scalar> vals (2);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      gblColInds[0] = (gblRow + 1) % gblNumCols;
      gblColInds[1] = (gblRow + 2) % gblNumCols;
      // Casting to MT first ensures that reasonable Scalar types will
      // be able to convert from an integer type.  std::complex<T>
      // might not have a constructor that takes int, for example.
      vals[0] = static_cast<Scalar> (static_cast<MT> (gblColInds[0]));
      vals[1] = static_cast<Scalar> (static_cast<MT> (gblColInds[1]));

      TEST_NOTHROW( A.insertGlobalValues (gblRow, gblColInds (), vals ()) );
    }

    TEST_ASSERT( ! A.isFillComplete () );

    out << "Test getGlobalRowCopy" << endl;

    // Make the arrays bigger than necessary, just to make sure that
    // the methods behave correctly.
    g_indices_type curGblColInds ("indices",5);
    values_type curVals ("values",5);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);

      size_t numEnt = 0;
      TEST_NOTHROW( A.getGlobalRowCopy (gblRow, curGblColInds, curVals, numEnt) );
      TEST_EQUALITY( numEnt, static_cast<size_t> (2) );
      if (numEnt != static_cast<size_t> (2)) {
        break; // avoid segfault on error
      }

      TEST_EQUALITY( curGblColInds[0], (gblRow + 1) % gblNumCols );
      TEST_EQUALITY( curGblColInds[1], (gblRow + 2) % gblNumCols );

      // Casting to MT first ensures that reasonable Scalar types will
      // be able to convert from an integer type.  std::complex<T>
      // might not have a constructor that takes int, for example.
      TEST_EQUALITY( curVals[0], static_cast<Scalar> (static_cast<MT> (curGblColInds[0])) );
      TEST_EQUALITY( curVals[1], static_cast<Scalar> (static_cast<MT> (curGblColInds[1])) );
    }

    // Test whether all processes passed the test.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (gblSuccess != 1) {
      out << "  Test FAILED on some process" << endl;
    } else {
      out << "  Test succeeded on all processes" << endl;
    }
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, GetLocalRowCopy, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::tuple;
    using std::endl;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;
    typedef typename crs_matrix_type::nonconst_local_inds_host_view_type l_indices_type;
    typedef typename crs_matrix_type::nonconst_values_host_view_type values_type;

    RCP<const Comm<int> > comm = getDefaultComm ();
    const int numProcs = comm->getSize ();

    out << "Tpetra::CrsMatrix getLocalRowCopy test" << endl;
    Teuchos::OSTab tab1 (out);

    out << "Create a row Map" << endl;

    const LO lclNumRows = 10;
    const GO gblNumRows = lclNumRows * numProcs;
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp (new map_type (static_cast<GST> (gblNumRows),
                         static_cast<size_t> (lclNumRows),
                         indexBase, comm));

    // We need a locally indexed matrix to test getLocalRowCopy, so we
    // need a column Map.  It can be the same as the row Map in this
    // case, since we're not testing global effects.
    RCP<const map_type> colMap = rowMap;
    const GO lclNumCols = lclNumRows;

    // Leave room for three locally owned entries per row.
    // We will only fill in two entries per row.
    crs_matrix_type A (rowMap, colMap, 3);

    out << "Fill the matrix" << endl;

    Teuchos::Array<LO> lclColInds (2);
    Teuchos::Array<Scalar> vals (2);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      lclColInds[0] = (lclRow + 1) % lclNumCols;
      lclColInds[1] = (lclRow + 2) % lclNumCols;
      // Casting to MT first ensures that reasonable Scalar types will
      // be able to convert from an integer type.  std::complex<T>
      // might not have a constructor that takes int, for example.
      vals[0] = static_cast<Scalar> (static_cast<MT> (lclColInds[0]));
      vals[1] = static_cast<Scalar> (static_cast<MT> (lclColInds[1]));

      TEST_NOTHROW( A.insertLocalValues (lclRow, lclColInds (), vals ()) );
    }

    TEST_ASSERT( ! A.isFillComplete () );

    out << "Test getLocalRowCopy" << endl;

    // Make the arrays bigger than necessary, just to make sure that
    // the methods behave correctly.
    l_indices_type curLclColInds ("indices",5);
    values_type curVals ("values",5);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      size_t numEnt = 0;
      TEST_NOTHROW( A.getLocalRowCopy (lclRow, curLclColInds, curVals, numEnt) );
      TEST_EQUALITY( numEnt, static_cast<size_t> (2) );
      if (numEnt != static_cast<size_t> (2)) {
        break; // avoid segfault on error
      }

      TEST_EQUALITY( curLclColInds[0], (lclRow + 1) % lclNumCols );
      TEST_EQUALITY( curLclColInds[1], (lclRow + 2) % lclNumCols );

      // Casting to MT first ensures that reasonable Scalar types will
      // be able to convert from an integer type.  std::complex<T>
      // might not have a constructor that takes int, for example.
      TEST_EQUALITY( curVals[0], static_cast<Scalar> (static_cast<MT> (curLclColInds[0])) );
      TEST_EQUALITY( curVals[1], static_cast<Scalar> (static_cast<MT> (curLclColInds[1])) );
    }

    // Test whether all processes passed the test.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (gblSuccess != 1) {
      out << "  Test FAILED on some process" << endl;
    } else {
      out << "  Test succeeded on all processes" << endl;
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, GetGlobalRowCopy, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, GetLocalRowCopy, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}


