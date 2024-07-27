// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Details_getNumDiags.hpp"

namespace {

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, MultipleFillCompletes, LO, GO, Scalar, Node )
  {
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::tuple;
    using std::cerr;
    using std::endl;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> ST;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numImages = comm->getSize ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    if (myRank == 0) {
      cerr << "Tpetra,CrsMatrix,MultipleFillCompletes test" << endl;
    }

    //
    // test that an exception is thrown when we exceed statically allocated memory
    //

    // create a Map
    const size_t numLocal = 1; // change to 10
    RCP<const map_type> map =
      Tpetra::createContigMapWithNode<LO, GO, Node> (INVALID, numLocal, comm);
    RCP<ParameterList> params = parameterList ();
    {
      if (myRank == 0) {
        cerr << "  Create a matrix with room for 2 entries in each row" << endl;
      }
      MAT matrix (map, 2); // room for two on each row

      if (myRank == 0) {
        cerr << "  Test insertGlobalValues does not throw "
          "if not out of room" << endl;
      }
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        TEST_NOTHROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())) );
        TEST_NOTHROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())) );
      }

      if (myRank == 0) {
        cerr << "  Test that the matrix is not yet fill complete room" << endl;
      }
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      if (myRank == 0) {
        cerr << "  Test that fillComplete (with \"Optimize Storage\" false) "
          "does not throw" << endl;
      }
      params->set ("Optimize Storage", false);
      TEST_NOTHROW( matrix.fillComplete (params) );

      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      TEST_EQUALITY_CONST( matrix.isStorageOptimized(), false );

      if (myRank == 0) {
        cerr << "  Call resumeFill and test that the matrix has room "
          "for more entries in each row" << endl;
      }
      matrix.resumeFill (); // Now there is room for more entries

      if (myRank == 0) {
        cerr << "  Test that insertLocalValues does not throw" << endl;
      }
      for (LO r = 0; r < static_cast<LO> (numLocal); ++r) {
        TEST_NOTHROW( matrix.insertLocalValues (r, tuple (r), tuple (ST::one ())) );
      }

      if (myRank == 0) {
        cerr << "  Call fillComplete with \"Optimize Storage\" true" << endl;
      }
      params->set ("Optimize Storage", true);
      TEST_NOTHROW( matrix.fillComplete(params) );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      TEST_EQUALITY_CONST( matrix.isStorageOptimized(), true );

      // test that the matrix is 3*I
      if (myRank == 0) {
        cerr << "  Test that the matrix is 3*I" << endl;
      }
      TEST_EQUALITY( Tpetra::Details::getGlobalNumDiags (matrix), static_cast<GO> (numLocal*numImages) );
      TEST_EQUALITY( Tpetra::Details::getLocalNumDiags (matrix), static_cast<LO> (numLocal) );
      TEST_EQUALITY( matrix.getGlobalNumEntries(), numLocal*numImages );
      TEST_EQUALITY( matrix.getLocalNumEntries(), numLocal );
      for (LO r = 0; r < static_cast<LO> (numLocal); ++r) {
        typename MAT::local_inds_host_view_type inds;
        typename MAT::values_host_view_type vals;
        TEST_NOTHROW( matrix.getLocalRowView(r,inds,vals) );
        TEST_COMPARE_ARRAYS( inds, tuple<LO> (r) );
        TEST_COMPARE_ARRAYS( vals, tuple<Scalar> (static_cast<Scalar> (3.0)) );

      }
    }

    // Test whether all processes passed the test.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (myRank == 0) {
      if (gblSuccess != 1) {
        cerr << "  Test FAILED on some process" << endl;
      } else {
        cerr << "  Test succeeded on all processes" << endl;
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, MultipleFillCompletes, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
