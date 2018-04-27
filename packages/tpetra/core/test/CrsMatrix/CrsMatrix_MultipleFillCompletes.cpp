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

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_ConfigDefs.hpp>

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include "Tpetra_Details_getNumDiags.hpp"

namespace {

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, MultipleFillCompletes, LO, GO, Scalar, Node )
  {
    using Tpetra::TestingUtilities::getDefaultComm;
    using Tpetra::TestingUtilities::getNode;
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

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
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
      Tpetra::createContigMapWithNode<LO, GO, Node> (INVALID, numLocal,
                                                     comm, node);
    RCP<ParameterList> params = parameterList ();
    {
      if (myRank == 0) {
        cerr << "  Create a matrix with room for 2 entries in each row" << endl;
      }
      MAT matrix (map, 2, Tpetra::StaticProfile); // room for two on each row

      if (myRank == 0) {
        cerr << "  Test insertGlobalValues does not throw "
          "if not out of room" << endl;
      }
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        TEST_NOTHROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())) );
        TEST_NOTHROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())) );
      }

      if (myRank == 0) {
        cerr << "  Test that insertGlobalValues throws if out of room" << endl;
      }
      {
        GO r = map->getMinGlobalIndex();
        TEST_THROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())), std::runtime_error );
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
      TEST_EQUALITY( matrix.getNodeNumEntries(), numLocal );
      for (LO r = 0; r < static_cast<LO> (numLocal); ++r) {
        ArrayView<const LO> inds;
        ArrayView<const Scalar> vals;
        TEST_NOTHROW( matrix.getLocalRowView(r,inds,vals) );
        TEST_COMPARE_ARRAYS( inds, tuple<LO> (r) );
        TEST_COMPARE_ARRAYS( vals, tuple<Scalar> (static_cast<Scalar> (3.0)) );

        LO rawNumEnt = 0;
        const Scalar* rawVals = NULL;
        const LO* rawInds = NULL;
        TEST_NOTHROW( matrix.getLocalRowView (r, rawNumEnt, rawVals, rawInds) );
        TEST_EQUALITY( rawNumEnt, static_cast<LO> (1) );
        TEST_ASSERT( rawVals != NULL && rawInds != NULL );
        if (rawVals != NULL && rawInds != NULL) {
          TEST_EQUALITY( rawInds[0], r );
          TEST_EQUALITY( rawVals[0], static_cast<Scalar> (3.0) );
        }
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


