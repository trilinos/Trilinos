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

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"

namespace { // (anonymous)

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TransformLocalValues, Scalar, LO, GO, Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
    const Scalar ONE = static_cast<Scalar> (1.0);
    const Scalar FIVE = static_cast<Scalar> (5.0);
    const Scalar SIX = static_cast<Scalar> (6.0);
    int lclSuccess = 0;
    int gblSuccess = 0;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    Teuchos::OSTab tab0 (out);
    out << "Tpetra transformLocalValues test" << endl;
    Teuchos::OSTab tab1 (out);
    out << "Create Map and matrix" << endl;

    // We only need one row on each process.
    const size_t lclNumRows = 1;
    const Tpetra::global_size_t gblNumRows = lclNumRows * numProcs;
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

    crs_matrix_type matrix (rowMap, 1, Tpetra::StaticProfile);

    out << "Fill matrix by calling insertGlobalValues" << endl;
    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        matrix.insertGlobalValues (gblRow,
                                   Teuchos::tuple<GO> (gblRow),
                                   Teuchos::tuple<Scalar> (ONE));
      }
    }

    out << "Call fillComplete on the matrix" << endl;
    matrix.fillComplete ();
    TEST_ASSERT( matrix.isFillComplete () );

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    typedef typename crs_matrix_type::impl_scalar_type ST;
    typedef typename Node::device_type DD;
    typedef typename Kokkos::View<LO*, DD>::HostMirror::device_type HD;
    Kokkos::View<LO*, HD> lclInds ("lclInds", 1);
    Kokkos::View<ST*, HD> vals ("vals", 1);
    auto colMap = matrix.getColMap ();

    // Use transformLocalValues to change entries of the matrix.
    matrix.resumeFill ();
    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        const LO lclCol = colMap->getLocalElement (gblRow);
        lclInds(0) = lclCol;
        vals(0) = FIVE;
        const LO numChanged =
          matrix.template transformLocalValues<Kokkos::View<LO*, HD>, Kokkos::View<ST*, HD>, std::plus<ST> > (lclRow, lclInds, vals, std::plus<ST> ());
        TEST_EQUALITY( numChanged, static_cast<LO> (1) );
      }
    }
    matrix.fillComplete ();

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Make sure that the values actually got changed.
    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        const LO lclCol = colMap->getLocalElement (gblRow);

        Teuchos::ArrayView<const LO> lclIndsT;
        Teuchos::ArrayView<const Scalar> valsT;

        matrix.getLocalRowView (lclRow, lclIndsT, valsT);
        TEST_EQUALITY( lclIndsT[0], lclCol );
        TEST_EQUALITY( valsT[0], SIX );

        LO rawNumEnt = 0;
        const Scalar* rawValsT = NULL;
        const LO* rawLclIndsT = NULL;
        const LO err = matrix.getLocalRowView (lclRow, rawNumEnt, rawValsT, rawLclIndsT);
        TEST_EQUALITY( err, static_cast<LO> (0) );
        if (err == 0) {
          TEST_EQUALITY( rawNumEnt, static_cast<LO> (lclIndsT.size ()) );
          if (rawNumEnt == static_cast<LO> (lclIndsT.size ())) {
            TEST_ASSERT( rawLclIndsT != NULL );
            if (rawLclIndsT != NULL) {
              TEST_EQUALITY( rawLclIndsT[0], lclIndsT[0] );
            }
          }
          if (rawNumEnt == static_cast<LO> (valsT.size ())) {
            TEST_ASSERT( rawValsT != NULL );
            if (rawValsT != NULL) {
              TEST_EQUALITY( rawValsT[0], valsT[0] );
            }
          }
        }
      }
    }

    // Again, use transformLocalValues to change entries of the
    // matrix, but use a different binary function this time.
    matrix.resumeFill ();
    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        const LO lclCol = colMap->getLocalElement (gblRow);
        lclInds(0) = lclCol;
        vals(0) = FIVE;
        const LO numChanged =
          matrix.template transformLocalValues<Kokkos::View<LO*, HD>, Kokkos::View<ST*, HD>, std::minus<ST> > (lclRow, lclInds, vals, std::minus<ST> ());
        TEST_EQUALITY( numChanged, static_cast<LO> (1) );
      }
    }
    matrix.fillComplete ();

    // Make sure that the values actually got changed.
    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        const LO lclCol = colMap->getLocalElement (gblRow);

        Teuchos::ArrayView<const LO> lclIndsT;
        Teuchos::ArrayView<const Scalar> valsT;

        matrix.getLocalRowView (lclRow, lclIndsT, valsT);
        TEST_EQUALITY( lclIndsT[0], lclCol );
        TEST_EQUALITY( valsT[0], ONE );
      }
    }

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (myRank == 0) {
      out << "All done!  Test " << (gblSuccess ? "succeeded" : "FAILED") << "."
          << endl;
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TransformGlobalValues, Scalar, LO, GO, Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
    const Scalar ONE = static_cast<Scalar> (1.0);
    const Scalar TWO = static_cast<Scalar> (2.0);
    const Scalar FOUR = static_cast<Scalar> (4.0);
    const Scalar FIVE = static_cast<Scalar> (5.0);
    const Scalar SIX = static_cast<Scalar> (6.0);
    int lclSuccess = 0;
    int gblSuccess = 0;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    Teuchos::OSTab tab0 (out);
    out << "Tpetra transformGlobalValues test" << endl;
    Teuchos::OSTab tab1 (out);
    out << "Create Map and matrix" << endl;

    // We only need one row on each process.
    const size_t lclNumRows = 1;
    const Tpetra::global_size_t gblNumRows = lclNumRows * numProcs;
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

    crs_matrix_type matrix (rowMap, 1, Tpetra::StaticProfile);

    out << "Fill matrix by calling insertGlobalValues" << endl;
    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        matrix.insertGlobalValues (gblRow,
                                   Teuchos::tuple<GO> (gblRow),
                                   Teuchos::tuple<Scalar> (ONE));
      }
    }

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    out << "Modify matrix entries by calling transformGlobalValues "
      "(while matrix is still globally indexed)" << endl;

    typedef typename crs_matrix_type::impl_scalar_type ST;
    typedef typename Node::device_type DD;
    typedef typename Kokkos::View<LO*, DD>::HostMirror::device_type HD;
    Kokkos::View<GO*, HD> gblInds ("gblInds", 1);
    Kokkos::View<ST*, HD> vals ("vals", 1);

    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        const GO gblCol = gblRow;
        gblInds(0) = gblCol;
        vals(0) = FIVE;
        const LO numChanged =
          matrix.template transformGlobalValues<std::plus<ST>, HD> (gblRow, gblInds, vals, std::plus<ST> ());
        TEST_EQUALITY( numChanged, static_cast<LO> (1) );
      }
    }

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    out << "Before calling fillComplete, test the values" << endl;
    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        const LO gblCol = gblRow;

        Teuchos::ArrayView<const GO> gblIndsT;
        Teuchos::ArrayView<const Scalar> valsT;

        matrix.getGlobalRowView (gblRow, gblIndsT, valsT);
        TEST_EQUALITY( gblIndsT[0], gblCol );
        TEST_EQUALITY( valsT[0], SIX );
      }
    }

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    out << "Call fillComplete on the matrix" << endl;
    matrix.fillComplete ();
    TEST_ASSERT( matrix.isFillComplete () );

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    out << "Modify matrix entries by calling transformGlobalValues "
      "(while matrix is locally indexed)" << endl;

    Kokkos::View<LO*, HD> lclInds ("lclInds", 1);
    auto colMap = matrix.getColMap ();
    matrix.resumeFill ();
    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        const GO gblCol = gblRow;
        gblInds(0) = gblCol;
        vals(0) = TWO;
        const LO numChanged =
          matrix.template transformGlobalValues<std::minus<ST>, HD> (gblRow, gblInds, vals, std::minus<ST> ());
        TEST_EQUALITY( numChanged, static_cast<LO> (1) );
      }
    }
    matrix.fillComplete ();

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Make sure that the values actually got changed.
    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        const LO lclCol = colMap->getLocalElement (gblRow);

        Teuchos::ArrayView<const LO> lclIndsT;
        Teuchos::ArrayView<const Scalar> valsT;

        matrix.getLocalRowView (lclRow, lclIndsT, valsT);
        TEST_EQUALITY( lclIndsT[0], lclCol );
        TEST_EQUALITY( valsT[0], FOUR );
      }
    }

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (myRank == 0) {
      out << "All done!  Test " << (gblSuccess ? "succeeded" : "FAILED") << "."
          << endl;
    }
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, TransformLocalValues, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, TransformGlobalValues, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

} // namespace (anonymous)

