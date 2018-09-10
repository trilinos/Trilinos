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

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "TpetraNew_Map.hpp"
#include "TpetraNew_MultiVector.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_Details_gathervPrint.hpp"

namespace {

  // Macro used inside the test(s) below.  It tests for global error,
  // and if so, prints each process' error message and quits the test
  // early.
  //
  // 'out' only prints on Process 0.  It's really not OK for other
  // processes to print to stdout, but it usually works and we need to
  // do it for debugging.
#define TPETRA_MV_TEST_REPORT_GLOBAL_ERR( WHAT_STRING ) do {       \
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess)); \
    TEST_EQUALITY_CONST( gblSuccess, 1 );                               \
    if (gblSuccess != 1) {                                              \
      out << WHAT_STRING << " FAILED on one or more processes!" << endl; \
      Tpetra::Details::gathervPrint (out, errStrm.str () + "\n", *comm); \
      return;                                                           \
    }                                                                   \
  } while (false)

  //
  // UNIT TESTS
  //

  using GO = TpetraNew::Map::global_ordinal_type;

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiVector_new, subViewSomeZeroRows, S)
  {
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    using map_type = TpetraNew::Map;
    using MV = TpetraNew::MultiVector<S>;

    out << "Tpetra::MultiVector: Test correct dimensions of result of "
      "subView, when some processes have zero rows" << endl;

    int lclSuccess = 1;
    int gblSuccess = 1;
    std::ostringstream errStrm; // for error collection

    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    if (numProcs != 2) {
      out << "Test requires exactly two processes to run." << endl
          << "End Result: TEST FAILED" << endl;
      success = false;
      return;
    }
    //
    // 1 row on both processors
    //
    Teuchos::Array<GO> globalIDs (1, myRank);
    Teuchos::Array<size_t> columnVec (1, 0);
    const GO IGO = Teuchos::OrdinalTraits<GO>::invalid ();
    RCP<const map_type> rowMap1;

    try {
      rowMap1 = rcp (new map_type (IGO, globalIDs (), 0, comm));
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": rowMap1 constructor threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "rowMap1 constructor threw exception" );

    RCP<MV> multiVectorA;
    RCP<const MV> multiVectorB;

    try {
      multiVectorA = rcp (new MV (rowMap1, 1));
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": MV constructor threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "MV constructor threw exception" );

    try {
      multiVectorB = multiVectorA->subView (columnVec ());
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": subView threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "subView threw exception" );

    const size_t numVectorsA = multiVectorA->getNumVectors ();
    const size_t numVectorsB = multiVectorB->getNumVectors ();
    lclSuccess = (numVectorsA != numVectorsB) ? 0 : lclSuccess;
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "numVectorsA != numVectorsB" );

    //
    // 0 rows on proc 0, 2 rows on proc 1
    //
    if (myRank == 0) {
      globalIDs.resize (0);
    }
    else {
      globalIDs.resize (2);
      globalIDs[0] = 0;
      globalIDs[1] = 1;
    }
    RCP<const map_type> rowMap2 =
      rcp (new map_type (IGO, globalIDs (), 0, comm));
    RCP<MV> multiVectorC;
    RCP<const MV> multiVectorD;

    try {
      multiVectorC = rcp (new MV (rowMap2, 1));
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": MV constructor threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "MV constructor threw exception" );

    try {
      multiVectorD = multiVectorC->subView (columnVec ());
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": subView threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "subView threw exception" );

    const size_t numVectorsC = multiVectorC->getNumVectors ();
    const size_t numVectorsD = multiVectorD->getNumVectors ();
    lclSuccess = (numVectorsC != numVectorsD) ? 0 : lclSuccess;
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "numVectorsC != numVectorsD" );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiVector_new, subViewNonConstSomeZeroRows, S)
  {
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    using map_type = TpetraNew::Map;
    using MV = TpetraNew::MultiVector<S>;

    out << "Tpetra::MultiVector: Test correct dimensions of result of "
      "subViewNonConst, when some processes have zero rows" << endl;

    int lclSuccess = 1;
    int gblSuccess = 1;
    std::ostringstream errStrm; // for error collection

    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    if (numProcs != 2) {
      out << "Test requires exactly two processes to run." << endl
          << "End Result: TEST FAILED" << endl;
      success = false;
      return;
    }
    //
    // 1 row on both processors
    //
    Teuchos::Array<GO> globalIDs (1, myRank);
    Teuchos::Array<size_t> columnVec (1, 0);
    const GO IGO = Teuchos::OrdinalTraits<GO>::invalid ();
    RCP<const map_type> rowMap1;

    try {
      rowMap1 = rcp (new map_type (IGO, globalIDs (), 0, comm));
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": rowMap1 constructor threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "rowMap1 constructor threw exception" );

    RCP<MV> multiVectorA;
    RCP<MV> multiVectorB;

    try {
      multiVectorA = rcp (new MV (rowMap1, 1));
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": MV constructor threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "MV constructor threw exception" );

    try {
      multiVectorB = multiVectorA->subViewNonConst (columnVec ());
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": subViewNonConst threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "subViewNonConst threw exception" );

    const size_t numVectorsA = multiVectorA->getNumVectors ();
    const size_t numVectorsB = multiVectorB->getNumVectors ();
    lclSuccess = (numVectorsA != numVectorsB) ? 0 : lclSuccess;
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "numVectorsA != numVectorsB" );

    //
    // 0 rows on proc 0, 2 rows on proc 1
    //
    if (myRank == 0) {
      globalIDs.resize (0);
    }
    else {
      globalIDs.resize (2);
      globalIDs[0] = 0;
      globalIDs[1] = 1;
    }
    RCP<const map_type> rowMap2 =
      rcp (new map_type (IGO, globalIDs (), 0, comm));
    RCP<MV> multiVectorC;
    RCP<MV> multiVectorD;

    try {
      multiVectorC = rcp (new MV (rowMap2, 1));
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": MV constructor threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "MV constructor threw exception" );

    try {
      multiVectorD = multiVectorC->subViewNonConst (columnVec ());
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": subViewNonConst threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "subViewNonConst threw exception" );

    const size_t numVectorsC = multiVectorC->getNumVectors ();
    const size_t numVectorsD = multiVectorD->getNumVectors ();
    lclSuccess = (numVectorsC != numVectorsD) ? 0 : lclSuccess;
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "numVectorsC != numVectorsD" );
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiVector_new, subViewSomeZeroRows, SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiVector_new, subViewNonConstSomeZeroRows, SCALAR )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_S( UNIT_TEST_GROUP )

}

