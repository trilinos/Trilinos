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

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Core.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>

namespace { // (anonymous)

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
      for (int p = 0; p < numProcs; ++p) {                              \
        if (myRank == p && lclSuccess != 1) {                           \
          std::cout << errStrm.str () << std::flush;                    \
        }                                                               \
        comm->barrier ();                                               \
        comm->barrier ();                                               \
        comm->barrier ();                                               \
      }                                                                 \
      return;                                                           \
    }                                                                   \
  } while (false)

  //
  // UNIT TESTS
  //

  // Test for Bug 5474: getData and getDataNonConst should return
  // Teuchos::null if the vector is locally empty.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Vector, Bug5474_1, S, LO, GO, NODE)
  {
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    typedef Tpetra::Map<LO, GO, NODE> map_type;
    typedef Tpetra::Vector<S, LO, GO, NODE> vec_type;
    typedef Tpetra::global_size_t GST;

    out << "Tpetra::Vector: Bug 5474 test (getData and getDataNonConst "
      "should return Teuchos::null if the vector is locally empty)" << endl;

    int lclSuccess = 1;
    int gblSuccess = 1;
    std::ostringstream errStrm; // for error collection

    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    //
    // Create a Map which has zero rows on all processes
    //
    RCP<const map_type> rowMap;
    const GST gblNumRows = 0;
    const GO indexBase = 0;

    try {
      rowMap = rcp (new map_type (gblNumRows, indexBase, comm));
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": Map constructor threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "Map constructor threw exception" );
    TEST_EQUALITY( rowMap->getGlobalNumElements (), static_cast<GST> (0) );
    TEST_EQUALITY( rowMap->getNodeNumElements (), static_cast<size_t> (0) );

    vec_type x (rowMap);
    TEST_EQUALITY( x.getNumVectors (), static_cast<size_t> (1) );

    typedef typename Teuchos::ArrayRCP<S>::size_type size_type;
    Teuchos::ArrayRCP<const S> x_const;
    try {
      x_const = x.getData (0);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": Vector::getData threw exception: "
              << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "Vector::getData threw exception" );
    TEST_EQUALITY( x_const.size (), static_cast<size_type> (0) );

    Teuchos::ArrayRCP<S> x_nonconst;
    try {
      x_nonconst = x.getDataNonConst (0);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": Vector::getDataNonConst threw "
        "exception: " << e.what () << endl;
    }
    TPETRA_MV_TEST_REPORT_GLOBAL_ERR( "Vector::getDataNonConst threw exception" );
    TEST_EQUALITY( x_nonconst.size (), static_cast<size_type> (0) );
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Vector, Bug5474_1, SCALAR, LO, GO, NODE)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

