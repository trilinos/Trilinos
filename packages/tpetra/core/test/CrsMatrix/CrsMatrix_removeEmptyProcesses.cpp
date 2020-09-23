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

#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace { // (anonymous)

using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::MpiComm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using Teuchos::SerialComm;
using Teuchos::toString;
using Teuchos::tuple;
using std::cerr;
using std::endl;

typedef Tpetra::global_size_t GST;

// This test is only meaningful in an MPI build.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, RemoveEmptyProcesses, Scalar, LO, GO, Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    using std::cerr;
    using std::endl;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    const Scalar ONE = STS::one ();
    int lclSuccess = 0;
    int gblSuccess = 0;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int commSize = comm->getSize();

    if (myRank == 0) {
      cerr << "Tpetra RemoveEmptyProcesses test" << endl
           << "Create Map and matrix" << endl;
    }

    // generate problem
    LO nEle;
    if (myRank<commSize-1)
      nEle = 4;
    else
      nEle = 0;
    GO ngEle = 4*(commSize-1);
    RCP<const map_type> map = rcp (new map_type (ngEle, nEle, 0, comm));

      TEUCHOS_TEST_FOR_EXCEPTION(
        ! Kokkos::is_initialized (), std::logic_error,
        "Kokkos is not initialized!" );


    RCP<crs_matrix_type> matrix = rcp (new crs_matrix_type (map, 1));
    const LO NumMyElements = map->getNodeNumElements ();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList ();

    // Make the matrix the identity matrix.
    if (myRank == 0) {
      cerr << "Fill matrix by calling insertGlobalValues" << endl;
    }
    for (LO i = 0; i < NumMyElements; ++i) {
      matrix->insertGlobalValues (MyGlobalElements[i],
                                  Teuchos::tuple<GO> (MyGlobalElements[i]),
                                  Teuchos::tuple<Scalar> (ONE));
    }

    if (myRank == 0) {
      cerr << "Call fillComplete on the matrix" << endl;
    }

    matrix->fillComplete ();
    TEST_ASSERT( matrix->isFillComplete () );
    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );


    if (myRank == 0) {
      cerr << "Call removeEmptyProcessesInPlace" << endl;
    }

    RCP<const map_type> reducedMap = map->removeEmptyProcesses();
    matrix->removeEmptyProcessesInPlace(reducedMap);


    if (myRank == commSize-1) {
      TEST_ASSERT( reducedMap.is_null() );
      TEST_ASSERT( matrix->getMap().is_null() );
      TEST_ASSERT( matrix->getRowMap().is_null() );
      TEST_ASSERT( matrix->getColMap().is_null() );
      TEST_ASSERT( matrix->getRangeMap().is_null() );
      TEST_ASSERT( matrix->getDomainMap().is_null() );
    } else {
      TEST_ASSERT( matrix->getMap()->getComm()->getSize() == commSize-1 );
      TEST_ASSERT( matrix->getRowMap()->getComm()->getSize() == commSize-1 );
      TEST_ASSERT( matrix->getColMap()->getComm()->getSize() == commSize-1 );
      TEST_ASSERT( matrix->getRangeMap()->getComm()->getSize() == commSize-1 );
      TEST_ASSERT( matrix->getDomainMap()->getComm()->getSize() == commSize-1 );
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );


    if (myRank == 0) {
      cerr << "All done!  Test " << (gblSuccess ? "succeeded" : "FAILED") << "."
           << endl;
    }
  }

//
// Instantiations of tests
//
#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, RemoveEmptyProcesses, SCALAR, LO, GO, NODE )

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

} // namespace (anonymous)
