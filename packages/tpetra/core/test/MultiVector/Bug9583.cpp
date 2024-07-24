// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Core.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_CommHelpers.hpp"

namespace { // (anonymous)

  // Test for Bug #9583
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MultiVector, Bug9583, S, LO, GO, NODE)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Tpetra::Map<LO, GO, NODE> map_type;
    typedef Tpetra::MultiVector<S, LO, GO, NODE> vec_type;

    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();

    // Create a Map which has zero rows on ranks >= gblNumRows
    RCP<const map_type> rowMap;
    const int gblNumRows = 3;
    const GO indexBase = 0;

    size_t numLocalElements = 1;
    if(myRank >= gblNumRows) {
        numLocalElements = 0;
    }
    rowMap = rcp (new map_type (gblNumRows, numLocalElements, indexBase, comm));

    // Create distributed zero vectors
    vec_type x (rowMap, 1, true);
    vec_type y (rowMap, 1, true);

    // Create a local vector z -- one entry that will be same on all ranks
    const size_t numElm = 1;
    RCP<const map_type > rowMapLocal = 
        Tpetra::createLocalMapWithNode<LO,GO,NODE>( numElm, comm);
    vec_type z (rowMapLocal, 1, true);

    // Fill z with a non-zero number
    z.putScalar(static_cast<S>(999));

    // z = x^T * y
    z.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, x, y, 0.0);

    // Since x and y are zero vectors, we expect z to also be zero
    // With bug #9583, z has value 999 on ranks >= gblNumRows
    // With fix, z has value zero on all ranks.
    TEST_EQUALITY( z.getData(0)[0], static_cast<S>(0.0) )
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MultiVector, Bug9583, SCALAR, LO, GO, NODE)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)


