// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Core.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_UnitTestHarness.hpp>

namespace {

  using Teuchos::RCP;
  using Teuchos::Comm;
  using Teuchos::tuple;

  TEUCHOS_UNIT_TEST(CrsMatrix, BlankRowImport)
  {
    typedef Tpetra::Map<>                     map_type;
    typedef map_type::local_ordinal_type      LO;
    typedef map_type::global_ordinal_type     GO;
    typedef Tpetra::CrsMatrix<>::scalar_type  SC;
    typedef Tpetra::CrsMatrix<>               crs_matrix_type;
    typedef Tpetra::Import<LO, GO>            import_type;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    // We run this test explicitly in MPI mode with 2 processors as described in
    // the CMakeLists.txt file. This is just asserting that fact.
    TEST_EQUALITY_CONST(comm->getSize(), 2)

    const int rank = comm->getRank();
    RCP<const map_type> destRowMap, sourceRowMap;
    if (rank == 0) {
      sourceRowMap = Tpetra::createNonContigMap<LO, GO> (tuple<GO> (0), comm);
    }
    else {
      sourceRowMap = Tpetra::createNonContigMap<LO, GO> (tuple<GO> (1), comm);
    }
    destRowMap = Tpetra::createNonContigMap<LO, GO> (tuple<GO> (0, 1), comm);

    RCP<crs_matrix_type> srcMat = Tpetra::createCrsMatrix<SC>(sourceRowMap, 1);
    if (rank == 0) {
      srcMat->insertGlobalValues (static_cast<GO> (0), tuple<GO> (0), tuple<SC> (1.0));
    }
    srcMat->fillComplete();
    /*
       srcMat = [1 ] // proc 0
                [  ] // proc 1
     */
    if (rank == 0) {
      TEST_EQUALITY_CONST( srcMat->getNumEntriesInGlobalRow(0), static_cast<size_t> (1) );
    } else {
      TEST_EQUALITY_CONST( srcMat->getNumEntriesInGlobalRow(1), static_cast<size_t> (0) );
    }

    RCP<crs_matrix_type> dstMat = Tpetra::createCrsMatrix<SC> (destRowMap, 1);
    RCP<const import_type> importer = Tpetra::createImport (sourceRowMap, destRowMap);
    // global row 1 in srcMat is empty: this is a null communication to dstMat
    dstMat->doImport (*srcMat, *importer, Tpetra::INSERT);
    /*
       dstMat_p0 = [1 ]
                   [  ]
       dstMat_p1 = [1 ]
                   [  ]
    */
    TEST_EQUALITY_CONST( dstMat->getNumEntriesInGlobalRow(0), static_cast<size_t> (1) );
    TEST_EQUALITY_CONST( dstMat->getNumEntriesInGlobalRow(1), static_cast<size_t> (0) );
  }

}


