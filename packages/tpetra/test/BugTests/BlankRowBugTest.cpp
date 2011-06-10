#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp>
namespace {

  using Teuchos::RCP;
  using Teuchos::Comm;
  using Tpetra::tuple;

  TEUCHOS_UNIT_TEST(CrsMatrix, BlankRowImport)
  {
    typedef Tpetra::Map<int,int>                    Map;
    typedef Tpetra::CrsMatrix<double,int,int> CrsMatrix;
    typedef Tpetra::Import<int,int>              Import;
    RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    // We run this test explicitly in MPI mode with 2 processors as described in
    // the CMakeLists.txt file. This is just asserting that fact.
    TEST_EQUALITY_CONST(comm->getSize(), 2)

    const int rank = comm->getRank();
    RCP<const Map> destRowMap, sourceRowMap;
    if (rank ==0) {
      sourceRowMap = Tpetra::createNonContigMap<int,int>( tuple<int>(0), comm );
    } else {
      sourceRowMap = Tpetra::createNonContigMap<int,int>( tuple<int>(1), comm );
    }
    destRowMap   = Tpetra::createNonContigMap<int,int>( tuple<int>(0,1), comm );

    RCP<CrsMatrix> srcMat = Tpetra::createCrsMatrix<double>(sourceRowMap);
    if (rank == 0) {
      srcMat->insertGlobalValues(0, tuple<int>(0), tuple<double>(1.0) );
    }
    srcMat->fillComplete();
    /* 
       srcMat = [1 ] // proc 0
                [  ] // proc 1
     */
    if (rank == 0) {
      TEST_EQUALITY_CONST( srcMat->getNumEntriesInGlobalRow(0), 1 );
    } else {
      TEST_EQUALITY_CONST( srcMat->getNumEntriesInGlobalRow(1), 0 );
    }

    RCP<CrsMatrix> dstMat = Tpetra::createCrsMatrix<double>(destRowMap);
    RCP<const Import> importer = Tpetra::createImport(sourceRowMap, destRowMap);
    // global row 1 in srcMat is empty: this is a null communication to dstMat
    dstMat->doImport(*srcMat, *importer, Tpetra::INSERT);
    /* 
       dstMat_p0 = [1 ]
                   [  ]
       dstMat_p1 = [1 ]
                   [  ]
    */
    TEST_EQUALITY_CONST( dstMat->getNumEntriesInGlobalRow(0), 1 );
    TEST_EQUALITY_CONST( dstMat->getNumEntriesInGlobalRow(1), 0 );
  }

}
