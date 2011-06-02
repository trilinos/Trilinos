#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include "Tpetra_Map_def.hpp"
#include "Tpetra_Directory_def.hpp"
#include "Tpetra_CrsGraph_def.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_MatrixIO_def.hpp"
#include "Tpetra_MultiVector_def.hpp"
#include "Tpetra_Vector_def.hpp"
#endif

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;
  using Teuchos::null;
  using Tpetra::Import;
  using Tpetra::Map;
  using Tpetra::CrsGraph;
  using Tpetra::CrsMatrix;
  using Tpetra::Vector;
  using Tpetra::global_size_t;
  using Tpetra::DefaultPlatform;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignord and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  RCP<DefaultPlatform::DefaultPlatformType::NodeType> getDefaultNode()
  {
    return DefaultPlatform::getDefaultPlatform().getNode();
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST( readHBMatrix, Bug5072_ReadOneRowMPI )
  {
    // failure reading 1x4 matrix under MPI
    typedef int                       LO;
    typedef int                       GO;
    typedef DefaultPlatform::DefaultPlatformType::NodeType Node;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    RCP<Node>             node = getDefaultNode();
    RCP<const CrsMatrix<double,int> > readMatrix, testMatrix;
    {
      // this is what the file looks like: 1 3 4 9
      RCP<const Map<int> > rng = Tpetra::createUniformContigMap<int,int>(1,comm);
      RCP<const Map<int> > dom = Tpetra::createUniformContigMap<int,int>(4,comm);
      RCP<CrsMatrix<double,int> > A = Tpetra::createCrsMatrix<double>(rng);
      if (myImageID == 0) {
        A->insertGlobalValues( 0, Teuchos::tuple<int>(0,1,2,3), Teuchos::tuple<double>(1.0,3.0,4.0,9.0) );
      }
      A->fillComplete(dom,rng,Tpetra::DoOptimizeStorage);
      testMatrix = A;
    }
    {
      RCP<CrsMatrix<double,int> > A;
      Tpetra::Utils::readHBMatrix("addA2.hb", comm, node, A);
      readMatrix = A;
    }
    // test that *readMatrix == *testMatrix 
    TEST_EQUALITY( testMatrix->getNodeNumRows(), readMatrix->getNodeNumRows() );
    TEST_EQUALITY( testMatrix->getNodeNumCols(), readMatrix->getNodeNumCols() );
    TEST_EQUALITY( testMatrix->getNodeNumEntries(), readMatrix->getNodeNumEntries() );
    if (success) {
      Teuchos::ArrayView<const int>    rowinds1, rowinds2;
      Teuchos::ArrayView<const double> rowvals1, rowvals2;
      for (int r=0; r < (int)testMatrix->getNodeNumRows(); ++r ) {
        testMatrix->getLocalRowView(r, rowinds1, rowvals1);  
        readMatrix->getLocalRowView(r, rowinds2, rowvals2);  
        TEST_COMPARE_ARRAYS( rowinds1, rowinds2 );
        TEST_COMPARE_ARRAYS( rowvals1, rowvals2 );
      }
    }
    // All procs fail if any proc fails 
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  // TEUCHOS_UNIT_TEST( Map, Bug4756_UnsignedGlobalOrdinal )
  // {
  //   // test bug where unsigned LO clashes with signed GO
  //   typedef unsigned int LO;
  //   // this still has to be bigger than LO
  //   typedef     long int GO;
  //   // this test assumes that global_size_t (default: size_t) is larger than GO=="long int", which should be true on 64-bit builds.
  //   // create a comm  
  //   RCP<const Comm<int> > comm = getDefaultComm();
  //   const global_size_t numGlobal = comm->getSize();
  //   RCP<const Map<LO,GO> > map = rcp(new Map<LO,GO>(numGlobal,0,comm) );
  //   RCP<const CrsGraph<LO,GO> > graph = rcp(new CrsGraph<LO,GO>(map,0,Tpetra::DynamicProfile) );
  //   TEST_EQUALITY_CONST( graph != null, true );
  //   // All procs fail if any proc fails 
  //   int globalSuccess_int = -1;
  //   reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  //   TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  // }

  ////
  TEUCHOS_UNIT_TEST( DistObject, Bug5129_OverlyStrictMapComparison )
  {
    // test bug where map test checks that maps are the same object, instead of checking that maps are equivalent
    typedef int LO;
    // this still has to be bigger than LO
    typedef int GO;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const global_size_t numGlobal = comm->getSize()*10;
    // create two separate, but identical, maps
    RCP<const Map<LO,GO> > mapImportIn  = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    RCP<const Map<LO,GO> > mapImportOut = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    RCP<const Map<LO,GO> > mapIn        = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    RCP<const Map<LO,GO> > mapOut       = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    TEST_EQUALITY_CONST( *mapImportIn  == *mapIn,  true );
    TEST_EQUALITY_CONST( *mapImportOut == *mapOut, true );
    TEST_EQUALITY_CONST( mapImportIn   == mapIn,  false );
    TEST_EQUALITY_CONST( mapImportOut  == mapOut, false );
    // create import, vectors from these maps
    RCP<const Import<LO,GO> > import = Tpetra::createImport(mapImportIn, mapImportOut);
    RCP<Vector<double,LO,GO> > vecIn = Tpetra::createVector<double>(mapIn);
    RCP<Vector<double,LO,GO> > vecOut = Tpetra::createVector<double>(mapOut);
    // do the import; under the bug, this should throw an exception
    TEST_NOTHROW( vecOut->doImport( *vecIn, *import, Tpetra::REPLACE ) )
    // All procs fail if any proc fails 
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

}
