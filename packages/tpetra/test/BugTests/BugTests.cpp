#include <Teuchos_UnitTestHarness.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include "Tpetra_Map_def.hpp"
#include "Tpetra_Directory_def.hpp"
#include "Tpetra_CrsGraph_def.hpp"
#endif

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;
  using Teuchos::null;
  using Tpetra::Map;
  using Tpetra::CrsGraph;
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

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST( Map, Bug4756_UnsignedGlobalOrdinal )
  {
    // test bug where unsigned LO clashes with signed GO
    typedef unsigned int LO;
    // this still has to be bigger than LO
    typedef     long int GO;
    // this test assumes that global_size_t (default: size_t) is larger than GO=="long int", which should be true on 64-bit builds.
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const global_size_t numGlobal = comm->getSize();
    RCP<const Map<LO,GO> > map = rcp(new Map<LO,GO>(numGlobal,0,comm) );
    RCP<const CrsGraph<LO,GO> > graph = rcp(new CrsGraph<LO,GO>(map,0,Tpetra::DynamicProfile) );
    TEST_EQUALITY_CONST( graph != null, true );
    // All procs fail if any proc fails 
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

}
