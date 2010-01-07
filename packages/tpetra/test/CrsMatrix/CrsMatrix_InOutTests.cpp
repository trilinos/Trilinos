#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace {

  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;

  bool testMpi = true;
  string filedir;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected matrix files.");
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
    RCP<const Comm<int> > ret;
    if (testMpi) {
      ret = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    }
    else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  //
  // UNIT TESTS
  // 



  // 
  // INSTANTIATIONS
  //


  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_ORDINAL_SCALAR(LO, GO, SCALAR)

#ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#   define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
           UNIT_TEST_ORDINAL_SCALAR(LO, GO, float)
    UNIT_TEST_GROUP(int, int)
#else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD
#   define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
           UNIT_TEST_ORDINAL_SCALAR(LO, GO, float)
    typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt, int)
    UNIT_TEST_GROUP_ORDINAL(int, int)
#endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
