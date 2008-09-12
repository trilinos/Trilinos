#include "Teuchos_UnitTestHarness.hpp"

#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Tpetra::Map;
  using Tpetra::Import;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

#define PRINT_VECTOR(v) \
   { \
     out << #v << ": "; \
     copy(v.begin(), v.end(), ostream_iterator<Ordinal>(out," ")); \
     out << endl; \
   }

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

  template<class Ordinal>
  RCP<const Platform<Ordinal> > getDefaultPlatform()
  {
    if (testMpi) {
      return DefaultPlatform<Ordinal>::getPlatform();
    }
    return rcp(new Tpetra::SerialPlatform<Ordinal>());
  }

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Import, basic, Ordinal )
  {
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // create a platform  
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    // create a comm  
    //RCP<Comm<Ordinal> > comm = platform.createComm();
    //const int numImages = comm->getSize();
    //const int myImageID = comm->getRank();
    // create Maps
    Map<Ordinal> source(as<Ordinal>(-1),as<Ordinal>(10),ZERO,*platform),
                 target(as<Ordinal>(-1),as<Ordinal>(5) ,ZERO,*platform);
    // create Import object
    Import<Ordinal> importer(source, target);
    
    Ordinal same = importer.getNumSameIDs();
    Ordinal permute = importer.getNumPermuteIDs();
    Ordinal remote = importer.getNumRemoteIDs();
    Ordinal sum = same + permute + remote;
    Ordinal expectedSum = target.getNumMyEntries();
    TEST_EQUALITY( sum, expectedSum );
  }



  //
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
// # define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Import, basic, ORDINAL )

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Import, basic, ORDINAL )

    typedef short int ShortInt;
    UNIT_TEST_GROUP_ORDINAL(ShortInt)
    typedef long int LongInt;
    UNIT_TEST_GROUP_ORDINAL(LongInt)
#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}

