#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"

namespace Teuchos {
  template <>
    ScalarTraits<int>::magnitudeType
    relErr( const int &s1, const int &s2 )
    {
      typedef ScalarTraits<int> ST;
      return ST::magnitude(s1-s2);
    }

  template <>
    ScalarTraits<char>::magnitudeType
    relErr( const char &s1, const char &s2 )
    {
      typedef ScalarTraits<char> ST;
      return ST::magnitude(s1-s2);
    }
}

namespace {

  using Teuchos::OrdinalTraits;
  using Teuchos::RCP;
  using Teuchos::Comm;
  using Tpetra::Platform;
  using Tpetra::DefaultPlatform;
  using Teuchos::rcp;

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

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Platform, basic, Ordinal )
  {
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    TEST_EQUALITY( myImageID < numImages, true );
  }


  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Platform, basic, ORDINAL )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
     UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD
     UNIT_TEST_GROUP_ORDINAL(int)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)
     typedef long int LongInt;   UNIT_TEST_GROUP_ORDINAL(LongInt)
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt; UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
