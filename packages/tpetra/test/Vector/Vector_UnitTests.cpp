#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_Vector.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;
  using Tpetra::Vector;
  using Teuchos::OrdinalTraits;
  using std::endl;

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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, basic, Ordinal, Scalar )
  {
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal NEGONE = ZERO - OrdinalTraits<Ordinal>::one();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    Map<Ordinal> map(NEGONE,indexBase,numLocal,platform);
    Vector<Ordinal,Scalar> vec(map);
    out << vec << endl;
  }



  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  #define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Vector, basic, ORDINAL, double )

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

      // FINISH: add complex tests

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Vector, basic, ORDINAL, char ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Vector, basic, ORDINAL, int ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Vector, basic, ORDINAL, float ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Vector, basic, ORDINAL, double )

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
