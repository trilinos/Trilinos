#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_DefaultPlatform.hpp>

#include "Amesos2_TpetraMultiVecAdapter.hpp"
#include "Amesos2_Util_is_same.hpp"

namespace {

using std::cout;
using std::endl;
using std::string;

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using Teuchos::ScalarTraits;
using Teuchos::OrdinalTraits;
using Teuchos::FancyOStream;
using Teuchos::VerboseObjectBase;

using Tpetra::global_size_t;
using Tpetra::MultiVector;
using Tpetra::DefaultPlatform;
using Tpetra::Map;

using Amesos::MultiVecAdapter;

using Amesos::Util::is_same;


typedef DefaultPlatform::DefaultPlatformType::NodeType Node;

bool testMpi = false;

// Where to look for input files
string filedir;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("filedir",&filedir,"Directory of matrix files.");
  clp.addOutputSetupOptions(true);
  clp.setOption("test-mpi", "test-serial", &testMpi,
    "Test Serial by default (for now) or force MPI test.  In a serial build,"
    " this option is ignored and a serial comm is always used." );
}

RCP<const Comm<int> > getDefaultComm()
{
  RCP<const Comm<int> > ret;
  if( testMpi ){
    ret = DefaultPlatform::getDefaultPlatform().getComm();
  } else {
    ret = rcp(new Teuchos::SerialComm<int>());
  }
  return ret;
}

RCP<FancyOStream> getDefaultOStream()
{
  return( VerboseObjectBase::getDefaultOStream() );
}

/*
 * UNIT TESTS
 */

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Initialization, SCALAR, LO, GO )
{
  /* Test correct initialization of the MultiVecAdapter
   *
   * - All Constructors
   * - Correct initialization of class members
   * - Correct typedefs ( using Amesos::is_same<> )
   */
  typedef ScalarTraits<SCALAR> ST;
  typedef MultiVector<SCALAR,LO,GO,Node> MV;
  typedef MultiVecAdapter<MV> ADAPT;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm = getDefaultComm();
  const size_t numprocs = comm->getSize();
  const size_t rank     = comm->getRank();
  // create a Map
  const size_t numLocal = 10;
  RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

  RCP<MV> mv = rcp(new MV(map,11));
  mv->randomize();
  // RCP<FancyOStream> os = getDefaultOStream();
  // mv->describe(*os,Teuchos::VERB_EXTREME);

  RCP<ADAPT> adapter = rcp(new ADAPT(mv));
  // Test copy constructor
  RCP<ADAPT> adapter2 = rcp(new ADAPT(*adapter));

  // Check that the values remain the same (more comprehensive test of get1dView elsewhere...
  TEST_EQUALITY( mv->get1dViewNonConst(),      adapter->get1dViewNonConst() );
  TEST_EQUALITY( adapter->get1dViewNonConst(), adapter2->get1dViewNonConst() );

  // The following should all pass at compile time
  TEST_ASSERT( (is_same<SCALAR,        typename ADAPT::scalar_type>::value) );
  TEST_ASSERT( (is_same<LO,            typename ADAPT::local_ordinal_type>::value) );
  TEST_ASSERT( (is_same<GO,            typename ADAPT::global_ordinal_type>::value) );
  TEST_ASSERT( (is_same<Node,          typename ADAPT::node_type>::value) );
  TEST_ASSERT( (is_same<global_size_t, typename ADAPT::global_size_type>::value) );
  TEST_ASSERT( (is_same<MV,            typename ADAPT::multivec_type>::value) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Dimensions, SCALAR, LO, GO )
{
  // Test that the dimensions reported by the adapter match those as reported
  // by the Tpetra::MultiVector
  typedef ScalarTraits<SCALAR> ST;
  typedef MultiVector<SCALAR,LO,GO,Node> MV;
  typedef MultiVecAdapter<MV> ADAPT;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm = getDefaultComm();
  const size_t numprocs = comm->getSize();
  const size_t rank     = comm->getRank();
  // create a Map
  const size_t numLocal = 10;
  RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

  RCP<MV> mv = rcp(new MV(map,11));
  mv->randomize();
  // RCP<FancyOStream> os = getDefaultOStream();
  // mv->describe(*os,Teuchos::VERB_EXTREME);

  RCP<ADAPT> adapter = rcp(new ADAPT(mv));

  TEST_EQUALITY( mv->getLocalLength(),  adapter->getLocalLength()      );
  TEST_EQUALITY( mv->getNumVectors(),   adapter->getLocalNumVectors()  );
  TEST_EQUALITY( mv->getNumVectors(),   adapter->getGlobalNumVectors() );
  TEST_EQUALITY( mv->getGlobalLength(), adapter->getGlobalLength()     );
  TEST_EQUALITY( mv->getStride(),       adapter->getStride()           );
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Copy, SCALAR, LO, GO )
{
  /* Test the get1dCopy() method of MultiVecAdapter.  We can check against a
   * known multivector and also check against what is returned by the
   * Tpetra::MultiVector.
   */
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, View, SCALAR, LO, GO )
{
  /* Test the get1dViewNonConst() method of MultiVecAdapter.  We can check against a
   * known multivector and also check against what is returned by the
   * Tpetra::MultiVector.
   */
  typedef ScalarTraits<SCALAR> ST;
  typedef typename ScalarTraits<SCALAR>::magnitudeType MAG;
  typedef MultiVector<SCALAR,LO,GO,Node> MV;
  typedef MultiVecAdapter<MV> ADAPT;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const MAG M0 = ScalarTraits<MAG>::zero();
  const SCALAR S0 = ST::zero();
  RCP<const Comm<int> > comm = getDefaultComm();
  const size_t numprocs = comm->getSize();
  const size_t rank     = comm->getRank();
  // create a Map
  const size_t numVectors = 7;
  const size_t numLocal = 13;
  RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

  RCP<MV> mv = rcp(new MV(map,numVectors));
  mv->randomize();
  // RCP<FancyOStream> os = getDefaultOStream();
  // mv->describe(*os,Teuchos::VERB_EXTREME);

  RCP<ADAPT> adapter = rcp(new ADAPT(mv));

  // Check that the values remain the same
  TEST_EQUALITY( mv->get1dViewNonConst(), adapter->get1dViewNonConst() );

  // Check that direct changes to the wrapped object are reflected in the adapter
  mv->randomize();
  TEST_EQUALITY( mv->get1dViewNonConst(), adapter->get1dViewNonConst() );

  // check that 1dView and 1dCopy have the same values
  ArrayRCP<SCALAR> view;
  view = adapter->get1dViewNonConst();

  // clear view, ensure that mv is zero
  std::fill(view.begin(), view.end(), S0);
  view = Teuchos::null;
  Array<MAG> norms(numVectors), zeros(numVectors,M0);
  mv->norm2(norms());
  TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
}

/* Also test the updateValues() method, once it is implemented.  It should
 * take a representation from either the get1dCopy or get2dCopy methods and
 * place the values back into the matrix (either serial or distributed)
 */


/*
 * Instantiations
 */

#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)\
  typedef std::complex<float> ComplexFloat; \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexFloat)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)\
  typedef std::complex<double> ComplexDouble; \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexDouble) 
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#endif

// Uncomment this for really fast development cycles but make sure to comment
// it back again before checking in so that we can test all the types.
// #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Initialization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Dimensions, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Copy, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, View, SCALAR, LO, GO ) 
  
#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO ) \
  UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, double) \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT( LO, GO )
UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO ) \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, float)  \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, double) \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO) \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO,GO)

UNIT_TEST_GROUP_ORDINAL(int)

typedef long int LongInt;
UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongInt )
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long int LongLongInt;
UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongLongInt )
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

} // end anonymous namespace
