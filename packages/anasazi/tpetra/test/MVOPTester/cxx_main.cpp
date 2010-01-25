#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicOutputManager.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Tpetra::MultiVector;
  using Tpetra::CrsMatrix;
  using std::endl;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Anasazi::OutputManager;
  using Anasazi::BasicOutputManager;
  using Anasazi::Warnings;
  using Teuchos::tuple;

  typedef DefaultPlatform::DefaultPlatformType::NodeType Node;

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
      DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  template<class Scalar, class O1, class O2>
  RCP<CrsMatrix<Scalar,O1,O2,Node> > constructDiagMatrix(const RCP<const Map<O1,O2,Node> > &map) 
  {
    RCP<CrsMatrix<Scalar,O1,O2,Node> > op = rcp( new CrsMatrix<Scalar,O1,O2,Node>(map,1) );
    for (size_t i=0; i<map->getNodeNumElements(); ++i) {
      op->insertGlobalValues(map->getGlobalElement(i),tuple(map->getGlobalElement(i)), tuple(ScalarTraits<Scalar>::one()));
    }
    op->fillComplete();
    return op;
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, MVTestDist, O1, O2, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,O1,O2> MV;
    const O2 dim = 500;
    const Teuchos_Ordinal numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new BasicOutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<Map<O1,O2,Node> > map = rcp( new Map<O1,O2,Node>(dim,0,comm) );
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Anasazi::TestMultiVecTraits<Scalar,MV>(MyOM,mvec);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, MVTestLocal, O1, O2, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,O1,O2> MV;
    const O2 dim = 500;
    const Teuchos_Ordinal numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new BasicOutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<Map<O1,O2,Node> > map = rcp(new Map<O1,O2,Node>(dim,0,comm,Tpetra::LocallyReplicated) );
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Anasazi::TestMultiVecTraits<Scalar,MV>(MyOM,mvec);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, OPTestLocal, O1, O2, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,O1,O2> MV;
    typedef Tpetra::Operator<Scalar,O1,O2>    OP;
    const O2 dim = 500;
    const Teuchos_Ordinal numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new BasicOutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<Map<O1,O2,Node> > map = rcp(new Map<O1,O2,Node>(dim,0,comm,Tpetra::LocallyReplicated) );
    // create a CrsMatrix
    RCP<OP> op = constructDiagMatrix<Scalar,O1,O2>(map);
    // create a multivector
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Anasazi::TestOperatorTraits<Scalar,MV,OP>(MyOM,mvec,op);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, OPTestDist, O1, O2, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,O1,O2> MV;
    typedef Tpetra::Operator<Scalar,O1,O2>    OP;
    const O2 dim = 500;
    const Teuchos_Ordinal numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new BasicOutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<Map<O1,O2,Node> > map = rcp( new Map<O1,O2,Node>(dim,0,comm) );
    // create a CrsMatrix
    RCP<OP> op = constructDiagMatrix<Scalar,O1,O2>(map);
    // create a multivector
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Anasazi::TestOperatorTraits<Scalar,MV,OP>(MyOM,mvec,op);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  //
  // INSTANTIATIONS
  //

#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(O1, O2)\
     typedef std::complex<float> ComplexFloat; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(O1, O2, ComplexFloat)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(O1, O2)\
     typedef std::complex<double> ComplexDouble; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(O1, O2, ComplexDouble)
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(O1, O2)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(O1, O2)
#endif

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( O1, O2, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, MVTestDist, O1, O2, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, MVTestLocal, O1, O2, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, OPTestDist, O1, O2, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, OPTestLocal, O1, O2, SCALAR )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
    UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( O1, O2 ) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(O1, O2) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(O1, O2, double)

     UNIT_TEST_GROUP_ORDINAL(int)
# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( O1, O2 ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(O1, O2, float) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(O1, O2, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(O1, O2) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(O1, O2)

     typedef short int ShortInt;
     UNIT_TEST_GROUP_ORDINAL(ShortInt)

     UNIT_TEST_GROUP_ORDINAL(int)

     typedef long int LongInt;
     UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongInt )
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt;
        UNIT_TEST_GROUP_ORDINAL_ORDINAL( int,LongLongInt)
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
