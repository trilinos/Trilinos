#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "BelosConfigDefs.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosOutputManager.hpp"

namespace {

  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Tpetra::LocallyReplicated;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::MultiVector;
  using Tpetra::CrsMatrix;
  using Tpetra::Operator;
  using Tpetra::global_size_t;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using std::endl;
  using Belos::OutputManager;
  using Belos::Warnings;

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
    RCP<const Comm<int> > ret;
    if (testMpi) {
      ret = DefaultPlatform::getDefaultPlatform().getComm();
    }
    else {
     ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  template<class Scalar, class LO, class GO>
  RCP<CrsMatrix<Scalar,LO,GO> > constructDiagMatrix(const RCP<const Map<LO,GO> > &map) 
  {
    // create identity matrix
    RCP<CrsMatrix<Scalar,LO,GO> > op = rcp( new CrsMatrix<Scalar,LO,GO>(map,1) );
    for (size_t lid=0; lid < map->getNodeNumElements(); ++lid) {
      op->insertGlobalValues(map->getGlobalElement(lid),tuple(map->getGlobalElement(lid)), tuple(ScalarTraits<Scalar>::one()));
    }
    op->fillComplete();
    return op;
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, MVTestDist, LO, GO, Scalar )
  {
    typedef MultiVector<Scalar,LO,GO> MV;
    const global_size_t dim = 500;
    const size_t numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcpFromRef(out)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<const Map<LO,GO> > map = rcp( new Map<LO,GO>(dim,static_cast<GO>(0),comm) );
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestMultiVecTraits<Scalar,MV>(MyOM,mvec);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, MVTestLocal, LO, GO, Scalar )
  {
    typedef MultiVector<Scalar,LO,GO> MV;
    const global_size_t dim = 500;
    const size_t numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcpFromRef(out)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<const Map<LO,GO> > map = rcp(new Map<LO,GO>(dim,static_cast<GO>(0),comm,LocallyReplicated) );
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestMultiVecTraits<Scalar,MV>(MyOM,mvec);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, OPTestLocal, LO, GO, Scalar )
  {
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef Operator<Scalar,LO,GO>    OP;
    const global_size_t dim = 500;
    const size_t numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcpFromRef(out)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map (local)
    RCP<const Map<LO,GO> > map = rcp( new Map<LO,GO>(dim,static_cast<GO>(0),comm,LocallyReplicated) );
    // create a CrsMatrix
    RCP<OP> op = constructDiagMatrix<Scalar,LO,GO>(map);
    // create a multivector
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestOperatorTraits<Scalar,MV,OP>(MyOM,mvec,op);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, OPTestDist, LO, GO, Scalar )
  {
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef Operator<Scalar,LO,GO>    OP;
    const global_size_t dim = 500;
    const size_t numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcpFromRef(out)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<const Map<LO,GO> > map = rcp( new Map<LO,GO>(dim,static_cast<GO>(0),comm) );
    // create a CrsMatrix
    RCP<OP> op = constructDiagMatrix<Scalar,LO,GO>(map);
    // create a multivector
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestOperatorTraits<Scalar,MV,OP>(MyOM,mvec,op);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  //
  // INSTANTIATIONS
  //

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
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, MVTestDist, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, MVTestLocal, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, OPTestDist, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, OPTestLocal, LO, GO, SCALAR )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
    UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO ) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, double)

     UNIT_TEST_GROUP_ORDINAL(int)
# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, float) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)

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
