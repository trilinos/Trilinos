#include <string>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include "Amesos2_Factory.hpp"
#include "Amesos2_Util_is_same.hpp"

namespace {

using std::cout;
using std::endl;
using std::string;

using Teuchos::as;
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
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;

using Amesos::MatrixAdapter;
using Amesos::MultiVecAdapter;
using Amesos::Superlu;

using Amesos::Util::is_same;


typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

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
    ret = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Superlu, Initialization, SCALAR, LO, GO )
{
  /* Test correct initialization of the MatrixAdapter
   *
   * - All Constructors
   * - Correct initialization of class members
   * - Correct typedefs ( using Amesos::is_same<> )
   */
  typedef ScalarTraits<SCALAR> ST;
  typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
  typedef MultiVector<SCALAR,LO,GO,Node> MV;
  typedef Superlu<MAT,MV> SOLVER;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm = getDefaultComm();
  //const size_t numprocs = comm->getSize();
  const size_t rank     = comm->getRank();
  // create a Map
  const size_t numLocal = 10;
  RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );
  RCP<MAT> eye = rcp( new MAT(map,1) );
  GO base = numLocal*rank;
  for( size_t i = 0; i < numLocal; ++i ){
    eye->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<SCALAR>(ST::one()));
  }
  eye->fillComplete();

  // Create X
  RCP<MV> X = rcp(new MV(map,11));
  X->randomize();

  // Create B
  RCP<MV> B = rcp(new MV(map,11));
  B->randomize();

  // Constructor from Factory
  RCP<Amesos::SolverBase> solver = Amesos::Factory<MAT,MV>::create("Superlu",eye,X,B);

  TEST_ASSERT( solver->getNumSymbolicFact() == 0 );
  TEST_ASSERT( solver->getNumNumericFact() == 0 );
  TEST_ASSERT( solver->getNumSolve() == 0 );

  // The following should all pass at compile time
  TEST_ASSERT( (is_same<MAT,typename SOLVER::matrix_type>::value) );
  TEST_ASSERT( (is_same<MV,typename SOLVER::vector_type>::value) );
  TEST_ASSERT( (is_same<SCALAR,typename SOLVER::scalar_type>::value) );
  TEST_ASSERT( (is_same<LO,typename SOLVER::local_ordinal_type>::value) );
  TEST_ASSERT( (is_same<GO,typename SOLVER::global_ordinal_type>::value) );
  TEST_ASSERT( (is_same<global_size_t,typename SOLVER::global_size_type>::value) );
  TEST_ASSERT( (is_same<Node,typename SOLVER::node_type>::value) );
}


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Superlu, SymbolicFactorization, SCALAR, LO, GO )
{
  typedef ScalarTraits<SCALAR> ST;
  typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
  typedef MultiVector<SCALAR,LO,GO,Node> MV;
  typedef Superlu<MAT,MV> SOLVER;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm = getDefaultComm();
  //const size_t numprocs = comm->getSize();
  const size_t rank     = comm->getRank();
  // create a Map
  const size_t numLocal = 10;
  RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );
  RCP<MAT> eye = rcp( new MAT(map,1) );
  GO base = numLocal*rank;
  for( size_t i = 0; i < numLocal; ++i ){
    eye->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<SCALAR>(ST::one()));
  }
  eye->fillComplete();

  // Create X
  RCP<MV> X = rcp(new MV(map,11));
  X->randomize();

  // Create B
  RCP<MV> B = rcp(new MV(map,11));
  B->randomize();

  // Constructor from Factory
  RCP<Amesos::SolverBase> solver = Amesos::Factory<MAT,MV>::create("Superlu",eye,X,B);

  solver->symbolicFactorization();
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Superlu, NumericFactorization, SCALAR, LO, GO )
{
  typedef ScalarTraits<SCALAR> ST;
  typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
  typedef MultiVector<SCALAR,LO,GO,Node> MV;
  typedef Superlu<MAT,MV> SOLVER;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm = getDefaultComm();
  //const size_t numprocs = comm->getSize();
  const size_t rank     = comm->getRank();
  // create a Map
  const size_t numLocal = 10;
  RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );
  RCP<MAT> eye = rcp( new MAT(map,1) );
  GO base = numLocal*rank;
  for( size_t i = 0; i < numLocal; ++i ){
    eye->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<SCALAR>(ST::one()));
  }
  eye->fillComplete();

  // Create X
  RCP<MV> X = rcp(new MV(map,11));
  X->randomize();

  // Create B
  RCP<MV> B = rcp(new MV(map,11));
  B->randomize();

  // Constructor from Factory
  RCP<Amesos::SolverBase> solver = Amesos::Factory<MAT,MV>::create("Superlu",eye,X,B);

  solver->symbolicFactorization().numericFactorization();

  // Good way to check the factors L and U?  Probs not, since they are private members
}


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Superlu, Solve, SCALAR, LO, GO )
{
  typedef ScalarTraits<SCALAR> ST;
  typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
  typedef MultiVector<SCALAR,LO,GO,Node> MV;
  typedef Superlu<MAT,MV> SOLVER;
  typedef typename ScalarTraits<SCALAR>::magnitudeType MAG;
  const MAG M0 = ScalarTraits<MAG>::zero();
    
  RCP<const Comm<int> > comm = getDefaultComm();
  //const size_t numprocs = comm->getSize();

  const size_t numVectors = 1;

  // create a Map
  global_size_t nrows = 6;
  RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(nrows,0,comm) );
  RCP<MAT> A = rcp( new MAT(map,3) ); // max of three entries in a row

  /* We will do two tests.
   *
   * 1) solving a system with known solution, for which we will be using the
   * following matrix:
   * 
   * [ [ 7,  0,  -3, 0,  -1, 0 ]
   *   [ 2,  8,  0,  0,  0,  0 ]
   *   [ 0,  0,  1,  0,  0,  0 ]
   *   [ -3, 0,  0,  5,  0,  0 ]
   *   [ 0,  -1, 0,  0,  4,  0 ]
   *   [ 0,  0,  0,  -2, 0,  6 ] ]
   *
   * and 2) solving a randomized system
   */
  // Construct matrix
  A->insertGlobalValues(0,tuple<GO>(0,2,4),tuple<SCALAR>(7,-3,-1));
  A->insertGlobalValues(1,tuple<GO>(0,1),tuple<SCALAR>(2,8));
  A->insertGlobalValues(2,tuple<GO>(2),tuple<SCALAR>(1));
  A->insertGlobalValues(3,tuple<GO>(0,3),tuple<SCALAR>(-3,5));
  A->insertGlobalValues(4,tuple<GO>(1,4),tuple<SCALAR>(-1,4));
  A->insertGlobalValues(5,tuple<GO>(3,5),tuple<SCALAR>(-2,6));
  A->fillComplete();

  // Create random X
  RCP<MV> X = rcp(new MV(map,numVectors));
  X->randomize();

  // RCP<FancyOStream> os = getDefaultOStream();
  // X->describe(*os,Teuchos::VERB_EXTREME);

  // Create B
  RCP<MV> B = rcp(new MV(map,tuple<SCALAR>(-7,18,3,17,18,28),nrows,numVectors));

  // Constructor from Factory
  RCP<Amesos::SolverBase> solver = Amesos::Factory<MAT,MV>::create("Superlu",A,X,B);

  solver->symbolicFactorization().numericFactorization().solve();

  /* Check X for correct solution
   *
   * Should be X = [1; 2; 3; 4; 5; 6]
   */
  Array<MAG> norms(numVectors), zeros(numVectors);
  std::fill(zeros.begin(),zeros.end(),M0);

  RCP<MV> xCheck = rcp(new MV(map,tuple<SCALAR>(1,2,3,4,5,6),nrows,numVectors));

  X->update(as<SCALAR>(-1),*xCheck,as<SCALAR>(1));
  X->norm2(norms);
  TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);

  // ETB: Should not use randomness in tests.  It would be better to find a
  // few good sample matrices to test with
  /*
   *  const size_t numVectors_r = 5;
   *
   *  // Create random X
   *  RCP<MV> X_r = rcp(new MV(A->getColMap(), numVectors));
   *  X_r->randomize();
   *
   *  // Create B from A and X
   *  RCP<MV> B_r = rcp(new MV(A->getRowMap(),numVectors));
   *  A->multiply(*X_r, *B_r, Teuchos::NO_TRANS, as<SCALAR>(1), as<SCALAR>(0)); // B_r = A * X_r
   *
   *  RCP<MV> xCheck_r = rcp(new MV(A->getColMap(), numVectors));
   *
   *  // Construct Solver from Factory
   *  RCP<Amesos::SolverBase> solver_r
   *    = Amesos::Factory<MAT,MV>::create("Superlu", A, xCheck_r, B_r);
   *
   *  solver_r->symbolicFactorization().numericFactorization().solve();
   *
   *  // Test equaility of original X and solution
   *  X_r->update(as<SCALAR>(-1),*xCheck_r,as<SCALAR>(1));
   *  X_r->norm2(norms);
   *  TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,0.5);
   */
}


/*
 * Instantiations
 */
#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)\
  typedef std::complex<float>  ComplexFloat;\
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexFloat)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)\
  typedef std::complex<double> ComplexDouble;\
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexDouble) 
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#endif

// Uncomment this for really fast development cycles but make sure to comment
// it back again before checking in so that we can test all the types.
// #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superlu, Initialization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superlu, SymbolicFactorization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superlu, NumericFactorization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superlu, Solve, SCALAR, LO, GO )


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
