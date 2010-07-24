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
#include <Tpetra_Map.hpp>

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
using Teuchos::rcpFromRef;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using Teuchos::ScalarTraits;
using Teuchos::OrdinalTraits;
using Teuchos::FancyOStream;
using Teuchos::VerboseObjectBase;
using Teuchos::ETransp;
using Teuchos::EUplo;
using Teuchos::LOWER_TRI;
using Teuchos::UPPER_TRI;
using Teuchos::CONJ_TRANS;
using Teuchos::TRANS;
using Teuchos::NO_TRANS;


using Tpetra::global_size_t;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;
using Tpetra::createContigMap;

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
  typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
  typedef ScalarTraits<SCALAR> ST;
  typedef MultiVector<SCALAR,LO,GO,Node> MV;
  typedef typename ST::magnitudeType Mag;
  typedef ScalarTraits<Mag> MT;
  const size_t numLocal = 13, numVecs = 7;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();
  // create a Map
  RCP<const Map<LO,GO,Node> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
  SCALAR SONE = ST::one();

  /* Create one of the following locally triangular matries:

     0  [1 2       ]
     1  [  1 3     ]
     .  [    .  .  ] = U
     n-2 [       1 n]
     n-1 [         1]

     0  [1           ]
     1  [2 1         ]
     .  [   .  .     ] = L
     n-2 [     n-1 1  ]
     n-1 [         n 1]

     Global matrices are diag(U,U,...,U) and diag(L,L,...,L)

     For each of these, we test with transpose and non-transpose solves.
  */

  MV X(map,numVecs), B(map,numVecs), Xhat(map,numVecs);
  X.setObjectLabel("X");
  B.setObjectLabel("B");
  Xhat.setObjectLabel("Xhat");
  X.randomize();
  for (size_t tnum=0; tnum < 4; ++tnum) {
    EUplo   uplo  = ((tnum & 1) == 1 ? UPPER_TRI  : LOWER_TRI);
    ETransp trans = ((tnum & 2) == 2 ? CONJ_TRANS : NO_TRANS);
    RCP<MAT> AMat;
    {
      // can let the matrix compute a column map
      AMat = rcp(new MAT(map,2));
      // fill the matrix
      if (uplo == UPPER_TRI) {
        for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
          if (gid == map->getMaxGlobalIndex()) {
            AMat->insertGlobalValues( gid, tuple<GO>(gid), tuple<SCALAR>(SONE) );
          }
          else {
            AMat->insertGlobalValues( gid, tuple<GO>(gid,gid+1), tuple<SCALAR>(SONE,as<GO>(gid+2)) );
          }
        }
      }
      else { // uplo == LOWER_TRI
        for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
          if (gid == map->getMinGlobalIndex()) {
            AMat->insertGlobalValues( gid, tuple<GO>(gid), tuple<SCALAR>(SONE) );
          }
          else {
            AMat->insertGlobalValues( gid, tuple<GO>(gid-1,gid), tuple<SCALAR>(as<GO>(gid+1),SONE) );
          }
        }
      }
      AMat->fillComplete();	// optimize storage
    }
    B.randomize();
    AMat->apply(X,B,trans,1.0,0.0);

    Xhat.randomize();

    // Solve A*Xhat = B for Xhat using the Superlu solver
    RCP<Amesos::SolverBase> solver
      = Amesos::Factory<MAT,MV>::create(
        "Superlu",
        AMat,
        rcpFromRef(Xhat),
        rcpFromRef(B) );

    Teuchos::ParameterList params;
    if( trans == CONJ_TRANS ){
      params.set("Trans","CONJ","Whether to solve with transpose");
    } else {			// trans == NO_TRANS
      params.set("Trans","NOTRANS","Whether to solve with transpose");
    }

    solver->setParameters( rcpFromRef(params) );
    solver->symbolicFactorization().numericFactorization().solve();

    // Check result of solve
    Xhat.update(-SONE,X,SONE);
    Array<Mag> errnrms(numVecs), normsB(numVecs), zeros(numVecs, MT::zero());
    Xhat.norm2(errnrms());
    B.norm2(normsB());
    Mag maxBnrm = *std::max_element( normsB.begin(), normsB.end() );
    TEST_COMPARE_FLOATING_ARRAYS( errnrms, zeros, maxBnrm );
  }
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
