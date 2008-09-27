#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"

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

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Tpetra::MultiVector;
  using std::endl;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Tpetra::CrsMatrix;

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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, BadCalls, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal NEGONE = ZERO - OrdinalTraits<Ordinal>::one();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    // create a random multivector
    MV mv1(map,1), mv2(map,2), mv3(map,3);
    // create the zero matrix
    CrsMatrix<Ordinal,Scalar> zero(map);
    TEST_THROW(zero.apply(mv1,mv1)            , std::runtime_error);
    zero.fillComplete();
    TEST_THROW(zero.submitEntry(0,0,ST::one()), std::runtime_error);
    TEST_THROW(zero.apply(mv2,mv1)            , std::runtime_error);
    TEST_THROW(zero.apply(mv2,mv3)            , std::runtime_error);
    // FINISH: add more tests here
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, TheEyeOfTruth, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal NEGONE = ZERO - OrdinalTraits<Ordinal>::one();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    const Ordinal numVecs  = 5;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    // create a random multivector
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.random();
    // create the identity matrix
    Ordinal base = numLocal*myImageID;
    CrsMatrix<Ordinal,Scalar> eye(map);
    for (int i=0; i<numLocal; ++i) {
      eye.submitEntry(base+i,base+i,ST::one());
    }
    eye.fillComplete();
    // test the properties
    TEST_EQUALITY(eye.getNumGlobalNonzeros()  , numImages*numLocal);
    TEST_EQUALITY(eye.getNumMyNonzeros()      , numLocal);
    TEST_EQUALITY(eye.getNumGlobalRows()      , numImages*numLocal);
    TEST_EQUALITY(eye.getNumMyRows()          , numLocal);
    TEST_EQUALITY(eye.getNumMyCols()          , numLocal);
    TEST_EQUALITY(eye.getNumGlobalDiagonals() , numImages*numLocal);
    TEST_EQUALITY(eye.getNumMyDiagonals()     , numLocal);
    TEST_EQUALITY(eye.getGlobalMaxNumEntries(), 1);
    TEST_EQUALITY(eye.getMyMaxNumEntries()    , 1);
    TEST_EQUALITY(eye.getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye.getRowMap().isSameAs(eye.getColMap())   , true);
    TEST_EQUALITY_CONST(eye.getRowMap().isSameAs(eye.getDomainMap()), true);
    TEST_EQUALITY_CONST(eye.getRowMap().isSameAs(eye.getRangeMap()) , true);
    // test the action
    mvres.random();
    eye.apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, TheEyeOfTruthDistAlloc, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal NEGONE = ZERO - OrdinalTraits<Ordinal>::one();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    const Ordinal numVecs  = 5;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    // create a random multivector
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.random();
    // create the identity matrix
    CrsMatrix<Ordinal,Scalar> eye(map);
    if (myImageID == 0) {
      for (int i=0; i<map.getNumGlobalEntries(); ++i) {
        eye.submitEntry(i,i,ST::one());
      }
    }
    eye.fillComplete();
    // test the properties
    TEST_EQUALITY(eye.getNumGlobalNonzeros()  , numImages*numLocal);
    TEST_EQUALITY(eye.getNumMyNonzeros()      , numLocal);
    TEST_EQUALITY(eye.getNumGlobalRows()      , numImages*numLocal);
    TEST_EQUALITY(eye.getNumMyRows()          , numLocal);
    TEST_EQUALITY(eye.getNumMyCols()          , numLocal);
    TEST_EQUALITY(eye.getNumGlobalDiagonals() , numImages*numLocal);
    TEST_EQUALITY(eye.getNumMyDiagonals()     , numLocal);
    TEST_EQUALITY(eye.getGlobalMaxNumEntries(), 1);
    TEST_EQUALITY(eye.getMyMaxNumEntries()    , 1);
    TEST_EQUALITY(eye.getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye.getRowMap().isSameAs(eye.getColMap())   , true);
    TEST_EQUALITY_CONST(eye.getRowMap().isSameAs(eye.getDomainMap()), true);
    TEST_EQUALITY_CONST(eye.getRowMap().isSameAs(eye.getRangeMap()) , true);
    // test the action
    mvres.random();
    eye.apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, SimpleEigTest, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal  ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    const Ordinal indexBase = ZERO;
    Map<Ordinal> map(NEGONE,ONE,indexBase,platform);
    // create a multivector ones(n,1)
    MV ones(map,ONE,false), threes(map,ONE,false);
    ones.putScalar(ST::one());
    /* create the following matrix:
       [2 1           ]
       [1 1 1         ]
       [  1 1 1       ]
       [   . . .      ]
       [     . . .    ]
       [       . . .  ]
       [         1 1 1]
       [           1 2]
     this matrix has an eigenvalue lambda=3, with eigenvector v = [1 ... 1]
    */
    Ordinal myNNZ;
    CrsMatrix<Ordinal,Scalar> A(map);
    if (myImageID == 0) {
      myNNZ = 2;
      Array<Scalar> vals(2); vals[0] = as<Scalar>(2.0)*ST::one(); vals[1]= ST::one();
      Array<Ordinal> cols(2); cols[0] = myImageID; cols[1] = myImageID+1;
      A.submitEntries(myImageID,cols(),vals());
    }
    else if (myImageID == numImages-1) {
      myNNZ = 2;
      Array<Scalar> vals(2); vals[0] = ST::one(); vals[1]= as<Scalar>(2.0)*ST::one();
      Array<Ordinal> cols(2); cols[0] = myImageID-1; cols[1] = myImageID;
      A.submitEntries(myImageID,cols(),vals());
    }
    else {
      myNNZ = 3;
      Array<Scalar> vals(3,ST::one());
      Array<Ordinal> cols(3); cols[0] = myImageID-1; cols[1] = myImageID; cols[2] = myImageID+1;
      A.submitEntries(myImageID,cols(),vals());
    }
    A.fillComplete();
    // test the properties
    TEST_EQUALITY(A.getNumGlobalNonzeros()   , 3*numImages-2);
    TEST_EQUALITY(A.getNumMyNonzeros()       , myNNZ);
    TEST_EQUALITY(A.getNumGlobalRows()       , numImages);
    TEST_EQUALITY_CONST(A.getNumMyRows()     , ONE);
    TEST_EQUALITY(A.getNumMyCols()           , myNNZ);
    TEST_EQUALITY(A.getNumGlobalDiagonals()  , numImages);
    TEST_EQUALITY_CONST(A.getNumMyDiagonals(), ONE);
    TEST_EQUALITY(A.getGlobalMaxNumEntries() , (numImages > 2 ? 3 : 2));
    TEST_EQUALITY(A.getMyMaxNumEntries()     , myNNZ);
    TEST_EQUALITY_CONST(A.getIndexBase()     , ZERO);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getColMap())   , false);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getDomainMap()), false);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getRangeMap()) , true);
    // test the action
    threes.random();
    A.apply(ones,threes);
    // now, threes should be 3*ones
    threes.update(as<Scalar>(-3.0)*ST::one(),ones,ST::one());
    Array<Mag> norms(1), zeros(1,MT::zero());
    threes.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, LessSimpleEigTest, Ordinal, Scalar )
  {
    // just like SimpleEigTest, but with a more interesting pattern and non-unique values
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal  ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;
    const Scalar  SONE = ST::one();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    const Ordinal indexBase = ZERO;
    Map<Ordinal> map(NEGONE,ONE,indexBase,platform);
    /* create the following matrix:
       [n n-1               ] [  1/n  ]     [2]
       [n n-1 n-2           ] [1/(n-1)]     [3]
       [  n-1 n-2 n-3       ] [1/(n-2)]     [3]
       [     .   .   .      ] [  ...  ]   = [.]
       [        .   .   .   ] [  1/3  ]     [3]
       [           3   2   1] [  1/2  ]     [3]
       [               2   1] [   1   ]     [2]
     this matrix has an eigenvalue lambda=3, with eigenvector v = [1 ... 1]
    */
    Ordinal myNNZ;
    CrsMatrix<Ordinal,Scalar> A(map);
    if (myImageID == 0) {
      myNNZ = 2;
      Array<Scalar> vals(2); vals[0] = as<Scalar>(numImages); vals[1] = vals[0] - SONE;
      Array<Ordinal> cols(2); cols[0] = myImageID; cols[1] = myImageID+1;
      A.submitEntries(myImageID,cols(),vals());
    }
    else if (myImageID == numImages-1) {
      myNNZ = 2;
      Array<Scalar> vals(2); vals[0] = as<Scalar>(2.0); vals[1] = SONE;
      Array<Ordinal> cols(2); cols[0] = myImageID-1; cols[1] = myImageID;
      A.submitEntries(myImageID,cols(),vals());
    }
    else {
      myNNZ = 3;
      Array<Scalar> vals(3); vals[0] = as<Scalar>(numImages - myImageID) + SONE; 
      vals[1] = vals[0] - SONE; vals[2] = vals[1] - SONE;
      Array<Ordinal> cols(3); cols[0] = myImageID-1; cols[1] = myImageID; cols[2] = myImageID+1;
      A.submitEntries(myImageID,cols(),vals());
    }
    A.fillComplete();
    // test the properties
    TEST_EQUALITY(A.getNumGlobalNonzeros()   , 3*numImages-2);
    TEST_EQUALITY(A.getNumMyNonzeros()       , myNNZ);
    TEST_EQUALITY(A.getNumGlobalRows()       , numImages);
    TEST_EQUALITY_CONST(A.getNumMyRows()     , ONE);
    TEST_EQUALITY(A.getNumMyCols()           , myNNZ);
    TEST_EQUALITY(A.getNumGlobalDiagonals()  , numImages);
    TEST_EQUALITY_CONST(A.getNumMyDiagonals(), ONE);
    TEST_EQUALITY(A.getGlobalMaxNumEntries() , (numImages > 2 ? 3 : 2));
    TEST_EQUALITY(A.getMyMaxNumEntries()     , myNNZ);
    TEST_EQUALITY_CONST(A.getIndexBase()     , ZERO);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getColMap())   , false);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getDomainMap()), false);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getRangeMap()) , true);
    // create a multivector: [1/n 1/(n-1) ... 2 1]
    MV fracs(map,ONE,false), res(map,ONE,false), ans(map,ONE,false);
    fracs.replaceMyValue(0,0, ST::one()/as<Scalar>(numImages-myImageID));
    if (myImageID == 0 || myImageID == (numImages-1)) {
      ans.replaceMyValue(0,0, as<Scalar>(2.0));
    }
    else {
      ans.replaceMyValue(0,0, as<Scalar>(3.0));
    }
    // test the action
    res.random();
    A.apply(fracs,res);
    // now, should res == ans
    res.update(-ST::one(),ans,ST::one());
    Array<Mag> norms(1), zeros(1,MT::zero());
    res.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, FullMatrix, Ordinal, Scalar )
  {
    // do a FEM-type communication, then apply to a MultiVector containing the identity
    // this will check more difficult communication and test multivector apply
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal  ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 3) return;
    // create a Map
    const Ordinal indexBase = ZERO;
    Map<Ordinal> map(NEGONE,ONE,indexBase,platform);
    /* create the following matrix:
    0  [1 .5           ]   [1  .5]
    1  [.5 2 .5        ]   [.5  1] + [1  .5]
    2  [  .5  2 .5     ]             [.5  1] + 
    3  [     .5        ] = 
       [           2 .5]
   n-1 [          .5  1]
    */
    Ordinal myNNZ;
    CrsMatrix<Ordinal,Scalar> A(map);
    MV mveye(map,numImages), mvans(map,numImages), mvres(map,numImages,false);
    if (myImageID != numImages-1) { // last image assigns none
      Array<Scalar> vals(3); vals[1] = as<Scalar>(0.5); vals[0] = vals[2] = ST::one();
      Array<Ordinal> cols(2); cols[0] = myImageID; cols[1] = myImageID + 1;
      A.submitEntries(myImageID  ,cols(),vals(0,2));
      A.submitEntries(myImageID+1,cols(),vals(1,2));
    }
    // divine myNNZ and build multivector with matrix
    mveye.replaceMyValue(0,myImageID,ST::one());
    out << "mveye: " << endl; mveye.printValues(out);
    if (myImageID == 0) {
      myNNZ = 2;
      mvans.replaceMyValue(0,0,as<Scalar>(1.0));
      mvans.replaceMyValue(0,1,as<Scalar>(0.5));
    }
    else if (myImageID == numImages-1) {
      myNNZ = 2;
      mvans.replaceMyValue(0,numImages-2,as<Scalar>(0.5));
      mvans.replaceMyValue(0,numImages-1,as<Scalar>(1.0));
    }
    else {
      myNNZ = 3;
      mvans.replaceMyValue(0,myImageID-1,as<Scalar>(0.5));
      mvans.replaceMyValue(0,myImageID  ,as<Scalar>(2.0));
      mvans.replaceMyValue(0,myImageID+1,as<Scalar>(0.5));
    }
    A.fillComplete();
    out << A << endl;
    out << "ColMap: " << A.getColMap() << endl;
    // test the properties
    TEST_EQUALITY(A.getNumGlobalNonzeros()   , 3*numImages-2);
    TEST_EQUALITY(A.getNumMyNonzeros()       , myNNZ);
    TEST_EQUALITY(A.getNumGlobalRows()       , numImages);
    TEST_EQUALITY_CONST(A.getNumMyRows()     , ONE);
    TEST_EQUALITY(A.getNumMyCols()           , myNNZ);
    TEST_EQUALITY(A.getNumGlobalDiagonals()  , numImages);
    TEST_EQUALITY_CONST(A.getNumMyDiagonals(), ONE);
    TEST_EQUALITY(A.getGlobalMaxNumEntries() , 3);
    TEST_EQUALITY(A.getMyMaxNumEntries()     , myNNZ);
    TEST_EQUALITY_CONST(A.getIndexBase()     , ZERO);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getColMap())   , false);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getDomainMap()), false);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getRangeMap()) , true);
    // test the action
    A.apply(mveye,mvres);
    out << "mvres: " << endl; mvres.printValues(out);
    mvres.update(-ST::one(),mvans,ST::one());
    Array<Mag> norms(numImages), zeros(numImages,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, ZeroMatrix, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal NEGONE = ZERO - OrdinalTraits<Ordinal>::one();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    const Ordinal numVecs  = 5;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    // create a random multivector
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.random();
    // create the zero matrix
    CrsMatrix<Ordinal,Scalar> zero(map);
    zero.fillComplete();
    mvres.random();
    zero.apply(mvrand,mvres);
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  // 
  // INSTANTIATIONS
  //

#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)\
     typedef std::complex<float> ComplexFloat; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexFloat)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)\
     typedef std::complex<double> ComplexDouble; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexDouble)
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
#endif

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, TheEyeOfTruth, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, TheEyeOfTruthDistAlloc, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, ZeroMatrix   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, BadCalls     , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, SimpleEigTest, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, LessSimpleEigTest, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, FullMatrix   , ORDINAL, SCALAR )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double)
     UNIT_TEST_GROUP_ORDINAL(int)
# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, char) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
     UNIT_TEST_GROUP_ORDINAL(int)

     typedef short int ShortInt;
     UNIT_TEST_GROUP_ORDINAL(ShortInt)
     typedef long int LongInt;
     UNIT_TEST_GROUP_ORDINAL(LongInt)
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt;
        UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
