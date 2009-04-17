#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#ifdef HAVE_TPETRA_TRIUTILS
#include <iohb.h>
#endif

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
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Tpetra::MultiVector;
  using std::endl;
  using std::swap;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::RowMatrix;
  using Tpetra::INSERT;
  using Tpetra::Import;
  using std::string;
  using Teuchos::tuple;

  bool testMpi = true;
  double errorTolSlack = 1e+1;
  string filedir;

#define ARRAYVIEW_TO_ARRAY(Type, arr, av) \
  { \
    ArrayView<Type> av2 = av; \
    arr.resize(av2.size()); \
    arr.assign(av2.begin(),av2.end()); \
  }

#define STD_TESTS(matrix) \
  { \
    RCP<const Comm<int> > STCOMM = matrix.getComm(); \
    ArrayView<const GO> STMYGIDS = matrix.getRowMap().getMyGlobalEntries(); \
    Teuchos_Ordinal STMAX = 0; \
    for (Teuchos_Ordinal STR=0; STR<matrix.numLocalRows(); ++STR) { \
      TEST_EQUALITY( matrix.numEntriesForMyRow(STR), matrix.numEntriesForGlobalRow( STMYGIDS[STR] ) ); \
      STMAX = std::max( STMAX, matrix.numEntriesForMyRow(STR) ); \
    } \
    TEST_EQUALITY( matrix.myMaxNumRowEntries(), STMAX ); \
    Teuchos_Ordinal STGMAX; \
    reduceAll( *STCOMM, Teuchos::REDUCE_MAX, STMAX, &STGMAX ); \
    TEST_EQUALITY( matrix.globalMaxNumRowEntries(), STGMAX ); \
  }


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
    RCP<Platform<double> > plat;
    if (testMpi) {
      plat = DefaultPlatform<double>::getPlatform();
    }
    else {
      plat = rcp(new Tpetra::SerialPlatform<double>());
    }
    return plat->getComm();
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, BadCalls, LO, GO, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Teuchos_Ordinal numLocal = 10;
    Map<LO,GO> map(INVALID,numLocal,0,comm);
    // create a random multivector
    MV mv1(map,1), mv2(map,2), mv3(map,3);
    // create the zero matrix
    RCP<RowMatrix<Scalar,LO,GO> > zero;
    {
      RCP<CrsMatrix<Scalar,LO,GO> > zero_crs = rcp( new CrsMatrix<Scalar,LO,GO>(map,0) );
      TEST_THROW(zero_crs->apply(mv1,mv1)            , std::runtime_error);
#   if defined(HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS)
      TEST_THROW(zero_crs->insertGlobalValues(map.getMinGlobalIndex(),tuple<GO>(0),tuple<Scalar>(ST::one())), std::runtime_error);
#   endif
      TEST_EQUALITY_CONST( zero_crs->isStaticGraph(), false );
      zero_crs->fillComplete();
      TEST_THROW(zero_crs->insertGlobalValues(0,tuple<GO>(0),tuple<Scalar>(ST::one())), std::runtime_error); // submit after fill
      zero = zero_crs;
    }
    TEST_THROW(zero->apply(mv2,mv1)            , std::runtime_error); // MVs have different number of vectors
    TEST_THROW(zero->apply(mv2,mv3)            , std::runtime_error); // MVs have different number of vectors
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, WithGraph, LO, GO, Scalar )
  {
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = size(*comm);
    // create a Map
    const Teuchos_Ordinal numLocal = 10;
    Map<LO,GO> map(INVALID,numLocal,0,comm);
    CrsGraph<LO,GO> graph(map,3,true);
    for (GO r=map.getMinGlobalIndex(); r<=map.getMaxGlobalIndex(); ++r) 
    {
      if (r == map.getMinAllGlobalIndex()) {
        graph.insertGlobalIndices(r,tuple(r,r+1));
      }
      else if (r == map.getMaxAllGlobalIndex()) {
        graph.insertGlobalIndices(r,tuple(r-1,r));
      }
      else {
        graph.insertGlobalIndices(r,tuple(r-1,r,r+1));
      }
    }
    graph.fillComplete();
    // create a matrix using the graph
    CrsMatrix<Scalar,LO,GO> matrix(graph);
    TEST_EQUALITY_CONST( matrix.isStaticGraph(), true );
    // insert throws exception: not allowed with static graph
    TEST_THROW( matrix.insertGlobalValues(map.getMinGlobalIndex(),tuple<GO>(map.getMinGlobalIndex()),tuple<Scalar>(ST::one())), std::runtime_error );
    // suminto and replace are allowed
    for (GO r=map.getMinGlobalIndex(); r<=map.getMaxGlobalIndex(); ++r) 
    {
      if (r == map.getMinAllGlobalIndex()) {
        matrix.replaceGlobalValues(r, tuple(r,r+1), tuple(ST::one(),ST::one()) );
      }
      else if (r == map.getMaxAllGlobalIndex()) {
        matrix.replaceGlobalValues(r, tuple(r-1,r), tuple(ST::one(),ST::one()) );
      }
      else {
        matrix.replaceGlobalValues(r, tuple(r-1,r,r+1), tuple(ST::one(),ST::one(),ST::one()) );
      }
    }
    for (GO r=map.getMinGlobalIndex(); r<=map.getMaxGlobalIndex(); ++r) 
    {
      matrix.sumIntoGlobalValues(r, tuple(r), tuple(ST::one()) );
    }
    matrix.fillComplete();
    TEST_EQUALITY( matrix.numMyDiagonals(), numLocal );
    TEST_EQUALITY( matrix.numGlobalDiagonals(), numImages*numLocal );
    TEST_EQUALITY( matrix.numGlobalEntries(), 3*numImages*numLocal - 2 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, ExceedStaticAlloc, LO, GO, Scalar )
  {
    // test that an exception is thrown when we exceed statically allocated memory
    typedef ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = size(*comm);
    // create a Map
    const Teuchos_Ordinal numLocal = 10;
    Map<LO,GO> map(INVALID,numLocal,0,comm);
    {
      CrsMatrix<Scalar,LO,GO> matrix(map,1,true);
      // room for one on each row
      for (GO r=map.getMinGlobalIndex(); r<=map.getMaxGlobalIndex(); ++r) 
      {
        matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
      }
      // no room for any more
      GO r = map.getMinGlobalIndex();
      TEST_THROW( matrix.insertGlobalValues( r, tuple(r+1), tuple(ST::one()) ), std::runtime_error );
    }
    if (numImages > 1) {
      // add too many entries globally
      CrsMatrix<Scalar,LO,GO> matrix(map,1,true);
      // room for one on each row
      for (GO r=map.getMinGlobalIndex(); r<=map.getMaxGlobalIndex(); ++r) 
      {
        matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
      }
      // always room for non-locals
      GO r = map.getMaxGlobalIndex() + 1;
      if (r > map.getMaxAllGlobalIndex()) r = map.getMinAllGlobalIndex();
      TEST_NOTHROW( matrix.insertGlobalValues( r, tuple(r), tuple(ST::one()) ) );
      // after communicating non-locals, failure trying to add them
      TEST_THROW( matrix.globalAssemble(), std::runtime_error );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, MultipleFillCompletes, LO, GO, Scalar )
  {
    // test that an exception is thrown when we exceed statically allocated memory
    typedef ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = size(*comm);
    // create a Map
    const Teuchos_Ordinal numLocal = 1; // change to 10
    Map<LO,GO> map(INVALID,numLocal,0,comm);
    {
      // room for two on each row
      CrsMatrix<Scalar,LO,GO> matrix(map,2,true);
      for (GO r=map.getMinGlobalIndex(); r<=map.getMaxGlobalIndex(); ++r) 
      {
        TEST_NOTHROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())) );
        TEST_NOTHROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())) );
      }
      // no room for more
      {
        GO r = map.getMinGlobalIndex();
        TEST_THROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())), std::runtime_error );
      }
      // fill complete adds them
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      TEST_NOTHROW( matrix.fillComplete() );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      // now there is room for more
      for (LO r=0; r<numLocal; ++r) 
      {
        TEST_NOTHROW( matrix.insertMyValues(r,tuple(r),tuple(ST::one())) );
      }
      TEST_NOTHROW( matrix.fillComplete() );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      // test that matrix is 3*I
      TEST_EQUALITY( matrix.numGlobalDiagonals(), numLocal*numImages );
      TEST_EQUALITY( matrix.numMyDiagonals(), numLocal );
      TEST_EQUALITY( matrix.numGlobalEntries(), numLocal*numImages );
      TEST_EQUALITY( matrix.numMyEntries(), numLocal );
      for (LO r=0; r<numLocal; ++r) {
        ArrayView<const LO> inds;
        ArrayView<const Scalar> vals;
        TEST_NOTHROW( matrix.extractMyRowConstView(r,inds,vals) );
        TEST_COMPARE_ARRAYS( inds, tuple<LO>(r) );
        TEST_COMPARE_ARRAYS( vals, tuple<Scalar>(as<Scalar>(3.0)) );
      }
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, CopiesAndViews, LO, GO, Scalar )
  {
    // test that an exception is thrown when we exceed statically allocated memory
    typedef ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = size(*comm);
    const int myImageID = rank(*comm);
    if (numImages < 2) return;
    // create a Map, one row per processor
    const Teuchos_Ordinal indexBase = 0;
    const Teuchos_Ordinal numLocal = 1;
    Map<LO,GO> rmap(INVALID,numLocal,indexBase,comm);
    GO myrowind = rmap.getGlobalIndex(0);
    // specify the column map to control ordering
    // construct tridiagonal graph
    Array<GO> ginds;
    Array<LO> linds;
    if (myImageID==0) {
      ARRAYVIEW_TO_ARRAY( GO, ginds, tuple(myrowind,myrowind+1) );
      ARRAYVIEW_TO_ARRAY( LO, linds, tuple(0,1) );
    }
    else if (myImageID==numImages-1) {
      ARRAYVIEW_TO_ARRAY( GO, ginds , tuple(myrowind-1,myrowind) );
      ARRAYVIEW_TO_ARRAY( LO, linds , tuple(0,1) );
    }
    else {
      ARRAYVIEW_TO_ARRAY( GO, ginds , tuple(myrowind-1,myrowind,myrowind+1) );
      ARRAYVIEW_TO_ARRAY( LO, linds , tuple(0,1,2) );
    }
    Array<Scalar> vals(ginds.size(),ST::one());
    Map<LO,GO> cmap(INVALID,ginds(),0,comm);
    for (int T=0; T<2; ++T) {
      bool StaticProfile = (T & 1) == 1;
      // bool OptimizeStorage = (T & 2) == 2;
      CrsMatrix<Scalar,LO,GO> matrix(rmap,cmap, ginds.size(),StaticProfile);   // only allocate as much room as necessary
      Array<GO> GCopy(4); Array<LO> LCopy(4); Array<Scalar> SCopy(4);
      ArrayView<const GO> CGView; ArrayView<const LO> CLView; ArrayView<const Scalar> CSView;
      ArrayView<GO> GView; ArrayView<LO> LView; ArrayView<Scalar> SView;
      Teuchos_Ordinal numentries;
      // at this point, the graph has not allocated data as global or local, so we can do views/copies for either local or global
      matrix.extractMyRowCopy(0,LCopy,SCopy,numentries);
      matrix.extractMyRowView(0,LView,SView);
      matrix.extractMyRowConstView(0,CLView,CSView);
      matrix.extractGlobalRowCopy(myrowind,GCopy,SCopy,numentries);
      matrix.extractGlobalRowView(myrowind,GView,SView);
      matrix.extractGlobalRowConstView(myrowind,CGView,CSView);
      // use multiple inserts: this illustrated an overwrite bug for column-map-specified graphs
      for (Teuchos_Ordinal j=0; j<(Teuchos_Ordinal)ginds.size(); ++j) {
        matrix.insertGlobalValues(myrowind,ginds(j,1),tuple(ST::one()));
      }
      TEST_EQUALITY( matrix.numEntriesForMyRow(0), matrix.getCrsGraph().numAllocatedEntriesForMyRow(0) ); // test that we only allocated as much room as necessary
      // if static graph, insert one additional entry on my row and verify that an exception is thrown
      if (StaticProfile) {
        TEST_THROW( matrix.insertGlobalValues(myrowind,arrayView(&myrowind,1),tuple(ST::one())), std::runtime_error );
      }
      matrix.fillComplete();
      // check that inserting global entries throws (inserting local entries is still allowed)
      {
        Array<GO> zero(0); Array<Scalar> vzero(0);
        TEST_THROW( matrix.insertGlobalValues(0,zero,vzero), std::runtime_error );
      }
      // check for throws and no-throws/values
      TEST_THROW( matrix.extractGlobalRowView(myrowind,GView,SView), std::runtime_error );
      TEST_THROW( matrix.extractGlobalRowConstView(myrowind,CGView,CSView), std::runtime_error );
      TEST_THROW( matrix.extractMyRowCopy(    0       ,LCopy(0,1),SCopy(0,1),numentries), std::runtime_error );
      TEST_THROW( matrix.extractGlobalRowCopy(myrowind,GCopy(0,1),SCopy(0,1),numentries), std::runtime_error );
      //
      TEST_NOTHROW( matrix.extractMyRowView(0,LView,SView) );
      TEST_NOTHROW( matrix.extractMyRowConstView(0,CLView,CSView) );
      TEST_COMPARE_ARRAYS( LView, CLView );
      TEST_COMPARE_ARRAYS( SView, CSView );
      TEST_COMPARE_ARRAYS( CLView, linds );
      TEST_COMPARE_ARRAYS( CSView, vals  );
      //
      TEST_NOTHROW( matrix.extractMyRowCopy(0,LCopy,SCopy,numentries) );
      TEST_COMPARE_ARRAYS( LCopy(0,numentries), linds );
      TEST_COMPARE_ARRAYS( SCopy(0,numentries), vals  );
      //
      TEST_NOTHROW( matrix.extractGlobalRowCopy(myrowind,GCopy,SCopy,numentries) );
      TEST_COMPARE_ARRAYS( GCopy(0,numentries), ginds );
      TEST_COMPARE_ARRAYS( SCopy(0,numentries), vals  );
      // 
      STD_TESTS(matrix);
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, TheEyeOfTruth, LO, GO, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a Map
    const LO numLocal = 10;
    const Teuchos_Ordinal numVecs  = 5;
    Map<LO,GO> map(INVALID,numLocal,0,comm);
    // create a random multivector
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.random();
    // create the identity matrix
    GO base = numLocal*myImageID;
    RCP<RowMatrix<Scalar,LO,GO> > eye;
    {
      RCP<CrsMatrix<Scalar,LO,GO> > eye_crs = rcp(new CrsMatrix<Scalar,LO,GO>(map,1));
      for (int i=0; i<numLocal; ++i) {
        eye_crs->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<Scalar>(ST::one()));
      }
      TEST_EQUALITY_CONST( eye_crs->isStaticGraph(), false );
      eye_crs->fillComplete();
      eye = eye_crs;
    }
    // test the properties
    TEST_EQUALITY(eye->numGlobalEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye->numMyEntries()      , numLocal);
    TEST_EQUALITY(eye->numGlobalRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->numLocalRows()          , numLocal);
    TEST_EQUALITY(eye->numLocalCols()          , numLocal);
    TEST_EQUALITY(eye->numGlobalDiagonals() , numImages*numLocal);
    TEST_EQUALITY(eye->numMyDiagonals()     , numLocal);
    TEST_EQUALITY(eye->globalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->myMaxNumRowEntries()    , 1);
    TEST_EQUALITY(eye->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye->getRowMap().isSameAs(eye->getColMap())   , true);
    TEST_EQUALITY_CONST(eye->getRowMap().isSameAs(eye->getDomainMap()), true);
    TEST_EQUALITY_CONST(eye->getRowMap().isSameAs(eye->getRangeMap()) , true);
    // test the action
    mvres.random();
    eye->apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, NonSquare, LO, GO, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const GO M = 3;
    const GO P = 5;
    const int N = comm->getSize();
    const int myImageID = comm->getRank();
    // create Maps
    // matrix is M*N-by-P
    //                  col
    //            0        1                  P-1
    //    0  [0        MN              ... (P-1)MN     ]
    //    .  [...      ...                 ...         ]
    //    0  [M-1      MN+M-1              (P-1)MN+M-1 ]
    //p   1  [M        MN+M                            ]
    //r   .  [...      ...                             ] = [A_ij], where A_ij = i+jMN
    //o   1  [2M-1     MN+2M-1                         ]
    //c   .  [...                                      ]
    //   N-1 [(N-1)M   MN+(N-1)(M-1)                   ]    
    //    .  [...      ...                             ]
    //   N-1 [MN-1     MN+MN-1                         ]
    // 
    // row map, range map is [0,M-1] [M,2M-1] [2M,3M-1] ... [MN-M,MN-1]
    // domain map will be map for X (lclmap)
    // 
    // input multivector X is not distributed:
    //               
    //   X = [  0    P    ...  (numVecs-1)P ]
    //       [ ...  ....  ...       ...     ] = [X_ji], where X_ij = i+jP
    //       [ P-1  2P-1  ...   numVecs*P-1 ]
    //
    // the result of the non-transpose multiplication should be 
    //                              P-1
    // (A*X)_ij = sum_k A_ik X_kj = sum (i+kMN)(k+jP) = jiP^2 + (i+jMNP)(P^2-P)/2 + MNP(P-1)(2P-1)/6
    //                              k=0
    // 
    // 
    // 
    const Teuchos_Ordinal numVecs  = 3;
    Map<LO,GO> rowmap(INVALID,M,0,comm);
    Map<LO,GO> lclmap(P,0,comm,true);
    // create the matrix
    CrsMatrix<Scalar,LO,GO> A(rowmap,P);
    for (GO i=0; i<M; ++i) {
      for (GO j=0; j<P; ++j) {
        A.insertGlobalValues( M*myImageID+i, tuple<GO>(j), tuple<Scalar>(M*myImageID+i + j*M*N) );
      }
    }
    // call fillComplete()
    TEST_EQUALITY_CONST( A.isStaticGraph(), false );
    A.fillComplete(lclmap,rowmap);
    // build the input multivector X
    MV X(lclmap,numVecs);
    for (GO i=0; i<P; ++i) {
      for (GO j=0; j<numVecs; ++j) {
        X.replaceGlobalValue(i,j,as<Scalar>(i+j*P));
      }
    }
    // build the expected output multivector B
    MV Bexp(rowmap,numVecs), Bout(rowmap,numVecs);
    for (GO i=myImageID*M; i<myImageID*M+M; ++i) {
      for (GO j=0; j<numVecs; ++j) {
        Bexp.replaceGlobalValue(i,j,as<Scalar>(j*i*P*P + (i+j*M*N*P)*(P*P-P)/2 + M*N*P*(P-1)*(2*P-1)/6));
      }
    }
    // test the action
    Bout.random();
    A.apply(X,Bout);
    Bout.update(-ST::one(),Bexp,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    Bout.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, Transpose, LO, GO, Scalar )
  {
    // this is the same matrix as in test NonSquare, but we will apply the transpose
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const GO M = 3;
    const GO P = 5;
    const int N = comm->getSize();
    const int myImageID = comm->getRank();
    // create Maps
    // matrix is M*N-by-P
    //                  col
    //            0        1                  P-1
    //    0  [0        MN              ... (P-1)MN     ]
    //    .  [...      ...                 ...         ]
    //    0  [M-1      MN+M-1              (P-1)MN+M-1 ]
    //p   1  [M        MN+M                            ]
    //r   .  [...      ...                             ] = [A_ij], where A_ij = i+jMN
    //o   1  [2M-1     MN+2M-1                         ]
    //c   .  [...                                      ]
    //   N-1 [(N-1)M   MN+(N-1)(M-1)                   ]    
    //    .  [...      ...                             ]
    //   N-1 [MN-1     MN+MN-1                         ]
    // 
    // row map, range map is [0,M-1] [M,2M-1] [2M,3M-1] ... [MN-M,MN-1]
    // domain map will be a non-distributed map for a vector of length P
    // 
    // input multivector is 
    //  col      0            1         ...        numVecs-1
    //     0 [0         MN                    (numVecs-1)MN     ]
    // p   . [...       ...                   ...               ]
    // r   0 [M-1       MN+M-1                (numVecs-1)MN+M-1 ]
    // o   1 [M         MN+M                                    ]
    // c   . [...       ...                                     ] = [X_ij], where X_ij = i+jMN
    //     1 [2M-1      MN+2M-1                                 ]
    //     . [...       ...                                     ]
    //    N-1[(N-1)M    MN+(N-1)(M-1)                           ]
    //     . [...       ...                                     ]
    //    N-1[MN-1      MN+MN-1                                 ]
    //
    // output multivector is not-distributed
    // the result of the transpose multiplication should be 
    //              MN-1              MN-1
    // (A^T*X)_ij = sum_k A_ki X_kj = sum (k+iMN)(k+jMN) 
    //              k=0               k=0
    //   MN-1
    // = sum k(i+j)MN + ij(MN)(MN) + k^2 = (i+j)(MN)^2(MN-1)/2 + ij(MN)^3 + (MN)(MN-1)(2MN-1)/6
    //   k=0
    // 
    const Teuchos_Ordinal numVecs  = 3;
    Map<LO,GO> rowmap(INVALID,M,0,comm);
    Map<LO,GO> lclmap(P,0,comm,true);
    // create the matrix
    CrsMatrix<Scalar,LO,GO> A(rowmap,P);
    for (int i=0; i<M; ++i) {
      for (int j=0; j<P; ++j) {
        A.insertGlobalValues( M*myImageID+i, tuple<GO>(j), tuple<Scalar>(M*myImageID+i + j*M*N) );
      }
    }
    // call fillComplete()
    TEST_EQUALITY_CONST( A.isStaticGraph(), false );
    A.fillComplete(lclmap,rowmap);
    out << "A: " << endl << A << endl;
    // build the input multivector X
    MV X(rowmap,numVecs);
    for (int i=myImageID*M; i<myImageID*M+M; ++i) {
      for (int j=0; j<numVecs; ++j) {
        X.replaceGlobalValue(i,j,as<Scalar>( i + j*M*N ) );
      }
    }
    // build the expected output multivector B
    MV Bexp(lclmap,numVecs), Bout(lclmap,numVecs);
    for (int i=0; i<P; ++i) {
      for (int j=0; j<numVecs; ++j) {
        Bexp.replaceGlobalValue(i,j,as<Scalar>( (i+j)*(M*N)*(M*N)*(M*N-1)/2 + i*j*(M*N)*(M*N)*(M*N) + (M*N)*(M*N-1)*(2*M*N-1)/6 ));
      }
    }
    // test the action
    Bout.random();
    A.apply(X,Bout,CONJ_TRANS);

    Bout.update(-ST::one(),Bexp,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    Bout.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, DomainRange, LO, GO, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages == 1) return;
    // create Maps
    // matrix is:
    //  proc
    //    0  [1 1               ]
    //    0  [1 1 1             ]
    //    1  [  1 1 1           ]
    //    1  [    1 1           ]
    //         ................
    //   n-1 [             1 1 1]
    //   n-1 [               1 1]
    //                              proc
    // input multivector will be      0  [0    2N   4N   6N    8N  ]
    //                                0  [1    2N+1 4N+1 6N+1  8N+1]
    //                                1  [2    2N+2 4N+2 6N+2  8N+2]
    //                                1  [3    2N+3 4N+3 6N+3  8N+3]                  (i,j) = 2*N*j+i
    //                                   [                         ]
    //                               n-1 [2N-2 4N-2 6N-2 8N-2 10N-2]
    //                               n-1 [2N-1 4N-1 6N-1 8N-1 10N-1]
    //                
    //                              proc                                                 
    // output multivector will be     0  [1    4N+1  8N+1 12N+1 16N+1]                    i=0:         (i,j)+(i+1,j) = 2Nj+i+2Nj+i+1 = 4Nj+2i+1 = 4Nj+1
    //                                1  [3    6N+3 12N+3 18N+3 24N+3]                 i=2n-1: (i-1,j)+(i,j)         = 2Nj+i-1+2Nj+i = 4Nj+2i-1 = 4Nj+2(2N-1)-1 = 4N(j+1)-3
    //                                2  [6    6N+6 12N+6 18N+6 24N+6]                   else: (i-1,j)+(i,j)+(i+1,j) = 2Nj+i-1+2Nj+i+2Nj+i+1 = 6Nj+3i
    //                                   [                           ]
    //                               n-1 [                           ]
    //                                0  [                           ]
    //                                1  [                           ]
    //                                   [                           ]
    //                               n-1 [4N-3 8N-3 12N-3 16N-3 20N-3]
    //
    // row map is [0,1]   [2,3]     [4,5]     etc
    // col map is [0,1,2] [1,2,3,4] [3,4,5,6] etc     (assembled by CrsMatrix, we construct one only for comparison)
    // domain map will be equal to the row map
    // range  map will be [0,np] [1,np+1] [2,np+2]
    const Teuchos_Ordinal numVecs  = 5;
    Map<LO,GO> rowmap(INVALID,tuple<GO>(2*myImageID,2*myImageID+1),0,comm);
    Map<LO,GO> rngmap(INVALID,tuple<GO>(myImageID,numImages+myImageID),0,comm);
    RCP<RowMatrix<Scalar,LO,GO> > tri;
    {
      RCP<CrsMatrix<Scalar,LO,GO> > tri_crs = rcp(new CrsMatrix<Scalar,LO,GO>(rowmap,3) );
      Array<Scalar>  vals(3,ST::one());
      if (myImageID == 0) {
        Array<GO> cols( tuple<GO>(2*myImageID,2*myImageID+1,2*myImageID+2) );
        tri_crs->insertGlobalValues(2*myImageID  ,cols(0,2),vals(0,2));
        tri_crs->insertGlobalValues(2*myImageID+1,cols(0,3),vals(0,3));
      }
      else if (myImageID == numImages-1) {        
        Array<GO> cols( tuple<GO>(2*myImageID-1,2*myImageID,2*myImageID+1) );
        tri_crs->insertGlobalValues(2*myImageID  ,cols(0,3),vals(0,3));
        tri_crs->insertGlobalValues(2*myImageID+1,cols(1,2),vals(1,2));
      }
      else {
        Array<GO> cols( tuple<GO>(2*myImageID-1,2*myImageID,2*myImageID+1,2*myImageID+2) );
        tri_crs->insertGlobalValues(2*myImageID  ,cols(0,3),vals(0,3));
        tri_crs->insertGlobalValues(2*myImageID+1,cols(1,3),vals(0,3));
      }
      // call fillComplete(), specifying domain and range maps and requiring custom importer and exporter
      tri_crs->fillComplete(rowmap,rngmap);
      tri = tri_crs;
    }
    // test the properties
    TEST_EQUALITY(tri->numGlobalEntries()  , 6*numImages-2);          
    TEST_EQUALITY(tri->numMyEntries()      , (myImageID > 0 && myImageID < numImages-1) ? 6 : 5);
    TEST_EQUALITY(tri->numGlobalRows()      , 2*numImages);
    TEST_EQUALITY(tri->numLocalRows()          , 2);
    TEST_EQUALITY(tri->numLocalCols()          , (myImageID > 0 && myImageID < numImages-1) ? 4 : 3);
    TEST_EQUALITY(tri->numGlobalDiagonals() , 2*numImages);
    TEST_EQUALITY(tri->numMyDiagonals()     , 2);
    TEST_EQUALITY(tri->globalMaxNumRowEntries(), 3);
    TEST_EQUALITY(tri->myMaxNumRowEntries()    , 3);
    TEST_EQUALITY(tri->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(tri->getRowMap().isSameAs(rowmap), true);
    TEST_EQUALITY_CONST(tri->getRangeMap().isSameAs(rngmap), true);
    TEST_EQUALITY_CONST(tri->getDomainMap().isSameAs(rowmap), true);
    // build the input and corresponding output multivectors
    MV mvin(rowmap,numVecs), mvout(rngmap,numVecs), mvexp(rngmap,numVecs);
    for (int j=0; j<numVecs; ++j) {
      mvin.replaceMyValue(0,j,as<Scalar>(j*2*numImages+2*myImageID  )); // entry (2*myImageID  ,j)
      mvin.replaceMyValue(1,j,as<Scalar>(j*2*numImages+2*myImageID+1)); // entry (2*myImageID+1,j)
      // entry (myImageID,j)
      if (myImageID==0) {
        mvexp.replaceMyValue(0,j,as<Scalar>(4*numImages*j+1));
      }                                                                    
      else {                                                               
        mvexp.replaceMyValue(0,j,as<Scalar>(6*numImages*j+3*myImageID));
      }                                                                    
      // entry (numImages+myImageID,j)
      if (myImageID==numImages-1) {                                        
        mvexp.replaceMyValue(1,j,as<Scalar>(4*numImages*(j+1)-3));
      }
      else {
        mvexp.replaceMyValue(1,j,as<Scalar>(6*numImages*j+3*(numImages+myImageID)));
      }
    }
    // test the action
    mvout.random();
    tri->apply(mvin,mvout);
    mvout.update(-ST::one(),mvexp,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvout.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, TheEyeOfTruthDistAlloc, LO, GO, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a Map
    const LO numLocal = 10;
    const Teuchos_Ordinal numVecs  = 5;
    Map<LO,GO> map(INVALID,numLocal,0,comm);
    // create a random multivector
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.random();
    // create the identity matrix
    RCP<RowMatrix<Scalar,LO,GO> > eye;
    {
      RCP<CrsMatrix<Scalar,LO,GO> > eye_crs = rcp(new CrsMatrix<Scalar,LO,GO>(map,1) );
      if (myImageID == 0) {
        for (int i=0; i<map.getNumGlobalEntries(); ++i) {
          eye_crs->insertGlobalValues(i,tuple<GO>(i),tuple<Scalar>(ST::one()));
        }
      }
      eye_crs->fillComplete();
      eye = eye_crs;
    }
    // test the properties
    TEST_EQUALITY(eye->numGlobalEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye->numMyEntries()      , numLocal);
    TEST_EQUALITY(eye->numGlobalRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->numLocalRows()          , numLocal);
    TEST_EQUALITY(eye->numLocalCols()          , numLocal);
    TEST_EQUALITY(eye->numGlobalDiagonals() , numImages*numLocal);
    TEST_EQUALITY(eye->numMyDiagonals()     , numLocal);
    TEST_EQUALITY(eye->globalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->myMaxNumRowEntries()    , 1);
    TEST_EQUALITY(eye->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye->getRowMap().isSameAs(eye->getColMap())   , true);
    TEST_EQUALITY_CONST(eye->getRowMap().isSameAs(eye->getDomainMap()), true);
    TEST_EQUALITY_CONST(eye->getRowMap().isSameAs(eye->getRangeMap()) , true);
    // test the action
    mvres.random();
    eye->apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, SimpleEigTest, LO, GO, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const LO ONE = OrdinalTraits<LO>::one();
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    Map<LO,GO> map(INVALID,ONE,0,comm);
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
    Teuchos_Ordinal myNNZ;
    CrsMatrix<Scalar,LO,GO> A(map,3);
    if (myImageID == 0) {
      myNNZ = 2;
      Array<Scalar> vals(tuple<Scalar>(as<Scalar>(2)*ST::one(), ST::one()));
      Array<GO> cols(tuple<GO>(myImageID, myImageID+1));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    else if (myImageID == numImages-1) {
      myNNZ = 2;
      Array<Scalar> vals(tuple<Scalar>(ST::one(), as<Scalar>(2)*ST::one()));
      Array<GO> cols(tuple<GO>(myImageID-1,myImageID));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    else {
      myNNZ = 3;
      Array<Scalar> vals(3,ST::one());
      Array<GO> cols(tuple<GO>(myImageID-1, myImageID, myImageID+1));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    A.fillComplete();
    // test the properties
    TEST_EQUALITY(A.numGlobalEntries()   , 3*numImages-2);
    TEST_EQUALITY(A.numMyEntries()       , myNNZ);
    TEST_EQUALITY(A.numGlobalRows()       , numImages);
    TEST_EQUALITY_CONST(A.numLocalRows()     , ONE);
    TEST_EQUALITY(A.numLocalCols()           , myNNZ);
    TEST_EQUALITY(A.numGlobalDiagonals()  , numImages);
    TEST_EQUALITY_CONST(A.numMyDiagonals(), ONE);
    TEST_EQUALITY(A.globalMaxNumRowEntries() , (numImages > 2 ? 3 : 2));
    TEST_EQUALITY(A.myMaxNumRowEntries()     , myNNZ);
    TEST_EQUALITY_CONST(A.getIndexBase()     , 0);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getColMap())   , false);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getDomainMap()), true);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getRangeMap()) , true);
    // test the action
    threes.random();
    A.apply(ones,threes);
    // now, threes should be 3*ones
    threes.update(as<Scalar>(-3)*ST::one(),ones,ST::one());
    Array<Mag> norms(1), zeros(1,MT::zero());
    threes.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, FullMatrixTriDiag, LO, GO, Scalar )
  {
    // do a FEM-type communication, then apply to a MultiVector containing the identity
    // this will check more difficult communication and test multivector apply
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const LO ONE = OrdinalTraits<LO>::one();
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 3) return;
    // create a Map
    Map<LO,GO> map(INVALID,ONE,0,comm);

    // for debugging: Teuchos::VerboseObjectBase::setDefaultOStream(Teuchos::rcp(&out,false));
    
    /* create the following matrix:
    0  [2 1       ]   [2 1]
    1  [1 4 1     ]   [1 2] + [2 1]
    2  [  1 4 1   ]           [1 2] + 
    3  [    1     ] = 
       [       4 1]
   n-1 [       1 2]
    */
    Teuchos_Ordinal myNNZ;
    CrsMatrix<Scalar,LO,GO> A(map,4);
    MV mveye(map,numImages), mvans(map,numImages), mvres(map,numImages,false);
    if (myImageID != numImages-1) { // last image assigns none
      Array<Scalar> vals(tuple<Scalar>(as<Scalar>(2)*ST::one(),ST::one(),as<Scalar>(2)*ST::one()));
      Array<GO> cols(tuple<GO>(myImageID,myImageID + 1));
      A.insertGlobalValues(myImageID  ,cols(),vals(0,2));
      A.insertGlobalValues(myImageID+1,cols(),vals(1,2));
    }
    // divine myNNZ and build multivector with matrix
    mveye.replaceMyValue(0,myImageID,ST::one());
    if (myImageID == 0) {
      myNNZ = 2;
      mvans.replaceMyValue(0,0,as<Scalar>(2));
      mvans.replaceMyValue(0,1,as<Scalar>(1));
    }
    else if (myImageID == numImages-1) {
      myNNZ = 2;
      mvans.replaceMyValue(0,numImages-2,as<Scalar>(1));
      mvans.replaceMyValue(0,numImages-1,as<Scalar>(2));
    }
    else {
      myNNZ = 3;
      mvans.replaceMyValue(0,myImageID-1,as<Scalar>(1));
      mvans.replaceMyValue(0,myImageID  ,as<Scalar>(4));
      mvans.replaceMyValue(0,myImageID+1,as<Scalar>(1));
    }
    A.fillComplete();
    // test the properties
    TEST_EQUALITY(A.numGlobalEntries()   , 3*numImages-2);
    TEST_EQUALITY(A.numMyEntries()       , myNNZ);
    TEST_EQUALITY(A.numGlobalRows()       , numImages);
    TEST_EQUALITY_CONST(A.numLocalRows()     , ONE);
    TEST_EQUALITY(A.numLocalCols()           , myNNZ);
    TEST_EQUALITY(A.numGlobalDiagonals()  , numImages);
    TEST_EQUALITY_CONST(A.numMyDiagonals(), ONE);
    TEST_EQUALITY(A.globalMaxNumRowEntries() , 3);
    TEST_EQUALITY(A.myMaxNumRowEntries()     , myNNZ);
    TEST_EQUALITY_CONST(A.getIndexBase()     , 0);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getColMap())   , false);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getDomainMap()), true);
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getRangeMap()) , true);
    // test the action
    A.apply(mveye,mvres);
    mvres.update(-ST::one(),mvans,ST::one());
    Array<Mag> norms(numImages), zeros(numImages,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
#ifdef HAVE_TPETRA_TRIUTILS
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, FullMatrixComplex, LO, GO, Scalar )
  {
    // assumes that Scalar has a constructor of the form: Scalar(realpart,imagpart)
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    if (Teuchos::ScalarTraits<Scalar>::isOrdinal) return;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();

    int dim,dim2,nnz,info;
    int rnnzmax;
    double *dvals = NULL;
    int *colptr = NULL,
        *rowind = NULL;
    nnz = -1;
    string fn = filedir + "mhd1280b.cua";
    if (myImageID == 0) {
      info = readHB_newmat_double(fn.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
      // find maximum NNZ over all rows
      vector<int> rnnz(dim,0);
      for (int *ri=rowind; ri<rowind+nnz; ++ri) {
        ++rnnz[*ri-1];
      }
      rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
    }
    else {
      // address uninitialized data warnings
      dvals = NULL;
      colptr = NULL;
      rowind = NULL;
    }
    Teuchos::broadcast(*comm,0,&info);
    Teuchos::broadcast(*comm,0,&nnz);
    Teuchos::broadcast(*comm,0,&dim);
    Teuchos::broadcast(*comm,0,&rnnzmax);
    if (info == 0 || nnz < 0) {
      success = false;
      out << "Error reading \"" << fn << "\"" << endl;
      return;
    }
    // create map: partition matrix equally among all procs
    Map<LO,GO> map_shared(dim,0,comm), map_AllOnRoot(dim,(myImageID==0?dim:0),0,comm);
    CrsMatrix<Scalar,LO,GO> A_crs(map_shared,rnnzmax);
    // create a multivector with the entire matrix on Root, we will export it to the other procs
    MV A_mv(map_shared,dim), A_mv_AllOnRoot(map_AllOnRoot,dim), mvres(map_shared,dim), mveye(map_shared,dim);
    Import<LO,GO> AllFromRoot(map_AllOnRoot,map_shared);
    if (myImageID == 0) {
      // Root fills the CrsMatrix and the MV A_mv_AllOnRoot
      // HB format is compressed column. CrsMatrix is compressed row. Convert.
      const double *dptr = dvals;
      const int *rptr = rowind;
      for (int c=0; c<dim; ++c) {
        for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
          A_crs.insertGlobalValues(*rptr-1,tuple<GO>(c),tuple<Scalar>(Scalar(dptr[0],dptr[1])));
          A_mv_AllOnRoot.replaceGlobalValue(*rptr-1,c,Scalar(dptr[0],dptr[1]));
          ++rptr;
          dptr += 2;
        }
      }
    }
    // fillComplete() will distribute matrix entries from Root to all other procs
    A_crs.fillComplete();
    // doExport() will distribute MV entries from Root to all other procs
    A_mv.doImport(A_mv_AllOnRoot,AllFromRoot,INSERT);

    if (myImageID == 0) {
      // Clean up allocated memory.
      free( dvals );
      free( colptr );
      free( rowind );
    }

    // build identity MV
    for (GO j=0; j<map_shared.getNumMyEntries(); ++j) {
      GO gid = map_shared.getGlobalIndex(j);
      mveye.replaceGlobalValue(gid,gid,ST::one());
    }
    // test the properties
    TEST_EQUALITY(A_crs.numGlobalEntries()   , nnz);
    TEST_EQUALITY(A_crs.numGlobalRows()       , dim);
    TEST_EQUALITY_CONST(A_crs.getIndexBase()     , 0);
    TEST_EQUALITY_CONST(A_crs.getRowMap().isSameAs(A_crs.getRangeMap()) , true);
    // test the action
    A_crs.apply(mveye,mvres);
    mvres.update(-ST::one(),A_mv,ST::one());
    Array<Mag> norms(dim), zeros(dim,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }
#endif


  ////
#ifdef HAVE_TPETRA_TRIUTILS
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, PowerComplex, LO, GO, Scalar )
  {
    // assumes that Scalar has a constructor of the form: Scalar(realpart,imagpart)
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    if (Teuchos::ScalarTraits<Scalar>::isOrdinal) return;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();

    int dim,dim2,nnz,info;
    int rnnzmax;
    double *dvals = NULL;
    int *colptr = NULL,
        *rowind = NULL;
    nnz = -1;
    string fn = filedir + "mhd1280b.cua";
    if (myImageID == 0) {
      info = readHB_newmat_double(fn.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
      // find maximum NNZ over all rows
      vector<int> rnnz(dim,0);
      for (int *ri=rowind; ri<rowind+nnz; ++ri) {
        ++rnnz[*ri-1];
      }
      rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
    }
    else {
      // address uninitialized data warnings
      dvals = NULL;
      colptr = NULL;
      rowind = NULL;
    }
    Teuchos::broadcast(*comm,0,&info);
    Teuchos::broadcast(*comm,0,&nnz);
    Teuchos::broadcast(*comm,0,&dim);
    Teuchos::broadcast(*comm,0,&rnnzmax);
    if (info == 0 || nnz < 0) {
      success = false;
      out << "Error reading \"" << fn << "\"" << endl;
      return;
    }
    // create map: partition matrix equally among all procs
    Map<LO,GO> map_shared(dim,0,comm);
    CrsMatrix<Scalar,LO,GO> A_crs(map_shared,rnnzmax);
    if (myImageID == 0) {
      // Root fills the CrsMatrix and the MV A_mv_AllOnRoot
      // HB format is compressed column. CrsMatrix is compressed row. Convert.
      const double *dptr = dvals;
      const int *rptr = rowind;
      for (int c=0; c<dim; ++c) {
        for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
          A_crs.insertGlobalValues(*rptr-1,tuple<GO>(c),tuple<Scalar>(Scalar(dptr[0],dptr[1])));
          ++rptr;
          dptr += 2;
        }
      }
    }
    // fillComplete() will distribute matrix entries from Root to all other procs
    A_crs.fillComplete();
    if (myImageID == 0) {
      // Clean up allocated memory.
      free( dvals );
      free( colptr );
      free( rowind );
    }

    // simple power method
    RCP<MV> x = rcp(new MV(map_shared,1)), 
            r = rcp(new MV(map_shared,1));
    Scalar lam, lam_left;
    Mag nrm, nrm_left;
    x->putScalar(Scalar(1.0f,1.0f)); x->norm2(arrayView<Mag>(&nrm,1)); x->scale(MT::one()/nrm);
    for (int i=0; i<20; ++i) {
      A_crs.apply(*x,*r);                                         // r = A*x
      x->dot(*r,arrayView<Scalar>(&lam,1));                       // lambda = x'*r = x'*A*x
      x->update(ST::one(),*r,-lam);                               // x = r - x*lam = A*x - x*lam \doteq oldres
      r->norm2(arrayView<Mag>(&nrm,1)); r->scale(MT::one()/nrm);  // r = A*x / |A*x| \doteq newx
      swap(x,r);                                                  // x = newx; r = oldres
      r->norm2(arrayView<Mag>(&nrm,1));                           // nrm = |r| = |oldres|
      out << "i: " << i << "\t\tlambda: " << lam << "\t\t|r|: " << nrm << endl;
    }
    // check that the computed right eigenpair is also a left eigenpair (the matrix is Hermitian)
    A_crs.apply(*x,*r,CONJ_TRANS);
    x->dot(*r,arrayView<Scalar>(&lam_left,1));
    x->update(ST::one(),*r,-lam_left);  // x = A'*x - x*lam_left
    x->norm2(arrayView<Mag>(&nrm_left,1));
    out << "lam_left: " << lam_left << "\t\tnrm_left: " << nrm_left << endl;
    TEST_FLOATING_EQUALITY(lam, Scalar(70.322f,0.0f), as<Mag>(0.000001f));
    TEST_FLOATING_EQUALITY(lam_left, lam, as<Mag>(0.000001f));
    TEST_EQUALITY_CONST(nrm      < 0.0001f, true);
    TEST_EQUALITY_CONST(nrm_left < 0.0001f, true);
  }
#endif


  ////
#ifdef HAVE_TPETRA_TRIUTILS
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, FullMatrix, LO, GO, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    if (Teuchos::ScalarTraits<Scalar>::isOrdinal) return;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();

    int dim,dim2,nnz,info;
    int rnnzmax;
    double *dvals = NULL;
    int *colptr = NULL,
        *rowind = NULL;
    nnz = -1;
    string fn = filedir + "west0067.rua";
    if (myImageID == 0) {
      info = readHB_newmat_double(fn.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
      // find maximum NNZ over all rows
      vector<int> rnnz(dim,0);
      for (int *ri=rowind; ri<rowind+nnz; ++ri) {
        ++rnnz[*ri-1];
      }
      rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
    }
    else {
      // address uninitialized data warnings
      dvals = NULL;
      colptr = NULL;
      rowind = NULL;
    }
    Teuchos::broadcast(*comm,0,&info);
    Teuchos::broadcast(*comm,0,&nnz);
    Teuchos::broadcast(*comm,0,&dim);
    Teuchos::broadcast(*comm,0,&rnnzmax);
    if (info == 0 || nnz < 0) {
      success = false;
      out << "Error reading \"" << fn << "\"" << endl;
      return;
    }
    // create map: partition matrix equally among all procs
    Map<LO,GO> map_shared(dim,0,comm), map_AllOnRoot(dim,(myImageID==0?dim:0),0,comm);
    CrsMatrix<Scalar,LO,GO> A_crs(map_shared,rnnzmax);
    // create a multivector with the entire matrix on Root, we will export it to the other procs
    MV A_mv(map_shared,dim), A_mv_AllOnRoot(map_AllOnRoot,dim), mvres(map_shared,dim), mveye(map_shared,dim);
    Import<LO,GO> AllFromRoot(map_AllOnRoot,map_shared);
    if (myImageID == 0) {
      // Root fills the CrsMatrix and the MV A_mv_AllOnRoot
      // HB format is compressed column. CrsMatrix is compressed row. Convert.
      double *dptr = dvals;
      int *rptr = rowind;
      for (int c=0; c<dim; ++c) {
        for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
          Scalar s = as<Scalar>(*dptr);
          A_crs.insertGlobalValues(*rptr-1,tuple<GO>(c),tuple(s));
          A_mv_AllOnRoot.replaceGlobalValue(*rptr-1,c,s);
          ++rptr;
          ++dptr;
        }
      }
    }
    // fillComplete() will distribute matrix entries from Root to all other procs
    A_crs.fillComplete();
    // doExport() will distribute MV entries from Root to all other procs
    A_mv.doImport(A_mv_AllOnRoot,AllFromRoot,INSERT);

    if (myImageID == 0) {
      // Clean up allocated memory.
      free( dvals );
      free( colptr );
      free( rowind );
    }

    // build identity MV
    for (GO j=0; j<map_shared.getNumMyEntries(); ++j) {
      GO gid = map_shared.getGlobalIndex(j);
      mveye.replaceGlobalValue(gid,gid,ST::one());
    }
    // test the properties
    TEST_EQUALITY(A_crs.numGlobalEntries()   , nnz);
    TEST_EQUALITY(A_crs.numGlobalRows()       , dim);
    TEST_EQUALITY_CONST(A_crs.getIndexBase()     , 0);
    TEST_EQUALITY_CONST(A_crs.getRowMap().isSameAs(A_crs.getRangeMap()) , true);
    // test the action
    A_crs.apply(mveye,mvres);
    mvres.update(-ST::one(),A_mv,ST::one());
    Array<Mag> norms(dim), zeros(dim,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }
#endif


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, BadGID, LO, GO, Scalar )
  {
    // what happens when we call CrsMatrix::insertGlobalValues() for a row that isn't on the Map?
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    // create a Map
    const LO numLocal = 10;
    Map<LO,GO> map(INVALID,numLocal,0,comm);
    {
      // create the matrix
      CrsMatrix<Scalar,LO,GO> A(map,1);
      // add an entry off the map: row too high
      // this will only be off the map for the last node, for the others it will induce communication
      A.insertGlobalValues(map.getMaxGlobalIndex()+1,tuple<GO>(map.getIndexBase()),tuple<Scalar>(ST::one()));
      TEST_THROW(A.fillComplete(), std::runtime_error);
    }
    {
      // create the matrix
      CrsMatrix<Scalar,LO,GO> A(map,1);
      // add an entry off the map: row too high
      // this will only be off the map for the last node, for the others there is nothing
      if (myImageID == numImages-1) {
        A.insertGlobalValues(map.getMaxAllGlobalIndex()+1,tuple<GO>(map.getIndexBase()),tuple<Scalar>(ST::one()));
      }
      TEST_THROW(A.fillComplete(), std::runtime_error);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, ZeroMatrix, LO, GO, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const LO numLocal = 10;
    const Teuchos_Ordinal numVecs  = 5;
    Map<LO,GO> map(INVALID,numLocal,0,comm);
    // create a random multivector
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.random();
    // create the zero matrix
    CrsMatrix<Scalar,LO,GO> zero(map,0);
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

#ifdef HAVE_TPETRA_TRIUTILS
# define TRIUTILS_USING_TESTS(LO, GO,SCALAR) \
      /*TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, FullMatrix, LO, GO, SCALAR ) */
#else
# define TRIUTILS_USING_TESTS(LO, GO,SCALAR)
#endif

#ifdef HAVE_TPETRA_TRIUTILS
# define COMPLEX_TRIUTILS_USING_TESTS(LO, GO,SCALAR) \
      /* TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, FullMatrixComplex, LO, GO, SCALAR ) */ \
      /* TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, PowerComplex, LO, GO, SCALAR ) */
#else
# define COMPLEX_TRIUTILS_USING_TESTS(LO, GO,SCALAR)
#endif

#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)\
     typedef std::complex<float> ComplexFloat; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexFloat) \
     COMPLEX_TRIUTILS_USING_TESTS(LO, GO, ComplexFloat)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)\
     typedef std::complex<double> ComplexDouble; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexDouble) \
     COMPLEX_TRIUTILS_USING_TESTS(LO, GO, ComplexDouble)
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#endif

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, TheEyeOfTruth, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, TheEyeOfTruthDistAlloc, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, ZeroMatrix   , LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, BadCalls     , LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, SimpleEigTest, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, BadGID       , LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, FullMatrixTriDiag, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, DomainRange, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, NonSquare, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, Transpose, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, WithGraph, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, ExceedStaticAlloc, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, MultipleFillCompletes, LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, CopiesAndViews, LO, GO, SCALAR ) \
      TRIUTILS_USING_TESTS( LO, GO, SCALAR )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
    UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, double)
     UNIT_TEST_GROUP_ORDINAL(int)
# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, float)  \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)

     UNIT_TEST_GROUP_ORDINAL(int)

     typedef long int LongInt;
     UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongInt )
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt;
        UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongLongInt )
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
