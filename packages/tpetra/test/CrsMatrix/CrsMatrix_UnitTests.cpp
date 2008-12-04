#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Tuple.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
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
  using Teuchos::ArrayRCP;
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
  using Tpetra::INSERT;
  using Tpetra::Import;
  using std::string;
  using Teuchos::tuple;

  bool testMpi = true;
  double errorTolSlack = 1e+1;
  string filedir;

#define PRINT_VECTOR(v) \
   { \
     out << #v << ": "; \
     copy(v.begin(), v.end(), ostream_iterator<Ordinal>(out," ")); \
     out << endl; \
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
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
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
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
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
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, NonSquare, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const Ordinal M = 3;
    const Ordinal P = 5;
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
    const Ordinal indexBase = ZERO;
    const Ordinal numVecs  = 3;
    Map<Ordinal> rowmap(INVALID,M,indexBase,platform);
    Map<Ordinal> lclmap(P,indexBase,platform,true);
    // create the matrix
    CrsMatrix<Ordinal,Scalar> A(rowmap);
    for (int i=0; i<M; ++i) {
      for (int j=0; j<P; ++j) {
        A.submitEntry( M*myImageID+i, j, as<Scalar>(M*myImageID+i + j*M*N) );
      }
    }
    // call fillComplete()
    A.fillComplete(lclmap,rowmap);
    // build the input multivector X
    MV X(lclmap,numVecs);
    for (int i=0; i<P; ++i) {
      for (int j=0; j<numVecs; ++j) {
        X.replaceGlobalValue(i,j,as<Scalar>(i+j*P));
      }
    }
    // build the expected output multivector B
    MV Bexp(rowmap,numVecs), Bout(rowmap,numVecs);
    for (int i=myImageID*M; i<myImageID*M+M; ++i) {
      for (int j=0; j<numVecs; ++j) {
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, Transpose, Ordinal, Scalar )
  {
    // this is the same matrix as in test NonSquare, but we will apply the transpose
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const Ordinal M = 3;
    const Ordinal P = 5;
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
    const Ordinal indexBase = ZERO;
    const Ordinal numVecs  = 3;
    Map<Ordinal> rowmap(INVALID,M,indexBase,platform);
    Map<Ordinal> lclmap(P,indexBase,platform,true);
    // create the matrix
    CrsMatrix<Ordinal,Scalar> A(rowmap);
    for (int i=0; i<M; ++i) {
      for (int j=0; j<P; ++j) {
        A.submitEntry( M*myImageID+i, j, as<Scalar>(M*myImageID+i + j*M*N) );
      }
    }
    // call fillComplete()
    A.fillComplete(lclmap,rowmap);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, DomainRange, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
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
    const Ordinal indexBase = ZERO;
    const Ordinal numVecs  = 5;
    Map<Ordinal> rowmap(INVALID,tuple<Ordinal>(2*myImageID,2*myImageID+1),indexBase,platform);
    Map<Ordinal> rngmap(INVALID,tuple<Ordinal>(myImageID,numImages+myImageID),indexBase,platform);
    // create the tridiagonal matrix
    CrsMatrix<Ordinal,Scalar> tri(rowmap);
    Array<Scalar>  vals(3,ST::one());
    if (myImageID == 0) {
      Array<Ordinal> cols( tuple<Ordinal>(2*myImageID,2*myImageID+1,2*myImageID+2) );
      tri.submitEntries(2*myImageID  ,cols(0,2),vals(0,2));
      tri.submitEntries(2*myImageID+1,cols(0,3),vals(0,3));
    }
    else if (myImageID == numImages-1) {        
      Array<Ordinal> cols( tuple<Ordinal>(2*myImageID-1,2*myImageID,2*myImageID+1) );
      tri.submitEntries(2*myImageID  ,cols(0,3),vals(0,3));
      tri.submitEntries(2*myImageID+1,cols(1,2),vals(1,2));
    }
    else {
      Array<Ordinal> cols( tuple<Ordinal>(2*myImageID-1,2*myImageID,2*myImageID+1,2*myImageID+2) );
      tri.submitEntries(2*myImageID  ,cols(0,3),vals(0,3));
      tri.submitEntries(2*myImageID+1,cols(1,3),vals(0,3));
    }
    // call fillComplete(), specifying domain and range maps and requiring custom importer and exporter
    tri.fillComplete(rowmap,rngmap);
    // test the properties
    TEST_EQUALITY(tri.getNumGlobalNonzeros()  , 6*numImages-2);          
    TEST_EQUALITY(tri.getNumMyNonzeros()      , (myImageID > 0 && myImageID < numImages-1) ? 6 : 5);
    TEST_EQUALITY(tri.getNumGlobalRows()      , 2*numImages);
    TEST_EQUALITY(tri.getNumMyRows()          , 2);
    TEST_EQUALITY(tri.getNumMyCols()          , (myImageID > 0 && myImageID < numImages-1) ? 4 : 3);
    TEST_EQUALITY(tri.getNumGlobalDiagonals() , 2*numImages);
    TEST_EQUALITY(tri.getNumMyDiagonals()     , 2);
    TEST_EQUALITY(tri.getGlobalMaxNumEntries(), 3);
    TEST_EQUALITY(tri.getMyMaxNumEntries()    , 3);
    TEST_EQUALITY(tri.getIndexBase()          , 0);
    TEST_EQUALITY_CONST(tri.getRowMap().isSameAs(rowmap), true);
    TEST_EQUALITY_CONST(tri.getRangeMap().isSameAs(rngmap), true);
    TEST_EQUALITY_CONST(tri.getDomainMap().isSameAs(rowmap), true);
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
    tri.apply(mvin,mvout);
    mvout.update(-ST::one(),mvexp,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvout.norm2(norms());
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
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
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
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
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
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    const Ordinal indexBase = ZERO;
    Map<Ordinal> map(INVALID,ONE,indexBase,platform);
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
    TEST_EQUALITY_CONST(A.getRowMap().isSameAs(A.getDomainMap()), true);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, FullMatrixTriDiag, Ordinal, Scalar )
  {
    // do a FEM-type communication, then apply to a MultiVector containing the identity
    // this will check more difficult communication and test multivector apply
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal  ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 3) return;
    // create a Map
    const Ordinal indexBase = ZERO;
    Map<Ordinal> map(INVALID,ONE,indexBase,platform);

    // for debugging: Teuchos::VerboseObjectBase::setDefaultOStream(Teuchos::rcp(&out,false));
    
    /* create the following matrix:
    0  [2 1       ]   [2 1]
    1  [1 4 1     ]   [1 2] + [2 1]
    2  [  1 4 1   ]           [1 2] + 
    3  [    1     ] = 
       [       4 1]
   n-1 [       1 2]
    */
    Ordinal myNNZ;
    CrsMatrix<Ordinal,Scalar> A(map);
    MV mveye(map,numImages), mvans(map,numImages), mvres(map,numImages,false);
    if (myImageID != numImages-1) { // last image assigns none
      Array<Scalar> vals(3); vals[1] = as<Scalar>(1); vals[0] = vals[2] = as<Scalar>(2);
      Array<Ordinal> cols(2); cols[0] = myImageID; cols[1] = myImageID + 1;
      A.submitEntries(myImageID  ,cols(),vals(0,2));
      A.submitEntries(myImageID+1,cols(),vals(1,2));
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, FullMatrixComplex, Ordinal, Scalar )
  {
    // assumes that Scalar has a constructor of the form: Scalar(realpart,imagpart)
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    if (Teuchos::ScalarTraits<Scalar>::isOrdinal) return;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int myImageID = comm->getRank();

    int dim,dim2,nnz,info;
    double *dvals = NULL;
    int *colptr = NULL,
        *rowind = NULL;
    nnz = -1;
    string fn = filedir + "mhd1280b.cua";
    if (myImageID == 0) {
      info = readHB_newmat_double(fn.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
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
    if (info == 0 || nnz < 0) {
      success = false;
      out << "Error reading \"" << fn << "\"" << endl;
      return;
    }
    // create map: partition matrix equally among all procs
    Map<Ordinal> map_shared(dim,0,platform), map_AllOnRoot(dim,(myImageID==0?dim:0),0,platform);
    CrsMatrix<Ordinal,Scalar> A_crs(map_shared);
    // create a multivector with the entire matrix on Root, we will export it to the other procs
    MV A_mv(map_shared,dim), A_mv_AllOnRoot(map_AllOnRoot,dim), mvres(map_shared,dim), mveye(map_shared,dim);
    Import<Ordinal> AllFromRoot(map_AllOnRoot,map_shared);
    if (myImageID == 0) {
      // Root fills the CrsMatrix and the MV A_mv_AllOnRoot
      // HB format is compressed column. CrsMatrix is compressed row. Convert.
      const double *dptr = dvals;
      const int *rptr = rowind;
      for (int c=0; c<dim; ++c) {
        for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
          A_crs.submitEntry(*rptr-1,c,Scalar(dptr[0],dptr[1]));
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
    for (Ordinal j=0; j<map_shared.getNumMyEntries(); ++j) {
      Ordinal gid = map_shared.getGlobalIndex(j);
      mveye.replaceGlobalValue(gid,gid,ST::one());
    }
    // test the properties
    TEST_EQUALITY(A_crs.getNumGlobalNonzeros()   , nnz);
    TEST_EQUALITY(A_crs.getNumGlobalRows()       , dim);
    TEST_EQUALITY_CONST(A_crs.getIndexBase()     , ZERO);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, PowerComplex, Ordinal, Scalar )
  {
    // assumes that Scalar has a constructor of the form: Scalar(realpart,imagpart)
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    if (Teuchos::ScalarTraits<Scalar>::isOrdinal) return;
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int myImageID = comm->getRank();

    int dim,dim2,nnz,info;
    double *dvals = NULL;
    int *colptr = NULL,
        *rowind = NULL;
    nnz = -1;
    string fn = filedir + "mhd1280b.cua";
    if (myImageID == 0) {
      info = readHB_newmat_double(fn.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
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
    if (info == 0 || nnz < 0) {
      success = false;
      out << "Error reading \"" << fn << "\"" << endl;
      return;
    }
    // create map: partition matrix equally among all procs
    Map<Ordinal> map_shared(dim,0,platform);
    CrsMatrix<Ordinal,Scalar> A_crs(map_shared);
    if (myImageID == 0) {
      // Root fills the CrsMatrix and the MV A_mv_AllOnRoot
      // HB format is compressed column. CrsMatrix is compressed row. Convert.
      const double *dptr = dvals;
      const int *rptr = rowind;
      for (int c=0; c<dim; ++c) {
        for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
          A_crs.submitEntry(*rptr-1,c,Scalar(dptr[0],dptr[1]));
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, FullMatrix, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    if (Teuchos::ScalarTraits<Scalar>::isOrdinal) return;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int myImageID = comm->getRank();

    int dim,dim2,nnz,info;
    double *dvals = NULL;
    int *colptr = NULL,
        *rowind = NULL;
    nnz = -1;
    string fn = filedir + "west0067.rua";
    if (myImageID == 0) {
      info = readHB_newmat_double(fn.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
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
    if (info == 0 || nnz < 0) {
      success = false;
      out << "Error reading \"" << fn << "\"" << endl;
      return;
    }
    // create map: partition matrix equally among all procs
    Map<Ordinal> map_shared(dim,0,platform), map_AllOnRoot(dim,(myImageID==0?dim:0),0,platform);
    CrsMatrix<Ordinal,Scalar> A_crs(map_shared);
    // create a multivector with the entire matrix on Root, we will export it to the other procs
    MV A_mv(map_shared,dim), A_mv_AllOnRoot(map_AllOnRoot,dim), mvres(map_shared,dim), mveye(map_shared,dim);
    Import<Ordinal> AllFromRoot(map_AllOnRoot,map_shared);
    if (myImageID == 0) {
      // Root fills the CrsMatrix and the MV A_mv_AllOnRoot
      // HB format is compressed column. CrsMatrix is compressed row. Convert.
      const double *dptr = dvals;
      const int *rptr = rowind;
      for (int c=0; c<dim; ++c) {
        for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
          A_crs.submitEntry(*rptr-1,c,as<Scalar>(*dptr));
          A_mv_AllOnRoot.replaceGlobalValue(*rptr-1,c,as<Scalar>(*dptr));
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
    for (Ordinal j=0; j<map_shared.getNumMyEntries(); ++j) {
      Ordinal gid = map_shared.getGlobalIndex(j);
      mveye.replaceGlobalValue(gid,gid,ST::one());
    }
    // test the properties
    TEST_EQUALITY(A_crs.getNumGlobalNonzeros()   , nnz);
    TEST_EQUALITY(A_crs.getNumGlobalRows()       , dim);
    TEST_EQUALITY_CONST(A_crs.getIndexBase()     , ZERO);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, BadGID, Ordinal, Scalar )
  {
    // what happens when we call CrsMatrix::submitEntry() for a row that isn't on the Map?
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
    {
      // create the matrix
      CrsMatrix<Ordinal,Scalar> A(map);
      // add an entry off the map: row too high
      // this will only be off the map for the last node, for the others it will induce communication
      A.submitEntry(map.getMaxGlobalIndex()+1,map.getIndexBase(),ST::one());
      TEST_THROW(A.fillComplete(), std::runtime_error);
    }
    {
      // create the matrix
      CrsMatrix<Ordinal,Scalar> A(map);
      // add an entry off the map: row too high
      // this will only be off the map for the last node, for the others there is nothing
      if (myImageID == numImages-1) {
        A.submitEntry(map.getMaxAllGlobalIndex()+1,map.getIndexBase(),ST::one());
      }
      TEST_THROW(A.fillComplete(), std::runtime_error);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, ZeroMatrix, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Ordinal,Scalar> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    const Ordinal numVecs  = 5;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
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

#ifdef HAVE_TPETRA_TRIUTILS
# define TRIUTILS_USING_TESTS(ORDINAL,SCALAR) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, FullMatrix, ORDINAL, SCALAR )
#else
# define TRIUTILS_USING_TESTS(ORDINAL,SCALAR)
#endif

#ifdef HAVE_TPETRA_TRIUTILS
# define COMPLEX_TRIUTILS_USING_TESTS(ORDINAL,SCALAR) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, FullMatrixComplex, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, PowerComplex, ORDINAL, SCALAR )
#else
# define COMPLEX_TRIUTILS_USING_TESTS(ORDINAL,SCALAR)
#endif

#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)\
     typedef std::complex<float> ComplexFloat; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexFloat) \
     COMPLEX_TRIUTILS_USING_TESTS(ORDINAL, ComplexFloat)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)\
     typedef std::complex<double> ComplexDouble; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexDouble) \
     COMPLEX_TRIUTILS_USING_TESTS(ORDINAL, ComplexDouble)
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
#endif

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, TheEyeOfTruth, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, TheEyeOfTruthDistAlloc, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, ZeroMatrix   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, BadCalls     , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, SimpleEigTest, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, BadGID       , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, FullMatrixTriDiag, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, DomainRange, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, NonSquare, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, Transpose, ORDINAL, SCALAR ) \
      TRIUTILS_USING_TESTS(ORDINAL, SCALAR)

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double)
     UNIT_TEST_GROUP_ORDINAL(int)
# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, char)   \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int)    \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)  \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)  \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
     UNIT_TEST_GROUP_ORDINAL(int)

     typedef long int LongInt;
     UNIT_TEST_GROUP_ORDINAL(LongInt)
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt;
        UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
