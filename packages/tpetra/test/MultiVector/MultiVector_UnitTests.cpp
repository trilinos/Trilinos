#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

// FINISH: add test for MultiVector with a node containing zero local entries
// FINISH: add tests for local MultiVectors 

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

  using std::vector;
  using std::sort;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::SerialDenseMatrix;
  using Teuchos::as;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::arrayView;
  using Teuchos::tuple;
  using Teuchos::rcp;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::MultiVector;
  using std::endl;
  using std::copy;
  using std::ostream_iterator;
  using std::string;

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

  Node& getDefaultNode()
  {
    return DefaultPlatform::getDefaultPlatform().getNode();
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, basic, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal numLocal = 10;
    const Teuchos_Ordinal numVecs  = 5;
    Map<Ordinal> map(INVALID,numLocal,0,comm);
    MV mvec(node,map,numVecs,true);
    TEST_EQUALITY( mvec.getNumVectors(), numVecs );
    TEST_EQUALITY( mvec.getMyLength(), numLocal );
    TEST_EQUALITY( mvec.getGlobalLength(), numImages*numLocal );
    // we zeroed it out in the constructor; all norms should be zero
    Array<Magnitude> norms(numVecs), zeros(numVecs);
    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Magnitude>::zero());
    mvec.norm2(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    mvec.norm1(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    mvec.normInf(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    // print it
    out << mvec << endl;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstNumVecs, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal numLocal = 10;
    Map<Ordinal> map(INVALID,numLocal,0,comm);
    TEST_THROW(MV mvec(node,map,0), std::invalid_argument);
    TEST_THROW(MV mvec(node,map,-1), std::invalid_argument);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstLDA, Ordinal, Scalar )
  {
    // numlocal > LDA
    // ergo, the arrayview doesn't contain enough data to specify the entries
    // also, if bounds checking is enabled, check that bad bounds are caught
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Teuchos_Ordinal numVecs = 2;
    // multivector has two vectors, each proc having two values per vector
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    // we need 4 scalars to specify values on each proc
    Array<Scalar> values(4);
#ifdef HAVE_TPETRA_DEBUG
    // too small an ArrayView (less than 4 values) is met with an exception, if debugging is on
    TEST_THROW(MV mvec(node,map,values(0,3),2,numVecs), std::runtime_error);
    // it could also be too small for the given LDA: 
    TEST_THROW(MV mvec(node,map,values(),2+1,numVecs), std::runtime_error);
    // too small for number of entries in a Vector
    TEST_THROW(V   vec(node,map,values(0,1)), std::runtime_error);
#endif
    // LDA < numLocal throws an exception anytime
    TEST_THROW(MV mvec(node,map,values(0,4),1,numVecs), std::runtime_error);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, LabeledObject, LO, GO, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    Node &node = getDefaultNode();
    // create Map
    Map<LO,GO> map(INVALID,3,0,comm);
    // test labeling
    const string lbl("mvecA");
    MV mvecA(node,map,2);
    string desc1 = mvecA.description();
    if (myImageID==0) out << desc1 << endl;
    mvecA.setObjectLabel(lbl);
    string desc2 = mvecA.description();
    if (myImageID==0) out << desc2 << endl;
    if (myImageID==0) {
      TEST_EQUALITY( mvecA.getObjectLabel(), lbl );
    }
    // test describing at different verbosity levels
    if (myImageID==0) out << "Describing with verbosity VERB_DEFAULT..." << endl;
    mvecA.describe(out);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_NONE..." << endl;
    mvecA.describe(out,VERB_NONE);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_LOW..." << endl;
    mvecA.describe(out,VERB_LOW);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_MEDIUM..." << endl;
    mvecA.describe(out,VERB_MEDIUM);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_HIGH..." << endl;
    mvecA.describe(out,VERB_HIGH);
    comm->barrier();
    comm->barrier();
    if (myImageID==0) out << "Describing with verbosity VERB_EXTREME..." << endl;
    mvecA.describe(out,VERB_EXTREME);
    comm->barrier();
    comm->barrier();
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadMultiply, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    const Scalar S1 = ScalarTraits<Scalar>::one(),
                 S0 = ScalarTraits<Scalar>::zero();
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    {
      // create local Maps
      Map<Ordinal> map3l(3,0,comm,true),
                   map2l(2,0,comm,true);
      MV mvecA(node,map3l,2),
         mvecB(node,map2l,3),
         mvecD(node,map2l,2);
      // failures, 8 combinations:
      // [NTC],[NTC]: A,B don't match
      // [NTC],[NTC]: C doesn't match A,B
      TEST_THROW( mvecD.multiply(NO_TRANS  ,NO_TRANS  ,S1,mvecA,mvecA,S0), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(NO_TRANS  ,CONJ_TRANS,S1,mvecA,mvecB,S0), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(CONJ_TRANS,NO_TRANS  ,S1,mvecB,mvecA,S0), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(CONJ_TRANS,CONJ_TRANS,S1,mvecB,mvecB,S0), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(NO_TRANS  ,NO_TRANS  ,S1,mvecA,mvecB,S0), std::runtime_error);   // 2x2: 3x2 x 2x3
      TEST_THROW( mvecD.multiply(NO_TRANS  ,CONJ_TRANS,S1,mvecA,mvecA,S0), std::runtime_error);   // 2x2: 3x2 x 2x3
      TEST_THROW( mvecD.multiply(CONJ_TRANS,NO_TRANS  ,S1,mvecB,mvecB,S0), std::runtime_error);   // 2x2: 3x2 x 2x3
      TEST_THROW( mvecD.multiply(CONJ_TRANS,CONJ_TRANS,S1,mvecB,mvecA,S0), std::runtime_error);   // 2x2: 3x2 x 2x3
    }
    // case 2: C(local) = A^T(distr) * B  (distr)  : one of these
    {
      Map<Ordinal> map3n(INVALID,3,0,comm),
                   map2n(INVALID,2,0,comm),
                   map2l(2,0,comm,true),
                   map3l(3,0,comm,true);
      MV mv3nx2(node,map3n,2),
         mv2nx2(node,map2n,2),
         mv2lx2(node,map2l,2),
         mv2lx3(node,map2l,3),
         mv3lx2(node,map3l,2),
         mv3lx3(node,map3l,3);
      // non-matching input lengths
      TEST_THROW( mv2lx2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv2nx2,S0), std::runtime_error);   // (2 x 3n) x (2n x 2) not compat
      TEST_THROW( mv2lx2.multiply(CONJ_TRANS,NO_TRANS,S1,mv2nx2,mv3nx2,S0), std::runtime_error);   // (2 x 2n) x (3n x 2) not compat
      // non-matching output size
      TEST_THROW( mv3lx3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx2,S0), std::runtime_error);   // (2 x 3n) x (3n x 2) doesn't fit 3x3
      TEST_THROW( mv3lx2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx2,S0), std::runtime_error);   // (2 x 3n) x (3n x 2) doesn't fit 3x2
      TEST_THROW( mv2lx3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx2,S0), std::runtime_error);   // (2 x 3n) x (3n x 2) doesn't fit 2x3
    }
    // case 3: C(distr) = A  (distr) * B^X(local)  : two of these
    {
      Map<Ordinal> map3n(INVALID,3,0,comm),
                   map2n(INVALID,2,0,comm),
                   map2l(2,0,comm,true),
                   map3l(3,0,comm,true);
      MV mv3nx2(node,map3n,2),
         mv2nx2(node,map2n,2),
         mv2x3(node,map2l,3),
         mv3x2(node,map3l,2);
      // non-matching input lengths
      TEST_THROW( mv3nx2.multiply(NO_TRANS,CONJ_TRANS,S1,mv3nx2,mv2x3,S0), std::runtime_error);   // (3n x 2) x (3 x 2) (trans) not compat
      TEST_THROW( mv3nx2.multiply(NO_TRANS,NO_TRANS  ,S1,mv3nx2,mv3x2,S0), std::runtime_error);   // (3n x 2) x (3 x 2) (nontrans) not compat
      // non-matching output sizes
      TEST_THROW( mv3nx2.multiply(NO_TRANS,CONJ_TRANS,S1,mv3nx2,mv3x2,S0), std::runtime_error);   // (3n x 2) x (2 x 3) doesn't fit 3nx2
      TEST_THROW( mv3nx2.multiply(NO_TRANS,NO_TRANS  ,S1,mv3nx2,mv2x3,S0), std::runtime_error);   // (3n x 2) x (2 x 3) doesn't fit 3nx2
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, Multiply, Ordinal, Scalar )
  {
    using Teuchos::View;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    Node &node = getDefaultNode();
    // create a Map
    Map<Ordinal> map3n(INVALID,3,0,comm),
                 map2n(INVALID,2,0,comm),
                 lmap3(3,0,comm,true),
                 lmap2(2,0,comm,true);
    const Scalar S1 = ScalarTraits<Scalar>::one(),
                 S0 = ScalarTraits<Scalar>::zero();
    const Mag    M0 = ScalarTraits<Mag>::zero();
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    // deterministic input/output
    {
      MV mv3x2l(node,lmap3,2),
         mv2x3l(node,lmap2,3),
         mv2x2l(node,lmap2,2),
         mv3x3l(node,lmap3,3);
      // fill multivectors with ones
      mv3x2l.putScalar(ScalarTraits<Scalar>::one());
      mv2x3l.putScalar(ScalarTraits<Scalar>::one());
      // fill expected answers Array
      Teuchos::Array<Scalar> check2(4,3); // each entry (of four) is the product [1 1 1]*[1 1 1]' = 3
      Teuchos::Array<Scalar> check3(9,2); // each entry (of nine) is the product [1 1]*[1 1]' = 2
      // test
      ArrayRCP<const Scalar> tmpView;
      mv3x3l.multiply(NO_TRANS  ,NO_TRANS  ,S1,mv3x2l,mv2x3l,S0);
      tmpView = mv3x3l.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView(0,9),check3,M0);
      mv2x2l.multiply(NO_TRANS  ,CONJ_TRANS,S1,mv2x3l,mv2x3l,S0);
      tmpView = mv2x2l.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView(0,4),check2,M0);
      mv2x2l.multiply(CONJ_TRANS,NO_TRANS  ,S1,mv3x2l,mv3x2l,S0);
      tmpView = mv2x2l.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView(0,4),check2,M0);
      mv3x3l.multiply(CONJ_TRANS,CONJ_TRANS,S1,mv2x3l,mv3x2l,S0);
      tmpView = mv3x3l.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView(0,9),check3,M0);
    }
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    // random input/output
    {
      Array<Scalar>     tmvCopy1(6), tmvCopy2(6);
      ArrayView<Scalar> sdmView(Teuchos::null);
      MV tmv3x2(node,lmap3,2),
         tmv2x3(node,lmap2,3),
         tmv2x2(node,lmap2,2),
         tmv3x3(node,lmap3,3);
      // fill multivectors with random, get copy of contents
      tmv3x2.random();  tmv3x2.get1dCopy(tmvCopy1(),3); 
      tmv2x3.random();  tmv2x3.get1dCopy(tmvCopy2(),2);
      // point SerialDenseMatrices at copies
      SerialDenseMatrix<int,Scalar> sdm3x2(View,tmvCopy1.getRawPtr(),3,3,2);
      SerialDenseMatrix<int,Scalar> sdm2x3(View,tmvCopy2.getRawPtr(),2,2,3);
      // space for answers
      SerialDenseMatrix<int,Scalar> sdm2x2(2,2), sdm3x3(3,3);
      // test: perform local Tpetra::MultiVector multiply and Teuchos::SerialDenseMatrix multiply, then check that answers are equivalent
      ArrayRCP<const Scalar> tmpView;
      {
        tmv3x3.multiply(NO_TRANS,NO_TRANS,S1,tmv3x2,tmv2x3,S0);
        sdm3x3.multiply(NO_TRANS,NO_TRANS,S1,sdm3x2,sdm2x3,S0);
        tmpView = tmv3x3.get1dView(); sdmView = arrayView(sdm3x3.values(),sdm3x3.numRows()*sdm3x3.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmpView,sdmView,ScalarTraits<Mag>::eps() * 10.);
      }
      {
        tmv2x2.multiply(NO_TRANS,CONJ_TRANS,S1,tmv2x3,tmv2x3,S0);
        sdm2x2.multiply(NO_TRANS,CONJ_TRANS,S1,sdm2x3,sdm2x3,S0);
        tmpView = tmv2x2.get1dView(); sdmView = arrayView(sdm2x2.values(),sdm2x2.numRows()*sdm2x2.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmpView,sdmView,ScalarTraits<Mag>::eps() * 10.);
      }
      {
        tmv2x2.multiply(CONJ_TRANS,NO_TRANS,S1,tmv3x2,tmv3x2,S0);
        sdm2x2.multiply(CONJ_TRANS,NO_TRANS,S1,sdm3x2,sdm3x2,S0);
        tmpView = tmv2x2.get1dView(); sdmView = arrayView(sdm2x2.values(),sdm2x2.numRows()*sdm2x2.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmpView,sdmView,ScalarTraits<Mag>::eps() * 10.);
      }
      {
        tmv3x3.multiply(CONJ_TRANS,CONJ_TRANS,S1,tmv2x3,tmv3x2,S0);
        sdm3x3.multiply(CONJ_TRANS,CONJ_TRANS,S1,sdm2x3,sdm3x2,S0);
        tmpView = tmv3x3.get1dView(); sdmView = arrayView(sdm3x3.values(),sdm3x3.numRows()*sdm3x3.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmpView,sdmView,ScalarTraits<Mag>::eps() * 10.);
      }
    }
    // case 2: C(local) = A^T(distr) * B  (distr)  : one of these
    {
      MV mv3nx2(node,map3n,2),
         mv3nx3(node,map3n,3),
         // locals
         mv2x2(node,lmap2,2),
         mv2x3(node,lmap2,3),
         mv3x2(node,lmap3,2),
         mv3x3(node,lmap3,3);
      // fill multivectors with ones
      mv3nx3.putScalar(ScalarTraits<Scalar>::one());
      mv3nx2.putScalar(ScalarTraits<Scalar>::one());
      // fill expected answers Array
      ArrayRCP<const Scalar> tmpView;
      Teuchos::Array<Scalar> check(9,3*numImages);
      // test
      mv2x2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx2,S0); 
      tmpView = mv2x2.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);
      mv2x3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx3,S0);
      tmpView = mv2x3.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);
      mv3x2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx3,mv3nx2,S0);
      tmpView = mv3x2.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);
      mv3x3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx3,mv3nx3,S0);
      tmpView = mv3x3.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);
    }
    // case 3: C(distr) = A  (distr) * B^X(local)  : two of these
    {
      // FINISH
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstAA, Ordinal, Scalar )
  {
    // constructor takes ArrayView<ArrayView<Scalar> A, NumVectors
    // A.size() == NumVectors
    // A[i].size() >= MyLength
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    // multivector has two vectors, each proc having two values per vector
    Map<Ordinal> map2(INVALID,numLocal  ,indexBase,comm),
                 map3(INVALID,numLocal+1,indexBase,comm);
    // we need 4 scalars to specify values on each proc
    Array<Scalar> values(4);
    Array<ArrayView<const Scalar> > arrOfarr(2,ArrayView<const Scalar>(Teuchos::null));
    Array<ArrayView<const Scalar> > emptyArr;
    arrOfarr[0] = values(0,2);
    arrOfarr[1] = values(2,2);
    // arrOfarr.size() == 0
    TEST_THROW(MV mvec(node,map2,emptyArr(),0), std::runtime_error);
#ifdef HAVE_TPETRA_DEBUG
    // individual ArrayViews could be too small
    TEST_THROW(MV mvec(node,map3,arrOfarr(),2), std::runtime_error);
#endif
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadDot, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    Map<Ordinal> map1(INVALID,1,indexBase,comm),
                 map2(INVALID,2,indexBase,comm);
    {
      MV mv12(node,map1,1),
         mv21(node,map2,1),
         mv22(node,map2,2);
      Array<Scalar> dots(2);
      // incompatible maps
      TEST_THROW(mv12.dot(mv21,dots()),std::runtime_error);
      // incompatible numvecs
      TEST_THROW(mv22.dot(mv21,dots()),std::runtime_error);
      // too small output array
#ifdef TEUCHOS_DEBUG
      TEST_THROW(mv22.dot(mv22,dots(0,1)),std::runtime_error);
#endif
    }
    {
      V v1(node,map1),
        v2(node,map2);
      // incompatible maps
      TEST_THROW(v1.dot(v2),std::runtime_error);
      TEST_THROW(v2.dot(v1),std::runtime_error);
      // wrong size output array through MultiVector interface
      Array<Scalar> dots(2);
#ifdef TEUCHOS_DEBUG
      TEST_THROW(v1.dot(v2,dots()),std::runtime_error);
      TEST_THROW(v2.dot(v1,dots()),std::runtime_error);
#endif
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, OrthoDot, Ordinal, Scalar )
  {
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Scalar S0 = ScalarTraits<Scalar>::zero();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Teuchos_Ordinal numVectors = 2;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    Array<Scalar> values(5);
    // values = {0, 1, 0, 1, 0}
    // values(0,4) = {0, 1, 0, 1} = [0 0]
    //                            = [1 1]
    // values(1,4) = {1, 0, 1, 0} = [1 1]
    //                            = [0 0]
    // these should be numerical orthogonal even in finite arithmetic
    values[0] = as<Scalar>(1);
    values[1] = as<Scalar>(0);
    values[2] = as<Scalar>(1);
    values[3] = as<Scalar>(0);
    values[4] = as<Scalar>(1);
    MV mvec1(node,map,values(0,4),2,numVectors),
       mvec2(node,map,values(1,4),2,numVectors);
    Array<Scalar> dots1(numVectors), dots2(numVectors), zeros(numVectors);
    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Scalar>::zero());
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
    TEST_EQUALITY_CONST( mvec1.getVector(0)->dot(*mvec2.getVector(0)), S0);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, ZeroScaleUpdate, Ordinal, Scalar )
  {
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Teuchos_Ordinal numVectors = 2;
    const Teuchos_Ordinal LDA = 2;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    Array<Scalar> values(6);
    // values = {1, 1, 2, 2, 4, 4}
    // values(0,4) = {1, 1, 2, 2} = [1 2]
    //                            = [1 2]
    // values(2,6) = {2, 2, 4, 4} = [2 4]
    //                            = [2 4]
    // a multivector A constructed from the first 
    // has values .5 of a multivector B constructed from the second
    // then 2*A - B = 0
    // we test both scale(), both update(), and norm()
    values[0] = as<Scalar>(1);
    values[1] = as<Scalar>(1);
    values[2] = as<Scalar>(2);
    values[3] = as<Scalar>(2);
    values[4] = as<Scalar>(4);
    values[5] = as<Scalar>(4);
    MV A(node,map,values(0,4),LDA,numVectors),
       B(node,map,values(2,4),LDA,numVectors);
    Array<Mag> norms(numVectors), zeros(numVectors);
    std::fill(zeros.begin(),zeros.end(),M0);
    //
    //      [.... ....]
    // A == [ones ones] 
    //      [.... ....]
    // 
    //      [.... ....]
    // B == [twos twos]
    //      [.... ....]
    //
    //   set A2 = A
    //   scale it by 2 in situ
    //   check that it equals B: subtraction in situ
    {
      MV A2(A);
      A2.scale(as<Scalar>(2));
      A2.update(as<Scalar>(-1),B,as<Scalar>(1));
      A2.norm2(norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }
    //   set A2 = A
    //   check that it equals B: scale,subtraction in situ
    {
      MV A2(A);
      A2.update(as<Scalar>(-1),B,as<Scalar>(2));
      A2.norm2(norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }
    //   set C random
    //   set it to zero by combination with A,B
    {
      MV C(node,map,numVectors);
      C.random();
      C.update(as<Scalar>(-1),B,as<Scalar>(2),A,as<Scalar>(0));
      C.norm2(norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }
    //   set C random
    //   scale it ex-situ
    //   check that it equals B: subtraction in situ
    {
      MV C(node,map,numVectors);
      C.scale(as<Scalar>(2),A);
      C.update(as<Scalar>(1),B,as<Scalar>(-1));
      C.norm2(norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, ZeroScaleUpdate, Ordinal, Scalar )
  {
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    Map<Ordinal> map(INVALID,2,0,comm);
    Array<Scalar> values(6);
    // values = {1, 1, 2, 2}
    // values(0,2) = {1, 1} = [1]
    //                      = [1]
    // values(2,2) = {2, 2} = [2]
    //                      = [2]
    // a vector A constructed from the first 
    // has values .5 of a vector B constructed from the second
    // thus 2*A - B = 0
    // we test both scale(), both update(), and norm()
    values[0] = as<Scalar>(1);
    values[1] = as<Scalar>(1);
    values[2] = as<Scalar>(2);
    values[3] = as<Scalar>(2);
    V A(node,map,values(0,2)),
      B(node,map,values(2,2));
    Mag norm;
    Array<Mag> norms(1);
    //
    //      [....]
    // A == [ones]
    //      [....]
    // 
    //      [....]
    // B == [twos]
    //      [....]
    //
    //   set A2 = A
    //   scale it by 2 in situ
    //   check that it equals B: subtraction in situ
    {
      V A2(A);
      A2.scale(as<Scalar>(2));
      A2.update(as<Scalar>(-1),B,as<Scalar>(1));
      norm = A2.norm2(); A2.norm2(norms());
      TEST_EQUALITY(norm,M0);
      TEST_EQUALITY(norm,norms[0]);
    }
    //   set A2 = A
    //   check that it equals B: scale,subtraction in situ
    {
      V A2(A);
      A2.update(as<Scalar>(-1),B,as<Scalar>(2));
      norm = A2.norm2(); A2.norm2(norms());
      TEST_EQUALITY(norm,M0);
      TEST_EQUALITY(norm,norms[0]);
    }
    //   set C random
    //   set it to zero by combination with A,B
    {
      V C(node,map);
      C.random();
      C.update(as<Scalar>(-1),B,as<Scalar>(2),A,as<Scalar>(0));
      norm = C.norm2(); C.norm2(norms());
      TEST_EQUALITY(norm,M0);
      TEST_EQUALITY(norm,norms[0]);
    }
    //   set C random
    //   scale it ex-situ
    //   check that it equals B: subtraction in situ
    {
      V C(node,map);
      C.random();
      C.scale(as<Scalar>(2),A);
      C.update(as<Scalar>(1),B,as<Scalar>(-1));
      norm = C.norm2(); C.norm2(norms());
      TEST_EQUALITY(norm,M0);
      TEST_EQUALITY(norm,norms[0]);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CopyConst, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Teuchos_Ordinal numVectors = 2;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    // create random MV
    MV morig(node,map,numVectors);
    morig.random();
    // copy it
    MV mcopy1(morig), mcopy2(morig);
    // verify that all three have identical values
    Array<Magnitude> norig(numVectors), ncopy1(numVectors), ncopy2(numVectors);
    morig.normInf(norig);
    mcopy1.normInf(ncopy1);
    mcopy2.normInf(ncopy2);
    TEST_COMPARE_FLOATING_ARRAYS(norig,ncopy1,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norig,ncopy2,M0);
    TEST_COMPARE_FLOATING_ARRAYS(ncopy1,ncopy2,M0);
    // modify all three
    morig.putScalar(as<Scalar>(0));
    mcopy1.putScalar(as<Scalar>(1));
    mcopy2.putScalar(as<Scalar>(2));
    // compute norms
    morig.normInf(norig);
    mcopy1.normInf(ncopy1);
    mcopy2.normInf(ncopy2);
    // check them
    bool local_success = true;
    for (Teuchos_Ordinal i=0; i<numVectors; ++i) {
      TEST_ARRAY_ELE_EQUALITY( norig,  i, as<Scalar>(0) );
      TEST_ARRAY_ELE_EQUALITY( ncopy1, i, as<Scalar>(1) );
      TEST_ARRAY_ELE_EQUALITY( ncopy2, i, as<Scalar>(2) );
    }
    success &= local_success;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, CopyConst, Ordinal, Scalar )
  {
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    Map<Ordinal> map(INVALID,2,0,comm);
    // create random MV
    V morig(node,map);
    morig.random();
    // copy it
    V mcopy1(morig), mcopy2(morig);
    // verify that all three have identical values
    Magnitude norig, ncopy1, ncopy2;
    norig = morig.normInf();
    ncopy1 = mcopy1.normInf();
    ncopy2 = mcopy2.normInf();
    TEST_EQUALITY(norig,ncopy1);
    TEST_EQUALITY(norig,ncopy2);
    TEST_EQUALITY(ncopy1,ncopy2);
    // modify all three
    morig.putScalar(as<Scalar>(0));
    mcopy1.putScalar(as<Scalar>(1));
    mcopy2.putScalar(as<Scalar>(2));
    // compute norms
    norig = morig.normInf();
    ncopy1 = mcopy1.normInf();
    ncopy2 = mcopy2.normInf();
    // check them
    TEST_EQUALITY(norig, as<Scalar>(0));
    TEST_EQUALITY(ncopy1,as<Scalar>(1));
    TEST_EQUALITY(ncopy2,as<Scalar>(2));
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, Indexing, Ordinal, Scalar )
  {
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    typedef ScalarTraits<Scalar>              SCT;
    typedef typename SCT::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    Map<Ordinal> map(INVALID,100,0,comm);
    // create two random Vector objects
    V v1(node,map), v2(node,map);
    v1.random();
    v2.random();
    // set values in both to 1.0
    // for the first, do via putScalar()
    // the the second, do manually, looping over all elements
    // verify that both have identical values
    v1.putScalar(SCT::one());
    {
      ArrayRCP<Scalar> view = v2.get1dViewNonConst();
      for (typename ArrayRCP<Scalar>::iterator v = view.begin(); v != view.end(); ++v) {
        *v = SCT::one();
      }
      view = Teuchos::null;
    }
    Magnitude err;
    // subtract v2 from v1; this should result in v1 = zeros
    v1.update(-1.0,v2,1.0);
    err = v1.norm2();
    TEST_EQUALITY_CONST(err,SCT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, SingleVecNormalize, Ordinal, Scalar )
  {
    // this documents a usage case in Anasazi::SVQBOrthoManager, which was failing
    // error turned out to be a neglected return in both implementations of update(), 
    // after passing the buck to scale() in the case of alpha==0 or beta==0 or gamma=0
    if (ScalarTraits<Scalar>::isOrdinal) return;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Magnitude M1  = ScalarTraits<Magnitude>::one();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 10;
    const Teuchos_Ordinal numVectors = 6;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    // create random MV
    MV mv(node,map,numVectors);
    mv.random();
    // compute the norms
    Array<Magnitude> norms(numVectors);
    mv.norm2(norms());
    for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
      // get a view of column j, normalize it using update()
      RCP<MV> mvj = mv.subViewNonConst(tuple<Teuchos_Ordinal>(j));
      switch (j) {
      case 0:
        mvj->scale( M1/norms[j] );
        break;
      case 1:
        mvj->update( M1/norms[j], *mvj, M0 );
        break;
      case 2:
        mvj->update( M0         , *mvj, M1/norms[j] );
        break;
      case 3:
        mvj->update( M0         , *mvj, M1/norms[j], *mvj, M0 );
        break;
      case 4:
        mvj->update( M1/norms[j], *mvj, M0         , *mvj, M0 );
        break;
      case 5:
        mvj->update( M0         , *mvj, M0         , *mvj, M1/norms[j] );
        break;
      }
    }
    mv.norm2(norms()); // should be all one now
    Array<Magnitude> ones(numVectors,M1);
    TEST_COMPARE_FLOATING_ARRAYS(norms,ones,ScalarTraits<Magnitude>::eps()*as<Magnitude>(10.));
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountDot, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Teuchos_Ordinal numVectors = 3;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    Array<Scalar> values(6);
    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    //                               [0 1 2]
    // dot(values,values) = [0*0+0*0 1+1*1 + 2*2+2*2] = [0 2 8]
    // summed over all procs, this is [0 2*nprocs 8*nprocs]
    values[0] = as<Scalar>(0);
    values[1] = as<Scalar>(0);
    values[2] = as<Scalar>(1);
    values[3] = as<Scalar>(1);
    values[4] = as<Scalar>(2);
    values[5] = as<Scalar>(2);
    MV mvec1(node,map,values(),2,numVectors),
       mvec2(node,map,values(),2,numVectors);
    Array<Scalar> dots1(numVectors), dots2(numVectors), answer(numVectors);
    answer[0] = as<Scalar>(0);
    answer[1] = as<Scalar>(2*numImages);
    answer[2] = as<Scalar>(8*numImages);
    // do the dots
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,answer,M0);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountDotNonTrivLDA, Ordinal, Scalar )
  {
    // same as CountDot, but the A,LDA has a non-trivial LDA (i.e., LDA != myLen)
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Teuchos_Ordinal numVectors = 3;
    const Teuchos_Ordinal LDA = 3;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    Array<Scalar> values(9);
    // A = {0, 0, -1, 1, 1, -1, 2, 2, -1} = [0   1  2]
    //                                      [0   1  2]
    //                                      [-1 -1 -1]
    // processed as a 2 x 3 with LDA==3, the result it
    //            values =       [0 1 2]
    //                           [0 1 2]
    // dot(values,values) = [0*0+0*0 1+1*1 + 2*2+2*2] = [0 2 8]
    // summed over all procs, this is [0 2*nprocs 8*nprocs]
    values[0] = as<Scalar>(0);
    values[1] = as<Scalar>(0);
    values[2] = as<Scalar>(-1);
    values[3] = as<Scalar>(1);
    values[4] = as<Scalar>(1);
    values[5] = as<Scalar>(-1);
    values[6] = as<Scalar>(2);
    values[7] = as<Scalar>(2);
    values[8] = as<Scalar>(-1);
    MV mvec1(node,map,values(),LDA,numVectors),
       mvec2(node,map,values(),LDA,numVectors);
    Array<Scalar> dots1(numVectors), dots2(numVectors), answer(numVectors);
    answer[0] = as<Scalar>(0);
    answer[1] = as<Scalar>(2*numImages);
    answer[2] = as<Scalar>(8*numImages);
    // do the dots
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,answer,M0);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountNorm1, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Ordinal numVectors = 3;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    Array<Scalar> values(6);
    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    //                               [0 1 2]
    // norm1(values) = [0 2 4]
    // over all procs, this is [0 2*nprocs 4*nprocs]
    values[0] = as<Scalar>(0);
    values[1] = as<Scalar>(0);
    values[2] = as<Scalar>(1);
    values[3] = as<Scalar>(1);
    values[4] = as<Scalar>(2);
    values[5] = as<Scalar>(2);
    MV mvec(node,map,values(),2,numVectors);
    Array<MT> norms(numVectors), answer(numVectors);
    answer[0] = as<MT>(0);
    answer[1] = as<MT>(2*numImages);
    answer[2] = as<MT>(4*numImages);
    // do the dots
    mvec.norm1(norms());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(norms,answer,M0);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountNormInf, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Ordinal numVectors = 3;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    Array<Scalar> values(6);
    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    //                               [0 1 2]
    // normInf(values) = [0 1 2]
    // over all procs, this is [0 1 2]
    values[0] = as<Scalar>(0);
    values[1] = as<Scalar>(0);
    values[2] = as<Scalar>(1);
    values[3] = as<Scalar>(1);
    values[4] = as<Scalar>(2);
    values[5] = as<Scalar>(2);
    MV mvec(node,map,values(),2,numVectors);
    Array<MT> norms(numVectors), answer(numVectors);
    answer[0] = as<MT>(0);
    answer[1] = as<MT>(1);
    answer[2] = as<MT>(2);
    // do the dots
    mvec.normInf(norms());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(norms,answer,M0);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, Norm2, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Ordinal numVectors = 2;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    MV mvec(node,map,2,numVectors);
    // randomize the multivector
    mvec.random();
    // take norms; they should not be zero
    Array<MT> normsRand(numVectors), normsZero(numVectors);
    mvec.norm2(normsRand());
    // zero the vector
    mvec.putScalar(ScalarTraits<Scalar>::zero());
    // take norms; they should be zero
    mvec.norm2(normsZero());
    // check the answers
    bool local_success = true;
    for (Teuchos_Ordinal i=0; i<numVectors; ++i) {
      TEST_ARRAY_ELE_INEQUALITY(normsRand,i,M0);
      TEST_ARRAY_ELE_EQUALITY(normsZero,i,M0);
    }
    success &= local_success;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, MinMaxMean, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Ordinal numLocal = 2;
    const Teuchos_Ordinal numVectors = 2;
    const Teuchos_Ordinal LDA = 2;
    Map<Ordinal> map(INVALID,numLocal,indexBase,comm);
    Array<Scalar> values(4);
    // on proc i of n:
    // values = {i, i+1, n-i, n-i-1} = [ i   n-i ]
    //                                 [i+1 n-i-1]
    // min values are [0,0], but from different procs
    // max values are [numImages,numImages], again from different procs
    // mean values are 
    //  1  n-1
    // --  sum i+i+1 == [(n-1)*n + n]/2n == (n-1)/2 
    // 2n  i=0
    //               == n*n/2n == n/2
    values[0] = as<Scalar>(myImageID);
    values[1] = as<Scalar>(myImageID+1);
    values[2] = as<Scalar>(numImages-myImageID);
    values[3] = as<Scalar>(numImages-myImageID-1);
    MV mvec(node,map,values(),LDA,numVectors);
    // FINISH this test
  }


#define PRINT_TYPE_AND_VALUE(val) { out << std::setw(30) << #val << std::setw(30) << Teuchos::typeName(val) << ": " << val << endl; }



  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadCombinations, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    Node &node = getDefaultNode();
    // create a Map
    const Ordinal indexBase = 0;
    const Scalar rnd = ScalarTraits<Scalar>::random();
    // two maps: one has two entires per node, the other disagrees on node 0
    Map<Ordinal> map1(INVALID,2,indexBase,comm),
                 map2(INVALID,(myImageID == 0 ? 1 : 2),indexBase,comm);
    // multivectors from different maps are incompatible for all ops
    // multivectors from the same map are compatible only if they have the same number of
    //    columns
    MV m1n1(node,map1,1), m1n2(node,map1,2), m2n2(node,map2,2), m1n2_2(node,map1,2);
    Array<Scalar> dots(1);
    Array<Mag>    norms(1);
    // FINISH: test multiply (both), reciprocalMultiply
    TEST_THROW(m1n2.dot(m1n1,dots()), std::runtime_error); // dot
    TEST_THROW(m1n2.dot(m2n2,dots()), std::runtime_error);
    TEST_THROW(m1n2.abs(m1n1), std::runtime_error);       // abs
    TEST_THROW(m1n2.abs(m2n2), std::runtime_error);
    TEST_THROW(m1n2.abs(m1n1), std::runtime_error);       // abs
    TEST_THROW(m1n2.abs(m2n2), std::runtime_error);
    TEST_THROW(m1n2.scale(rnd,m1n1), std::runtime_error); // abs
    TEST_THROW(m1n2.scale(rnd,m2n2), std::runtime_error);
    TEST_THROW(m1n2.update(rnd,m1n1,rnd), std::runtime_error); // update(alpha,A,beta)
    TEST_THROW(m1n2.update(rnd,m2n2,rnd), std::runtime_error);
    TEST_THROW(m1n2.update(rnd,m2n2  ,rnd,m1n2_2,rnd), std::runtime_error); // update(alpha,A,beta,B,gamma) // A incompat
    TEST_THROW(m1n2.update(rnd,m2n2  ,rnd,m1n2_2,rnd), std::runtime_error); // incompt is length            // A incompat
    TEST_THROW(m1n2.update(rnd,m1n2_2,rnd,m2n2  ,rnd), std::runtime_error);                                 // B incompat
    TEST_THROW(m1n2.update(rnd,m1n2_2,rnd,m2n2  ,rnd), std::runtime_error);                                 // B incompat
    TEST_THROW(m1n2.update(rnd,m2n2  ,rnd,m2n2  ,rnd), std::runtime_error);                                 // A,B incompat
    TEST_THROW(m1n2.update(rnd,m2n2  ,rnd,m2n2  ,rnd), std::runtime_error);                                 // A,B incompat
    TEST_THROW(m1n2.update(rnd,m1n1  ,rnd,m1n2_2,rnd), std::runtime_error); // incompt is numVecs           // A incompat
    TEST_THROW(m1n2.update(rnd,m1n1  ,rnd,m1n2_2,rnd), std::runtime_error);                                 // A incompat
    TEST_THROW(m1n2.update(rnd,m1n2_2,rnd,m1n1  ,rnd), std::runtime_error);                                 // B incompat
    TEST_THROW(m1n2.update(rnd,m1n2_2,rnd,m1n1  ,rnd), std::runtime_error);                                 // B incompat
    TEST_THROW(m1n2.update(rnd,m1n1  ,rnd,m1n1  ,rnd), std::runtime_error);                                 // A,B incompat
    TEST_THROW(m1n2.update(rnd,m1n1  ,rnd,m1n1  ,rnd), std::runtime_error);                                 // A,B incompat
    TEST_THROW(m1n1.normWeighted(m1n2,norms()), std::runtime_error);        // normWeighted
    TEST_THROW(m1n2.normWeighted(m2n2,norms()), std::runtime_error);
    TEST_THROW(m1n2.reciprocal(m1n1), std::runtime_error);                  // reciprocal
    TEST_THROW(m1n2.reciprocal(m2n2), std::runtime_error);
  }

// 
// INSTANTIATIONS
//

#ifdef HAVE_TEUCHOS_BLASFLOAT
#  define UNIT_TEST_GROUP_ORDINAL_FLOAT(ORDINAL)\
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
#else
#  define UNIT_TEST_GROUP_ORDINAL_FLOAT(ORDINAL)
#endif

#ifdef HAVE_TEUCHOS_COMPLEX
#  ifdef HAVE_TEUCHOS_BLASFLOAT
#    define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)\
       typedef std::complex<float> ComplexFloat; \
       UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexFloat)
#  else
#    define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)
#  endif
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
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, basic             , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadConstNumVecs   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadConstLDA       , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadConstAA        , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CopyConst         , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(      Vector, CopyConst         , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(      Vector, Indexing          , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, OrthoDot          , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountDot          , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountDotNonTrivLDA, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadDot            , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountNorm1        , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountNormInf      , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, Norm2             , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, ZeroScaleUpdate   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(      Vector, ZeroScaleUpdate   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadCombinations   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadMultiply       , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, SingleVecNormalize, ORDINAL, SCALAR ) \
      /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, Multiply          , ORDINAL, SCALAR )*/ \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, LabeledObject     , ORDINAL, ORDINAL, SCALAR ) 


#ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)
     UNIT_TEST_GROUP_ORDINAL(int)

#else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_FLOAT(ORDINAL) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)  \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
     UNIT_TEST_GROUP_ORDINAL(int)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)
     typedef long int LongInt;   UNIT_TEST_GROUP_ORDINAL(LongInt)
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt; UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#    endif

#endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
