#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

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
  using Teuchos::SerialDenseMatrix;
  using Teuchos::arrayView;
  using std::copy;
  using std::ostream_iterator;

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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, basic, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    const Ordinal numVecs  = 5;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
    MV mvec(map,numVecs,true);
    TEST_EQUALITY( mvec.numVectors(), numVecs );
    TEST_EQUALITY( mvec.myLength(), numLocal );
    TEST_EQUALITY( mvec.globalLength(), numImages*numLocal );
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
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
    TEST_THROW(MV mvec(map,ZERO), std::invalid_argument);
    TEST_THROW(MV mvec(map,INVALID), std::invalid_argument);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstLDA, Ordinal, Scalar )
  {
    // numlocal > LDA
    // ergo, the arrayview doesn't contain enough data to specify the entries
    // also, if bounds checking is enabled, check that bad bounds are caught
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef Tpetra::Vector<Ordinal,Scalar> V;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Ordinal TWO = ONE + ONE;
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVecs = TWO;
    // multivector has two vectors, each proc having two values per vector
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
    // we need 4 scalars to specify values on each proc
    Array<Scalar> values(4);
#ifdef TEUCHOS_DEBUG
    // too small an ArrayView (less than 4 values) is met with an exception, if debugging is on
    TEST_THROW(MV mvec(map,values(0,3),TWO,numVecs), std::runtime_error);
    // it could also be too small for the given LDA: 
    TEST_THROW(MV mvec(map,values(),TWO+ONE,numVecs), std::runtime_error);
    // too small for number of entries in a Vector
    TEST_THROW(V   vec(map,values(0,1)), std::runtime_error);
#endif
    // LDA < numLocal throws an exception anytime
    TEST_THROW(MV mvec(map,values(0,4),ONE,numVecs), std::invalid_argument);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadMultiply, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    const Scalar SONE = ScalarTraits<Scalar>::one(),
                SZERO = ScalarTraits<Scalar>::zero();
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    {
      // create local Maps
      Map<Ordinal> map3l(3*ONE,ZERO,platform,true),
                   map2l(2*ONE,ZERO,platform,true);
      MV mvecA(map3l,2*ONE),
         mvecB(map2l,3*ONE),
         mvecD(map2l,2*ONE);
      // failures, 8 combinations:
      // [NTC],[NTC]: A,B don't match
      // [NTC],[NTC]: C doesn't match A,B
      TEST_THROW( mvecD.multiply(NO_TRANS  ,NO_TRANS  ,SONE,mvecA,mvecA,SZERO), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(NO_TRANS  ,CONJ_TRANS,SONE,mvecA,mvecB,SZERO), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(CONJ_TRANS,NO_TRANS  ,SONE,mvecB,mvecA,SZERO), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(CONJ_TRANS,CONJ_TRANS,SONE,mvecB,mvecB,SZERO), std::runtime_error);   // 2x2: 3x2 x 3x2
      TEST_THROW( mvecD.multiply(NO_TRANS  ,NO_TRANS  ,SONE,mvecA,mvecB,SZERO), std::runtime_error);   // 2x2: 3x2 x 2x3
      TEST_THROW( mvecD.multiply(NO_TRANS  ,CONJ_TRANS,SONE,mvecA,mvecA,SZERO), std::runtime_error);   // 2x2: 3x2 x 2x3
      TEST_THROW( mvecD.multiply(CONJ_TRANS,NO_TRANS  ,SONE,mvecB,mvecB,SZERO), std::runtime_error);   // 2x2: 3x2 x 2x3
      TEST_THROW( mvecD.multiply(CONJ_TRANS,CONJ_TRANS,SONE,mvecB,mvecA,SZERO), std::runtime_error);   // 2x2: 3x2 x 2x3
    }
    // case 2: C(local) = A^T(distr) * B  (distr)  : one of these
    {
      Map<Ordinal> map3n(INVALID,3*ONE,ZERO,platform),
                   map2n(INVALID,2*ONE,ZERO,platform),
                   map2l(2*ONE,ZERO,platform,true),
                   map3l(3*ONE,ZERO,platform,true);
      MV mv3nx2(map3n,2*ONE),
         mv2nx2(map2n,2*ONE),
         mv2lx2(map2l,2*ONE),
         mv2lx3(map2l,3*ONE),
         mv3lx2(map3l,2*ONE),
         mv3lx3(map3l,3*ONE);
      // non-matching input lengths
      TEST_THROW( mv2lx2.multiply(CONJ_TRANS,NO_TRANS,SONE,mv3nx2,mv2nx2,SZERO), std::runtime_error);   // (2 x 3n) x (2n x 2) not compat
      TEST_THROW( mv2lx2.multiply(CONJ_TRANS,NO_TRANS,SONE,mv2nx2,mv3nx2,SZERO), std::runtime_error);   // (2 x 2n) x (3n x 2) not compat
      // non-matching output size
      TEST_THROW( mv3lx3.multiply(CONJ_TRANS,NO_TRANS,SONE,mv3nx2,mv3nx2,SZERO), std::runtime_error);   // (2 x 3n) x (3n x 2) doesn't fit 3x3
      TEST_THROW( mv3lx2.multiply(CONJ_TRANS,NO_TRANS,SONE,mv3nx2,mv3nx2,SZERO), std::runtime_error);   // (2 x 3n) x (3n x 2) doesn't fit 3x2
      TEST_THROW( mv2lx3.multiply(CONJ_TRANS,NO_TRANS,SONE,mv3nx2,mv3nx2,SZERO), std::runtime_error);   // (2 x 3n) x (3n x 2) doesn't fit 2x3
    }
    // case 3: C(distr) = A  (distr) * B^X(local)  : two of these
    {
      Map<Ordinal> map3n(INVALID,3*ONE,ZERO,platform),
                   map2n(INVALID,2*ONE,ZERO,platform),
                   map2l(2*ONE,ZERO,platform,true),
                   map3l(3*ONE,ZERO,platform,true);
      MV mv3nx2(map3n,2*ONE),
         mv2nx2(map2n,2*ONE),
         mv2x3(map2l,3*ONE),
         mv3x2(map3l,2*ONE);
      // non-matching input lengths
      TEST_THROW( mv3nx2.multiply(NO_TRANS,CONJ_TRANS,SONE,mv3nx2,mv2x3,SZERO), std::runtime_error);   // (3n x 2) x (3 x 2) (trans) not compat
      TEST_THROW( mv3nx2.multiply(NO_TRANS,NO_TRANS  ,SONE,mv3nx2,mv3x2,SZERO), std::runtime_error);   // (3n x 2) x (3 x 2) (nontrans) not compat
      // non-matching output sizes
      TEST_THROW( mv3nx2.multiply(NO_TRANS,CONJ_TRANS,SONE,mv3nx2,mv3x2,SZERO), std::runtime_error);   // (3n x 2) x (2 x 3) doesn't fit 3nx2
      TEST_THROW( mv3nx2.multiply(NO_TRANS,NO_TRANS  ,SONE,mv3nx2,mv2x3,SZERO), std::runtime_error);   // (3n x 2) x (2 x 3) doesn't fit 3nx2
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, Multiply, Ordinal, Scalar )
  {
    using Teuchos::View;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    // create a Map
    Map<Ordinal> map3n(INVALID,3*ONE,ZERO,platform),
                 map2n(INVALID,2*ONE,ZERO,platform),
                lmap3(3*ONE,ZERO,platform,true),
                lmap2(2*ONE,ZERO,platform,true);
    const Scalar SONE = ScalarTraits<Scalar>::one(),
                SZERO = ScalarTraits<Scalar>::zero();
    const Mag   MZERO = ScalarTraits<Mag>::zero();
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    // deterministic input/output
    {
      MV mv3x2l(lmap3,2*ONE),
         mv2x3l(lmap2,3*ONE),
         mv2x2l(lmap2,2*ONE),
         mv3x3l(lmap3,3*ONE);
      // fill multivectors with ones
      mv3x2l.putScalar(ScalarTraits<Scalar>::one());
      mv2x3l.putScalar(ScalarTraits<Scalar>::one());
      // fill expected answers Array
      ArrayView<const Scalar> tmpView(Teuchos::null); Ordinal dummy;
      Teuchos::Array<Scalar> check2(4,3*ONE); // each entry (of four) is the product [1 1 1]*[1 1 1]' = 3
      Teuchos::Array<Scalar> check3(9,2*ONE); // each entry (of nine) is the product [1 1]*[1 1]' = 2
      // test
      mv3x3l.multiply(NO_TRANS  ,NO_TRANS  ,SONE,mv3x2l,mv2x3l,SZERO);
      mv3x3l.extractConstView(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check3,MZERO);
      mv2x2l.multiply(NO_TRANS  ,CONJ_TRANS,SONE,mv2x3l,mv2x3l,SZERO);
      mv2x2l.extractConstView(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check2,MZERO);
      mv2x2l.multiply(CONJ_TRANS,NO_TRANS  ,SONE,mv3x2l,mv3x2l,SZERO);
      mv2x2l.extractConstView(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check2,MZERO);
      mv3x3l.multiply(CONJ_TRANS,CONJ_TRANS,SONE,mv2x3l,mv3x2l,SZERO);
      mv3x3l.extractConstView(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check3,MZERO);
    }
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    // random input/output
    {
      ArrayView<Scalar> tmvView(Teuchos::null), sdmView(Teuchos::null);
      Ordinal stride;
      MV tmv3x2(lmap3,2*ONE),
         tmv2x3(lmap2,3*ONE),
         tmv2x2(lmap2,2*ONE),
         tmv3x3(lmap3,3*ONE);
      // setup serial dense matrices 
      tmv3x2.extractView(tmvView,stride); SerialDenseMatrix<int,Scalar> sdm3x2(View,tmvView.getRawPtr(),stride,3*ONE,2*ONE);
      tmv2x3.extractView(tmvView,stride); SerialDenseMatrix<int,Scalar> sdm2x3(View,tmvView.getRawPtr(),stride,2*ONE,3*ONE);
      SerialDenseMatrix<int,Scalar> sdm2x2(2*ONE,2*ONE), sdm3x3(3*ONE,3*ONE);
      // fill multivectors with ones (fills sdm3x2 and sdm2x3 as well, because of view semantics)
      tmv3x2.random();
      tmv2x3.random();
      // test: perform MV multiply and SDM multiply, then check that answers are equivalent
      {
        tmv3x3.multiply(NO_TRANS,NO_TRANS,SONE,tmv3x2,tmv2x3,SZERO);
        sdm3x3.multiply(NO_TRANS,NO_TRANS,SONE,sdm3x2,sdm2x3,SZERO);
        tmv3x3.extractView(tmvView,stride); sdmView = arrayView(sdm3x3.values(),sdm3x3.numRows()*sdm3x3.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmvView,sdmView,MZERO);
      }
      {
        tmv2x2.multiply(NO_TRANS,CONJ_TRANS,SONE,tmv2x3,tmv2x3,SZERO);
        sdm2x2.multiply(NO_TRANS,CONJ_TRANS,SONE,sdm2x3,sdm2x3,SZERO);
        tmv2x2.extractView(tmvView,stride); sdmView = arrayView(sdm2x2.values(),sdm2x2.numRows()*sdm2x2.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmvView,sdmView,MZERO);
      }
      {
        tmv2x2.multiply(CONJ_TRANS,NO_TRANS,SONE,tmv3x2,tmv3x2,SZERO);
        sdm2x2.multiply(CONJ_TRANS,NO_TRANS,SONE,sdm3x2,sdm3x2,SZERO);
        tmv2x2.extractView(tmvView,stride); sdmView = arrayView(sdm2x2.values(),sdm2x2.numRows()*sdm2x2.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmvView,sdmView,MZERO);
      }
      {
        tmv3x3.multiply(CONJ_TRANS,CONJ_TRANS,SONE,tmv2x3,tmv3x2,SZERO);
        sdm3x3.multiply(CONJ_TRANS,CONJ_TRANS,SONE,sdm2x3,sdm3x2,SZERO);
        tmv3x3.extractView(tmvView,stride); sdmView = arrayView(sdm3x3.values(),sdm3x3.numRows()*sdm3x3.numCols());
        TEST_COMPARE_FLOATING_ARRAYS(tmvView,sdmView,MZERO);
      }
    }
    // case 2: C(local) = A^T(distr) * B  (distr)  : one of these
    {
      MV mv3nx2(map3n,2*ONE),
         mv3nx3(map3n,3*ONE),
         // locals
         mv2x2(lmap2,2*ONE),
         mv2x3(lmap2,3*ONE),
         mv3x2(lmap3,2*ONE),
         mv3x3(lmap3,3*ONE);
      // fill multivectors with ones
      mv3nx3.putScalar(ScalarTraits<Scalar>::one());
      mv3nx2.putScalar(ScalarTraits<Scalar>::one());
      // fill expected answers Array
      ArrayView<const Scalar> tmpView(Teuchos::null); Ordinal dummy;
      Teuchos::Array<Scalar> check(9,3*ONE*numImages);
      // test
      mv2x2.multiply(CONJ_TRANS,NO_TRANS,SONE,mv3nx2,mv3nx2,SZERO); 
      mv2x2.extractConstView(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),MZERO);
      mv2x3.multiply(CONJ_TRANS,NO_TRANS,SONE,mv3nx2,mv3nx3,SZERO);
      mv2x3.extractConstView(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),MZERO);
      mv3x2.multiply(CONJ_TRANS,NO_TRANS,SONE,mv3nx3,mv3nx2,SZERO);
      mv3x2.extractConstView(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),MZERO);
      mv3x3.multiply(CONJ_TRANS,NO_TRANS,SONE,mv3nx3,mv3nx3,SZERO);
      mv3x3.extractConstView(tmpView,dummy); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),MZERO);
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
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Ordinal TWO = ONE + ONE;
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    // multivector has two vectors, each proc having two values per vector
    Map<Ordinal> map2(INVALID,numLocal  ,indexBase,platform),
                 map3(INVALID,numLocal+1,indexBase,platform);
    // we need 4 scalars to specify values on each proc
    Array<Scalar> values(4);
    Array<ArrayView<const Scalar> > arrOfarr(2,ArrayView<const Scalar>(Teuchos::null));
    Array<ArrayView<const Scalar> > emptyArr;
    arrOfarr[0] = values(0,2);
    arrOfarr[1] = values(2,2);
#ifdef TEUCHOS_DEBUG
    // arrOfarr.size() == 0
    TEST_THROW(MV mvec(map2,emptyArr()), std::runtime_error);
    // individual ArrayViews could be too small
    TEST_THROW(MV mvec(map3,arrOfarr()), std::runtime_error);
#else
    (void)success;
    (void)out;
#endif
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadDot, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef Tpetra::Vector<Ordinal,Scalar> V;
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = as<Ordinal>(0);
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    Map<Ordinal> map1(INVALID,as<Ordinal>(1),indexBase,platform),
                 map2(INVALID,as<Ordinal>(2),indexBase,platform);
    {
      MV mv12(map1,as<Ordinal>(1)),
         mv21(map2,as<Ordinal>(1)),
         mv22(map2,as<Ordinal>(2));
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
      V v1(map1),
        v2(map2);
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
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Ordinal TWO = ONE + ONE;
    const Scalar SZERO = Teuchos::ScalarTraits<Scalar>::zero();
    const Mag MZERO = Teuchos::ScalarTraits<Mag>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVectors = TWO;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
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
    MV mvec1(map,values(0,4),TWO,numVectors),
       mvec2(map,values(1,4),TWO,numVectors);
    Array<Scalar> dots1(numVectors), dots2(numVectors), zeros(numVectors);
    std::fill(zeros.begin(),zeros.end(),Teuchos::ScalarTraits<Scalar>::zero());
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,MZERO);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,MZERO);
    TEST_EQUALITY_CONST( mvec1(0)->dot(*mvec2(0)), SZERO);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, ZeroScaleUpdate, Ordinal, Scalar )
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Ordinal TWO = ONE + ONE;
    const Mag MZERO = Teuchos::ScalarTraits<Mag>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVectors = TWO;
    const Ordinal LDA = TWO;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
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
    MV A(map,values(0,4),LDA,numVectors),
       B(map,values(2,4),LDA,numVectors);
    Array<Mag> norms(numVectors), zeros(numVectors);
    std::fill(zeros.begin(),zeros.end(),MZERO);
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
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MZERO);
    }
    //   set A2 = A
    //   check that it equals B: scale,subtraction in situ
    {
      MV A2(A);
      A2.update(as<Scalar>(-1),B,as<Scalar>(2));
      A2.norm2(norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MZERO);
    }
    //   set C random
    //   set it to zero by combination with A,B
    {
      MV C(map,numVectors);
      C.random();
      C.update(as<Scalar>(-1),B,as<Scalar>(2),A,as<Scalar>(0));
      C.norm2(norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MZERO);
    }
    //   set C random
    //   scale it ex-situ
    //   check that it equals B: subtraction in situ
    {
      MV C(map,numVectors);
      C.scale(as<Scalar>(2),A);
      C.update(as<Scalar>(1),B,as<Scalar>(-1));
      C.norm2(norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MZERO);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, ZeroScaleUpdate, Ordinal, Scalar )
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::Vector<Ordinal,Scalar> V;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Mag MZERO = Teuchos::ScalarTraits<Mag>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    Map<Ordinal> map(INVALID,2,0,platform);
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
    V A(map,values(0,2)),
      B(map,values(2,2));
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
      TEST_EQUALITY(norm,MZERO);
      TEST_EQUALITY(norm,norms[0]);
    }
    //   set A2 = A
    //   check that it equals B: scale,subtraction in situ
    {
      V A2(A);
      A2.update(as<Scalar>(-1),B,as<Scalar>(2));
      norm = A2.norm2(); A2.norm2(norms());
      TEST_EQUALITY(norm,MZERO);
      TEST_EQUALITY(norm,norms[0]);
    }
    //   set C random
    //   set it to zero by combination with A,B
    {
      V C(map);
      C.random();
      C.update(as<Scalar>(-1),B,as<Scalar>(2),A,as<Scalar>(0));
      norm = C.norm2(); C.norm2(norms());
      TEST_EQUALITY(norm,MZERO);
      TEST_EQUALITY(norm,norms[0]);
    }
    //   set C random
    //   scale it ex-situ
    //   check that it equals B: subtraction in situ
    {
      V C(map);
      C.random();
      C.scale(as<Scalar>(2),A);
      C.update(as<Scalar>(1),B,as<Scalar>(-1));
      norm = C.norm2(); C.norm2(norms());
      TEST_EQUALITY(norm,MZERO);
      TEST_EQUALITY(norm,norms[0]);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CopyConst, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Ordinal TWO = ONE + ONE;
    const Magnitude MZERO = ScalarTraits<Magnitude>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVectors = TWO;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
    // create random MV
    MV morig(map,numVectors);
    morig.random();
    // copy it
    MV mcopy1(morig), mcopy2(morig);
    // verify that all three have identical values
    Array<Magnitude> norig(numVectors), ncopy1(numVectors), ncopy2(numVectors);
    morig.normInf(norig);
    mcopy1.normInf(ncopy1);
    mcopy2.normInf(ncopy2);
    TEST_COMPARE_FLOATING_ARRAYS(norig,ncopy1,MZERO);
    TEST_COMPARE_FLOATING_ARRAYS(norig,ncopy2,MZERO);
    TEST_COMPARE_FLOATING_ARRAYS(ncopy1,ncopy2,MZERO);
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
    for (Ordinal i=ZERO; i<numVectors; ++i) {
      TEST_ARRAY_ELE_EQUALITY( norig,  i, as<Scalar>(0) );
      TEST_ARRAY_ELE_EQUALITY( ncopy1, i, as<Scalar>(1) );
      TEST_ARRAY_ELE_EQUALITY( ncopy2, i, as<Scalar>(2) );
    }
    success &= local_success;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, CopyConst, Ordinal, Scalar )
  {
    typedef Tpetra::Vector<Ordinal,Scalar> V;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    Map<Ordinal> map(INVALID,2,0,platform);
    // create random MV
    V morig(map);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, SingleVecNormalize, Ordinal, Scalar )
  {
    // this documents a usage case in Anasazi::SVQBOrthoManager, which was failing
    // error turned out to be a neglected return in both implementations of update(), 
    // after passing to buck to scale() in the case of alpha==0 or beta==0 or gamma=0
    if (as<Scalar>(1.0) / as<Scalar>(3.0) == ScalarTraits<Scalar>::zero()) return;
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Magnitude MONE  = ScalarTraits<Magnitude>::one();
    const Magnitude MZERO = ScalarTraits<Magnitude>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = as<Ordinal>(10);
    const Ordinal numVectors = as<Ordinal>(6);
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
    // create random MV
    MV mv(map,numVectors);
    mv.random();
    // compute the norms
    Array<Magnitude> norms(numVectors);
    mv.norm2(norms());
    for (Ordinal j=ZERO; j<numVectors; ++j) {
      // get a view of column j, normalize it using update()
      Array<Teuchos_Ordinal> ind(1,j);
      RCP<MV> mvj = mv.subView(ind());
      switch (j){
      case 0:
        mvj->scale( MONE/norms[j] );
        break;
      case 1:
        mvj->update( MONE/norms[j], *mvj, MZERO );
        break;
      case 2:
        mvj->update( MZERO        , *mvj, MONE/norms[j] );
        break;
      case 3:
        mvj->update( MZERO        , *mvj, MONE/norms[j], *mvj, MZERO );
        break;
      case 4:
        mvj->update( MONE/norms[j], *mvj, MZERO        , *mvj, MZERO );
        break;
      case 5:
        mvj->update( MZERO        , *mvj, MZERO        , *mvj, MONE/norms[j] );
        break;
      }
    }
    mv.norm2(norms()); // should be all one now
    Array<Magnitude> ones(numVectors,MONE);
    TEST_COMPARE_FLOATING_ARRAYS(norms,ones,ScalarTraits<Magnitude>::eps()*as<Magnitude>(10.));
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountDot, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal TWO = ONE+ONE;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Magnitude MZERO = ScalarTraits<Magnitude>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = ONE+ONE;
    const Ordinal numVectors = ONE+ONE+ONE;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
    Array<Scalar> values(6);
    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    //                               [0 1 2]
    // dot(values,values) = [0*0+0*0 1*1+1*1 + 2*2+2*2] = [0 2 8]
    // summed over all procs, this is [0 2*nprocs 8*nprocs]
    values[0] = as<Scalar>(0);
    values[1] = as<Scalar>(0);
    values[2] = as<Scalar>(1);
    values[3] = as<Scalar>(1);
    values[4] = as<Scalar>(2);
    values[5] = as<Scalar>(2);
    MV mvec1(map,values(),TWO,numVectors),
       mvec2(map,values(),TWO,numVectors);
    Array<Scalar> dots1(numVectors), dots2(numVectors), answer(numVectors);
    answer[0] = as<Scalar>(0);
    answer[1] = as<Scalar>(2*numImages);
    answer[2] = as<Scalar>(8*numImages);
    // do the dots
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,MZERO);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,answer,MZERO);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountDotNonTrivLDA, Ordinal, Scalar )
  {
    // same as CountDot, but the A,LDA has a non-trivial LDA (i.e., LDA != myLen)
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal TWO = ONE+ONE;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const Magnitude MZERO = ScalarTraits<Magnitude>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = ONE+ONE;
    const Ordinal numVectors = ONE+ONE+ONE;
    const Ordinal LDA = ONE+TWO;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
    Array<Scalar> values(9);
    // A = {0, 0, -1, 1, 1, -1, 2, 2, -1} = [0   1  2]
    //                                      [0   1  2]
    //                                      [-1 -1 -1]
    // processed as a 2 x 3 with LDA==3, the result it
    //            values =       [0 1 2]
    //                           [0 1 2]
    // dot(values,values) = [0*0+0*0 1*1+1*1 + 2*2+2*2] = [0 2 8]
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
    MV mvec1(map,values(),LDA,numVectors),
       mvec2(map,values(),LDA,numVectors);
    Array<Scalar> dots1(numVectors), dots2(numVectors), answer(numVectors);
    answer[0] = as<Scalar>(0);
    answer[1] = as<Scalar>(2*numImages);
    answer[2] = as<Scalar>(8*numImages);
    // do the dots
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,MZERO);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,answer,MZERO);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountNorm1, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal TWO = ONE+ONE;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT MZERO = ScalarTraits<MT>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = ONE+ONE;
    const Ordinal numVectors = ONE+ONE+ONE;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
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
    MV mvec(map,values(),TWO,numVectors);
    Array<MT> norms(numVectors), answer(numVectors);
    answer[0] = as<MT>(0);
    answer[1] = as<MT>(2*numImages);
    answer[2] = as<MT>(4*numImages);
    // do the dots
    mvec.norm1(norms());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(norms,answer,MZERO);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountNormInf, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal TWO = ONE+ONE;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT MZERO = ScalarTraits<MT>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = ONE+ONE;
    const Ordinal numVectors = ONE+ONE+ONE;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
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
    MV mvec(map,values(),TWO,numVectors);
    Array<MT> norms(numVectors), answer(numVectors);
    answer[0] = as<MT>(0);
    answer[1] = as<MT>(1);
    answer[2] = as<MT>(2);
    // do the dots
    mvec.normInf(norms());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(norms,answer,MZERO);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, Norm2, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal TWO = ONE+ONE;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT MZERO = ScalarTraits<MT>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVectors = TWO;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
    MV mvec(map,TWO,numVectors);
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
    for (Ordinal i=ZERO; i<numVectors; ++i) {
      TEST_ARRAY_ELE_INEQUALITY(normsRand,i,MZERO);
      TEST_ARRAY_ELE_EQUALITY(normsZero,i,MZERO);
    }
    success &= local_success;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, MinMaxMean, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal TWO = ONE+ONE;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    const MT MZERO = ScalarTraits<MT>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVectors = TWO;
    const Ordinal LDA = TWO;
    Map<Ordinal> map(INVALID,numLocal,indexBase,platform);
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
    values[3] = as<Scalar>(numImages-myImageID-ONE);
    MV mvec(map,values(),LDA,numVectors);
  }


#define PRINT_TYPE_AND_VALUE(val) { out << std::setw(30) << #val << std::setw(30) << Teuchos::typeName(val) << ": " << val << endl; }



  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadCombinations, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal TWO = ONE+ONE;
    const Ordinal INVALID = OrdinalTraits<Ordinal>::invalid();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int myImageID = comm->getRank();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Scalar rnd = ScalarTraits<Scalar>::random();
    // two maps: one has two entires per node, the other disagrees on node 0
    Map<Ordinal> map1(INVALID,TWO,indexBase,platform),
                 map2(INVALID,(myImageID == ZERO ? ONE : TWO),indexBase,platform);
    // multivectors from different maps are incompatible for all ops
    // multivectors from the same map are compatible only if they have the same number of
    //    columns
    MV m1n1(map1,ONE), m1n2(map1,TWO), m2n2(map2,TWO), m1n2_2(map1,TWO);
    Array<Scalar> dots(ONE);
    Array<Mag>    norms(ONE);
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


  /* TODO 
     Many constructors left to test
     MultiVector (const Map< Ordinal > &map, const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &arrayOfArrays, Ordinal numVectors)

     MultiVector<Ordinal,Scalar> subCopy(const Teuchos::Range1D &colRng) const;
     MultiVector<Ordinal,Scalar> subCopy(const Teuchos::ArrayView<Teuchos_Index> &cols) const;
     MultiVector<Ordinal,Scalar> subView(const Teuchos::Range1D &colRng);
     MultiVector<Ordinal,Scalar> subView(const Teuchos::ArrayView<Teuchos_Index> &cols);
     const MultiVector<Ordinal,Scalar> subViewConst(const Teuchos::Range1D &colRng) const;
     const MultiVector<Ordinal,Scalar> subViewConst(const Teuchos::ArrayView<Teuchos_Index> &cols) const;

     Mod routines left to test
     void replaceGlobalValue (Ordinal globalRow, Ordinal vectorIndex, const Scalar &value)
     void sumIntoGlobalValue (Ordinal globalRow, Ordinal vectorIndex, const Scalar &value)
     void replaceMyValue (Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue)
     void sumIntoMyValue (Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue)

     Arithmetic methods left to test:
     void multiply (Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const MultiVector< Ordinal, Scalar > &B, const Scalar &beta)
     void multiply (const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const MultiVector< Ordinal, Scalar > &B, const Scalar &beta)
     void reciprocalMultiply (const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const MultiVector< Ordinal, Scalar > &B, const Scalar &beta)
  */

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
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, Multiply          , ORDINAL, SCALAR ) 


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
