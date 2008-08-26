#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_MultiVector.hpp"

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
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Tpetra::MultiVector;
  using std::endl;
  using Teuchos::Array;

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
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal NEGONE = ZERO - OrdinalTraits<Ordinal>::one();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    const Ordinal numVecs  = 5;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    MultiVector<Ordinal,Scalar> mvec(map,numVecs);
    TEST_EQUALITY( mvec.numVectors(), numVecs );
    out << mvec << endl;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstNumVecs, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal NEGONE = ZERO - OrdinalTraits<Ordinal>::one();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = 10;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    TEST_THROW(MV mvec(map,ZERO), std::invalid_argument);
    TEST_THROW(MV mvec(map,NEGONE), std::invalid_argument);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstLDA, Ordinal, Scalar )
  {
    // numlocal > LDA
    // ergo, the arrayview doesn't contain enough data to specify the entries
    // also, if bounds checking is enabled, check that bad bounds are caught
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;
    const Ordinal TWO = ONE + ONE;
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVecs = TWO;
    // multivector has two vectors, each proc having two values per vector
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    // we need 4 scalars to specify values on each proc
    Array<Scalar> values(4);
#ifdef TEUCHOS_DEBUG
    // too small an ArrayView (less than 4 values) is met with an exception, if debugging is on
    TEST_THROW(MV mvec(map,values(0,3),TWO,numVecs), std::invalid_argument);
    // it could also be too small for the given LDA: 
    TEST_THROW(MV mvec(map,values(),TWO+ONE,numVecs), std::invalid_argument);
#endif
    // LDA < numLocal throws an exception anytime
    TEST_THROW(MV mvec(map,values(0,4),ONE,numVecs), std::invalid_argument);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, OrthoDot, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;
    const Ordinal TWO = ONE + ONE;
    const Scalar SZERO = ScalarTraits<Scalar>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVectors = TWO;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
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
    std::fill(zeros.begin(),zeros.end(),SZERO);
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,SZERO);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,SZERO);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CopyConst, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;
    const Ordinal TWO = ONE + ONE;
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVectors = TWO;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    // create random MV
    MV morig(map,numVectors);
    morig.random();
    // copy it
    MV mcopy1(morig), mcopy2(morig);
    // verify that all three have identical values
    // modify all three, check independence
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountDot, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal TWO = ONE+ONE;
    const Ordinal NEGONE = ZERO - ONE;
    const Scalar SZERO = ScalarTraits<Scalar>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = ONE+ONE;
    const Ordinal numVectors = ONE+ONE+ONE;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    Array<Scalar> values(6);
    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    //                               [0 1 2]
    // dot(values,values) = [0 2 4]
    // summed over all procs, this is [0 2*nprocs 4*nprocs]
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
    answer[2] = as<Scalar>(4*numImages);
    // do the dots
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,SZERO);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,answer,SZERO);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountNorm1, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal TWO = ONE+ONE;
    const Ordinal NEGONE = ZERO - ONE;
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
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
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
    const Ordinal NEGONE = ZERO - ONE;
    const MT MZERO = ScalarTraits<MT>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = ONE+ONE;
    const Ordinal numVectors = ONE+ONE+ONE;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    Array<Scalar> values(6);
    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    //                               [0 1 2]
    // normInf(values) = [0 1 2]
    // over all procs, this is [0 2 4]
    values[0] = as<Scalar>(0);
    values[1] = as<Scalar>(0);
    values[2] = as<Scalar>(1);
    values[3] = as<Scalar>(1);
    values[4] = as<Scalar>(2);
    values[5] = as<Scalar>(2);
    MV mvec(map,values(),TWO,numVectors);
    Array<MT> norms(numVectors), answer(numVectors);
    answer[0] = as<MT>(0);
    answer[1] = as<MT>(2);
    answer[2] = as<MT>(4);
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
    const Ordinal NEGONE = ZERO - ONE;
    const MT MZERO = ScalarTraits<MT>::zero();
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal = TWO;
    const Ordinal numVectors = TWO;
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform);
    MV mvec(map,TWO,numVectors);
    // randomize the multivector
    mvec.random();
    // take norms; they should not be zero
    Array<MT> normsRand(numVectors), normsZero(numVectors);
    mvec.norm2(normsRand());
    // zero the vector
    mvec.scale(ScalarTraits<Scalar>::zero());
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


  /* TODO 
     Many constructors left to test

     Mod routines left to test
     void replaceGlobalValue (Ordinal globalRow, Ordinal vectorIndex, const Scalar &value)
     void sumIntoGlobalValue (Ordinal globalRow, Ordinal vectorIndex, const Scalar &value)
     void replaceMyValue (Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue)
     void sumIntoMyValue (Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue)
     void putScalar (const Scalar &ScalarConstant)
     void random ()

     Arithmetic methods left to test:
     void dot (const MultiVector< Ordinal, Scalar > &A, Teuchos::Array< Scalar > &dots) const
     void abs (const MultiVector< Ordinal, Scalar > &A)
     void reciprocal (const MultiVector< Ordinal, Scalar > &A)
     void scale (const Scalar &alpha)
     void scale (const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A)
     void update (const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const Scalar &beta)
     void update (const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const Scalar &beta, const MultiVector< Ordinal, Scalar > &B, const Scalar &gamma)
     void norm1 (Teuchos::Array< Scalar > &norms) const
     void norm2 (Teuchos::Array< Scalar > &norms) const
     void normInf (Teuchos::Array< Scalar > &norms) const
     void normWeighted (const MultiVector< Ordinal, Scalar > &weights, Teuchos::Array< Scalar > &norms) const
     void minValue (Teuchos::Array< Scalar > &mins) const
     void maxValue (Teuchos::Array< Scalar > &maxs) const
     void meanValue (Teuchos::Array< Scalar > &means) const
     void multiply (Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const MultiVector< Ordinal, Scalar > &B, const Scalar &beta)
     void multiply (const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const MultiVector< Ordinal, Scalar > &B, const Scalar &beta)
     void reciprocalMultiply (const Scalar &alpha, const MultiVector< Ordinal, Scalar > &A, const MultiVector< Ordinal, Scalar > &B, const Scalar &beta)
  */

  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  #define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, basic, ORDINAL, double ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadConstNumVecs, ORDINAL, double ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadConstLDA, ORDINAL, double ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, OrthoDot, ORDINAL, double ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountDot, ORDINAL, double ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountNorm1, ORDINAL, double ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CountNormInf, ORDINAL, double ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, Norm2, ORDINAL, double )

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

      // FINISH: add complex tests

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, basic, ORDINAL, double )

    typedef short int ShortInt;
    UNIT_TEST_GROUP_ORDINAL(ShortInt)
    typedef long int LongInt;
    UNIT_TEST_GROUP_ORDINAL(LongInt)
#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
