// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_KokkosCounter.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <iterator>

// FINISH: add test for MultiVector with a node containing zero local entries
// FINISH: add tests for local MultiVectors


// Macro that marks a function as "possibly unused," in order to
// suppress build warnings.
#if ! defined(TRILINOS_UNUSED_FUNCTION)
#  if defined(__GNUC__) || (defined(__INTEL_COMPILER) && !defined(_MSC_VER))
#    define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#  elif defined(__clang__)
#    if __has_attribute(unused)
#      define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#    else
#      define TRILINOS_UNUSED_FUNCTION
#    endif // Clang has 'unused' attribute
#  elif defined(__IBMCPP__)
// IBM's C++ compiler for Blue Gene/Q (V12.1) implements 'used' but not 'unused'.
//
// http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp
#    define TRILINOS_UNUSED_FUNCTION
#  else // some other compiler
#    define TRILINOS_UNUSED_FUNCTION
#  endif
#endif // ! defined(TRILINOS_UNUSED_FUNCTION)


namespace Teuchos {
  template <>
  ScalarTraits<int>::magnitudeType
  relErr( const int &s1, const int &s2 ) {
    typedef ScalarTraits<int> ST;
    return ST::magnitude(s1-s2);
  }

  template <>
  ScalarTraits<char>::magnitudeType
  relErr( const char &s1, const char &s2 ) {
    typedef ScalarTraits<char> ST;
    return ST::magnitude(s1-s2);
  }
}

namespace {

  using Tpetra::TestingUtilities::getDefaultComm;

  using std::endl;
  using std::copy;
  using std::string;

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::arrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::null;
  using Teuchos::Range1D;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::OrdinalTraits;
  using Teuchos::outArg;
  using Teuchos::ScalarTraits;
  using Teuchos::SerialDenseMatrix;
  using Teuchos::Tuple;
  using Teuchos::tuple;
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
  using Tpetra::MultiVector;
  using Tpetra::global_size_t;
  using Tpetra::GloballyDistributed;

  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;

  double errorTolSlack = 1.0e+2;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  // no ScalarTraits<>::eps() for integer types
  template <class Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType testingTol() { return Teuchos::ScalarTraits<Scalar>::eps(); }
  template <>
  TRILINOS_UNUSED_FUNCTION int testingTol<int>() { return 0; }
  template <>
  TRILINOS_UNUSED_FUNCTION long testingTol<long>() { return 0; }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, NonMemberConstructors, LO, GO, Scalar , Node )
  {
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::Vector<Scalar,LO,GO,Node> V;
    constexpr bool debug = true;

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test: MultiVector, NonMemberConstructors" << endl;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Create a Map" << endl;
    auto comm = getDefaultComm ();
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid ();
    const size_t numLocal = 13;
    const size_t numVecs  = 7;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Create a MultiVector, and make sure that it has "
      "the right number of vectors (columns)" << endl;
    RCP<MV> mvec = Tpetra::createMultiVector<Scalar>(map,numVecs);
    TEST_EQUALITY(mvec->getNumVectors(), numVecs);

    myOut << "Create a Vector, and make sure that "
      "it has exactly one vector (column)" << endl;
    RCP<V> vec = Tpetra::createVector<Scalar>(map);
    TEST_EQUALITY_CONST(vec->getNumVectors(), 1);

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, basic, LO, GO, Scalar , Node )
  {
    using map_type = Tpetra::Map<LO, GO, Node>;
    using MV = Tpetra::MultiVector<Scalar, LO, GO, Node>;
    using vec_type = Tpetra::Vector<Scalar, LO, GO, Node>;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    constexpr bool debug = true;

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test: MultiVector, basic" << endl;
    Teuchos::OSTab tab0 (myOut);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid ();
    RCP<const Comm<int> > comm = getDefaultComm ();
    const int numImages = comm->getSize ();

    myOut << "Create Map" << endl;
    const size_t numLocal = 13;
    const size_t numVecs  = 7;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Test MultiVector's & Vector's default constructors" << endl;
    {
      MV defaultConstructedMultiVector;
      auto dcmv_map = defaultConstructedMultiVector.getMap ();
      TEST_ASSERT( dcmv_map.get () != nullptr );
      if (dcmv_map.get () != nullptr) {
        TEST_EQUALITY( dcmv_map->getGlobalNumElements (),
                       Tpetra::global_size_t (0) );
      }
      vec_type defaultConstructedVector;
      auto dcv_map = defaultConstructedVector.getMap ();
      TEST_ASSERT( dcv_map.get () != nullptr );
      if (dcv_map.get () != nullptr) {
        TEST_EQUALITY( dcv_map->getGlobalNumElements (),
                       Tpetra::global_size_t (0) );
      }
    }

    myOut << "Test MultiVector's usual constructor" << endl;
    RCP<MV> mvec;
    TEST_NOTHROW( mvec = rcp (new MV (map, numVecs, true)) );
    if (mvec.is_null ()) {
      myOut << "MV constructor threw an exception: returning" << endl;
      return;
    }
    TEST_EQUALITY( mvec->getNumVectors(), numVecs );
    TEST_EQUALITY( mvec->getLocalLength(), numLocal );
    TEST_EQUALITY( mvec->getGlobalLength(), numImages*numLocal );

    myOut << "Test that all norms are zero" << endl;
    Array<Magnitude> norms(numVecs), zeros(numVecs);
    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Magnitude>::zero());
    TEST_NOTHROW( mvec->norm2(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    TEST_NOTHROW( mvec->norm1(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    TEST_NOTHROW( mvec->normInf(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    // print it
    myOut << *mvec << endl;

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, large, LO, GO, Scalar , Node )
  {
    using map_type = Tpetra::Map<LO, GO, Node>;
    using MV = Tpetra::MultiVector<Scalar, LO, GO, Node>;
    using vec_type = Tpetra::Vector<Scalar, LO, GO, Node>;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    constexpr bool debug = true;

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test: MultiVector, basic" << endl;
    Teuchos::OSTab tab0 (myOut);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid ();
    RCP<const Comm<int> > comm = getDefaultComm ();
    const int numImages = comm->getSize ();

    myOut << "Create Map" << endl;
    const size_t numLocal = 15000;
    const size_t numVecs  = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Test MultiVector's & Vector's default constructors" << endl;
    {
      MV defaultConstructedMultiVector;
      auto dcmv_map = defaultConstructedMultiVector.getMap ();
      TEST_ASSERT( dcmv_map.get () != nullptr );
      if (dcmv_map.get () != nullptr) {
        TEST_EQUALITY( dcmv_map->getGlobalNumElements (),
                       Tpetra::global_size_t (0) );
      }
      vec_type defaultConstructedVector;
      auto dcv_map = defaultConstructedVector.getMap ();
      TEST_ASSERT( dcv_map.get () != nullptr );
      if (dcv_map.get () != nullptr) {
        TEST_EQUALITY( dcv_map->getGlobalNumElements (),
                       Tpetra::global_size_t (0) );
      }
    }

    myOut << "Test MultiVector's usual constructor" << endl;
    RCP<MV> mvec;
    TEST_NOTHROW( mvec = rcp (new MV (map, numVecs, true)) );
    if (mvec.is_null ()) {
      myOut << "MV constructor threw an exception: returning" << endl;
      return;
    }
    TEST_EQUALITY( mvec->getNumVectors(), numVecs );
    TEST_EQUALITY( mvec->getLocalLength(), numLocal );
    TEST_EQUALITY( mvec->getGlobalLength(), numImages*numLocal );

    myOut << "Test that all norms are zero" << endl;
    Array<Magnitude> norms(numVecs), zeros(numVecs);
    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Magnitude>::zero());
    TEST_NOTHROW( mvec->norm2(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    TEST_NOTHROW( mvec->norm1(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    TEST_NOTHROW( mvec->normInf(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    // print it
    myOut << *mvec << endl;

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, BadConstLDA, LO, GO, Scalar , Node )
  {
    // numlocal > LDA
    // ergo, the arrayview doesn't contain enough data to specify the entries
    // also, if bounds checking is enabled, check that bad bounds are caught
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

    out << "Test: MultiVector, BadConstLDA" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 2;
    const size_t numVecs = 2;
    // multivector has two vectors, each proc having two values per vector
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    // we need 4 scalars to specify values on each proc
    Array<Scalar> values(4);
#ifdef HAVE_TPETRA_DEBUG
    typedef Tpetra::Vector<Scalar,LO,GO,Node>       V;
    // too small an ArrayView (less than 4 values) is met with an exception, if debugging is on
    TEST_THROW(MV mvec(map,values(0,3),2,numVecs), std::invalid_argument);
    // it could also be too small for the given LDA:
    TEST_THROW(MV mvec(map,values(),2+1,numVecs), std::invalid_argument);
    // too small for number of entries in a Vector
    TEST_THROW(V   vec(map,values(0,1)), std::invalid_argument);
#endif
    // LDA < numLocal throws an exception anytime
    TEST_THROW(MV mvec(map,values(0,4),1,numVecs), std::runtime_error);

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, NonContigView, LO, GO, Scalar , Node )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::Vector<Scalar,LO,GO,Node> V;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;

    out << "Test: MultiVector, NonContigView: Scalar="
        << Teuchos::TypeNameTraits<Scalar>::name() << endl;
    Teuchos::OSTab tab0 (out);

    const Mag tol = errorTolSlack * errorTolSlack * testingTol<Scalar>();   // extra slack on this test; dots() seem to be a little sensitive for single precision types
    out << "tol: " << tol << endl;

    const Mag M0  = ScalarTraits<Mag>::zero();
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 53; // making this larger reduces the change that A below will have no non-zero entries, i.e., that C = abs(A) is still equal to A (we assume it is not)
    const size_t numVecs = 7;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    //
    // we will create a non-contig subview of the vector; un-viewed vectors should not be changed
    Tuple<size_t,4> inView1 = tuple<size_t>(1,4,3,2);
    Tuple<size_t,3> exView1 = tuple<size_t>(0,5,6);
    Tuple<size_t,4> inView2 = tuple<size_t>(6,0,4,3);
    Tuple<size_t,4> exView2 = tuple<size_t>(1,2,5,7);
    const size_t numView = 4;
    TEUCHOS_TEST_FOR_EXCEPTION(numView != as<size_t>(inView1.size()), std::logic_error, "Someone ruined a test invariant.");
    TEUCHOS_TEST_FOR_EXCEPTION(numView != as<size_t>(inView1.size()), std::logic_error, "Someone ruined a test invariant.");
    TEUCHOS_TEST_FOR_EXCEPTION(numView != as<size_t>(inView2.size()), std::logic_error, "Someone ruined a test invariant.");
    {
      // test dot, all norms, randomize
      MV mvOrig1(map,numVecs), mvOrig2(map,numVecs+1), mvWeights(map,numVecs);
      mvWeights.randomize();
      RCP<const MV> mvW1 = mvWeights.subView(tuple<size_t>(0));
      RCP<const MV> mvSubWeights = mvWeights.subView(inView1);
      mvOrig1.randomize();
      mvOrig2.randomize();
      //
      Array<Mag> nOrig2(numVecs), nOrig1(numVecs), nOrigI(numVecs);
      Array<Scalar> meansOrig(numVecs), dotsOrig(numView);
      mvOrig1.norm1(nOrig1());
      mvOrig1.norm2(nOrig2());
      mvOrig1.normInf(nOrigI());
      mvOrig1.meanValue(meansOrig());
      for (size_t j=0; j < numView; ++j) {
        RCP<const V> v1 = mvOrig1.getVector(inView1[j]),
          v2 = mvOrig2.getVector(inView2[j]);
        dotsOrig[j] = v1->dot(*v2);
      }
      // create the views, compute and test
      RCP<      MV> mvView1 = mvOrig1.subViewNonConst(inView1);
      TEST_ASSERT( ! mvView1.is_null() );
      if (! mvView1.is_null() ) {
        {
          auto mvView1_d = mvView1->getLocalViewDevice(Tpetra::Access::ReadOnly);
          auto mvOrig1_d = mvOrig1.getLocalViewDevice(Tpetra::Access::ReadOnly);
          TEST_ASSERT( mvView1_d.data() == mvOrig1_d.data() );
        }
        {
          auto mvView1_h = mvView1->getLocalViewHost(Tpetra::Access::ReadOnly);
          auto mvOrig1_h = mvOrig1.getLocalViewHost(Tpetra::Access::ReadOnly);
          TEST_ASSERT( mvView1_h.data() == mvOrig1_h.data() );
        }
      }
      RCP<const MV> mvView2 = mvOrig2.subView(inView2);
      TEST_ASSERT( ! mvView2.is_null() );
      if (! mvView2.is_null() ) {
        {
          auto mvView2_lcl = mvView2->getLocalViewDevice(Tpetra::Access::ReadOnly);
          auto mvOrig2_lcl = mvOrig2.getLocalViewDevice(Tpetra::Access::ReadOnly);
          TEST_ASSERT( mvView2_lcl.data() == mvOrig2_lcl.data() );
        }
        {
          auto mvView2_h = mvView2->getLocalViewHost(Tpetra::Access::ReadOnly);
          auto mvOrig2_h = mvOrig2.getLocalViewHost(Tpetra::Access::ReadOnly);
          TEST_ASSERT( mvView2_h.data() == mvOrig2_h.data() );
        }
      }
      Array<Mag> nView2(numView), nView1(numView), nViewI(numView);
      Array<Scalar> meansView(numView), dotsView(numView);
      mvView1->norm1(nView1());
      mvView1->norm2(nView2());
      mvView1->normInf(nViewI());
      mvView1->meanValue(meansView());
      mvView1->dot( *mvView2, dotsView() );

      for (size_t j=0; j < numView; ++j) {
        TEST_FLOATING_EQUALITY(nOrig1[inView1[j]],  nView1[j],  tol);
      }
      for (size_t j=0; j < numView; ++j) {
        TEST_FLOATING_EQUALITY(nOrig2[inView1[j]],  nView2[j],  tol);
      }
      for (size_t j=0; j < numView; ++j) {
        TEST_FLOATING_EQUALITY(nOrigI[inView1[j]],  nViewI[j],  tol);
      }
      for (size_t j=0; j < numView; ++j) {
        TEST_FLOATING_EQUALITY(meansOrig[inView1[j]], meansView[j], tol);
      }
      for (size_t j=0; j < numView; ++j) {
        TEST_FLOATING_EQUALITY(dotsOrig[j], dotsView[j], tol);
      }

      int lclSuccess = success ? 1 : 0;
      int gblSuccess = 0;
      using Teuchos::outArg;
      using Teuchos::reduceAll;
      using Teuchos::REDUCE_MIN;
      reduceAll(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (! success) {
        return;
      }

      // randomize the view, compute view one-norms, test difference
      mvView2 = Teuchos::null;
      mvView1->randomize();
      Array<Mag> nView1_aft(numView);
      mvView1->norm1(nView1_aft());
      for (size_t j=0; j < numView; ++j) {
        TEST_INEQUALITY(nView1[j], nView1_aft[j]);
      }
      // release the view, test that viewed columns changed, others didn't
      mvView1 = Teuchos::null;
      Array<Mag> nOrig1_aft(numVecs);
      mvOrig1.norm1(nOrig1_aft());
      for (size_t j=0; j < as<size_t>(inView1.size()); ++j) {
        TEST_INEQUALITY(nOrig1[inView1[j]], nOrig1_aft[inView1[j]]);
      }
      for (size_t j=0; j < as<size_t>(exView1.size()); ++j) {
        TEST_FLOATING_EQUALITY(nOrig1[exView1[j]], nOrig1_aft[exView1[j]], tol);
      }
    }
    {
      MV mvOrigA(map,numVecs), mvOrigB(map,numVecs), mvOrigC(map,numVecs+1);
      // we don't know what the distribution is, so that we don't know that abs(A) != A
      // therefore, do some manipulation
      mvOrigA.randomize();
      mvOrigB.randomize();
      // A = A - B, then new (indendent) values for B
      mvOrigA.update(as<Scalar>(-1),mvOrigB, as<Scalar>(1));
      mvOrigB.randomize();
      mvOrigC.randomize();
      Array<Mag> nrmOrigA(numVecs), nrmOrigB(numVecs), nrmOrigC(numVecs+1);
      mvOrigA.norm2(nrmOrigA());
      mvOrigB.norm2(nrmOrigB());
      mvOrigC.norm2(nrmOrigC());
      RCP<MV> mvViewA = mvOrigA.subViewNonConst(inView1);
      RCP<MV> mvViewB = mvOrigB.subViewNonConst(inView1);
      RCP<MV> mvViewC = mvOrigC.subViewNonConst(inView2);
      // set C = abs(A)
      {
        Array<Scalar> mnA_bef(inView1.size()), mnC_bef(inView1.size()),
                      mnA_aft(inView1.size()), mnC_aft(inView1.size());
        mvViewA->meanValue(mnA_bef());
        mvViewC->meanValue(mnC_bef());
        mvViewC->abs(*mvViewA);
        mvViewA->meanValue(mnA_aft());
        mvViewC->meanValue(mnC_aft());
        for (size_t j=0; j < as<size_t>(inView1.size()); ++j) {
          TEST_FLOATING_EQUALITY(mnA_bef[j], mnA_aft[j], tol);
          TEST_INEQUALITY(mnC_bef[j], mnC_aft[j]);
        }
      }
      // then set A = B = C
      // good excuse for some double views
      // use full views of C and B for this, check means before and after
      // to make sure that only A and B change.
      {
        Array<Scalar> A_bef(inView1.size()), B_bef(inView1.size()), C_bef(inView2.size());
        mvViewA->meanValue(A_bef());
        mvViewB->meanValue(B_bef());
        mvViewC->meanValue(C_bef());
        RCP<MV> doubleViewA = mvViewA->subViewNonConst(Range1D(0,inView1.size()-1));
        RCP<MV> doubleViewB = mvViewB->subViewNonConst(Range1D(0,inView1.size()-1));
        RCP<const MV> doubleViewC = mvViewC->subView(Range1D(0,inView1.size()-1));
        //(*doubleViewA) = (*doubleViewB) = (*doubleViewC);
        deep_copy((*doubleViewB),(*doubleViewC));
        deep_copy((*doubleViewA),(*doubleViewB));
        doubleViewA = Teuchos::null;
        doubleViewB = Teuchos::null;
        doubleViewC = Teuchos::null;
        Array<Scalar> A_aft(inView1.size()), B_aft(inView1.size()), C_aft(inView2.size());
        mvViewA->meanValue(A_aft());
        mvViewB->meanValue(B_aft());
        mvViewC->meanValue(C_aft());
        for (size_t j=0; j < as<size_t>(inView1.size()); ++j) {
          TEST_FLOATING_EQUALITY(C_bef[j], C_aft[j], tol);
          TEST_FLOATING_EQUALITY(C_bef[j], B_aft[j], tol);
          TEST_FLOATING_EQUALITY(C_bef[j], A_aft[j], tol);
          TEST_INEQUALITY(A_bef[j], A_aft[j]);
          TEST_INEQUALITY(B_bef[j], B_aft[j]);
        }
      }
      {
        TEUCHOS_TEST_FOR_EXCEPTION(inView1.size() != 4, std::logic_error, "Someone ruined a test invariant.");
        Tuple<size_t,4> reorder = tuple<size_t>(3,1,0,2);
        RCP<MV> dvA = mvViewA->subViewNonConst(reorder);
        RCP<MV> dvB = mvViewB->subViewNonConst(reorder);
        RCP<MV> dvC = mvViewC->subViewNonConst(reorder);
        // C == B == A
        //   C *= 2                ->  C == 2*A == 2*B            scale(alpha)
        dvC->scale( as<Scalar>(2) );
        //   A = -C + 2*A          ->  C == 2*B, A == 0           update(alpha,mv,beta)
        dvA->update(as<Scalar>(-1),*dvC, as<Scalar>(2));
        //   C = 2*A + 2*B - .5*C ->   C == B, A == 0,            update(alpha,mv,beta,mv,gamma)
        dvC->update(as<Scalar>(2),*dvA, as<Scalar>(2), *dvB, as<Scalar>(-.5));
        //   B = 0.5              ->   B = 0.5, A == 0,           putScalar(alpha)
        dvB->putScalar( as<Scalar>(0.5) );
        //   C.recip(B)           ->   C = 2, B == 0.5, A == 0,   reciprocal(mv)
        dvC->reciprocal(*dvB);
        //   B = C/2              ->   A == 0, B == 1, C == 2
        dvB->scale(as<Mag>(0.5),*dvC);
        dvA = Teuchos::null;
        dvB = Teuchos::null;
        dvC = Teuchos::null;
        Array<Mag> nrmA(4), nrmB(4), nrmC(4);
        mvViewA->norm1(nrmA()); // norm1(0)   = 0
        mvViewB->norm1(nrmB()); // norm1(1.0) = N
        mvViewC->norm1(nrmC()); // norm1(2.0) = 2 * N
        const Mag  OneN = as<Mag>(mvViewA->getGlobalLength());
        const Mag  TwoN = OneN + OneN;
        for (size_t j=0; j < 4; ++j) {
          TEST_FLOATING_EQUALITY( nrmA[j],    M0, tol );
          TEST_FLOATING_EQUALITY( nrmB[j],  OneN, tol );
          TEST_FLOATING_EQUALITY( nrmC[j],  TwoN, tol );
        }
      }
      // done with these views; clear them, ensure that only the viewed
      // vectors changed in the original multivectors
      mvViewA = Teuchos::null;
      mvViewB = Teuchos::null;
      mvViewC = Teuchos::null;
      Array<Mag> nrmOrigA_aft(numVecs), nrmOrigB_aft(numVecs), nrmOrigC_aft(numVecs+1);
      mvOrigA.norm2(nrmOrigA_aft());
      mvOrigB.norm2(nrmOrigB_aft());
      mvOrigC.norm2(nrmOrigC_aft());
      for (size_t j=0; j < as<size_t>(inView1.size()); ++j) {
        TEST_INEQUALITY(nrmOrigA[inView1[j]], nrmOrigA_aft[inView1[j]]);
        TEST_INEQUALITY(nrmOrigB[inView1[j]], nrmOrigB_aft[inView1[j]]);
        TEST_INEQUALITY(nrmOrigC[inView2[j]], nrmOrigC_aft[inView2[j]]);
      }
      for (size_t j=0; j < as<size_t>(exView1.size()); ++j) {
        TEST_FLOATING_EQUALITY(nrmOrigA[exView1[j]], nrmOrigA_aft[exView1[j]], tol);
        TEST_FLOATING_EQUALITY(nrmOrigB[exView1[j]], nrmOrigB_aft[exView1[j]], tol);
      }
      for (size_t j=0; j < as<size_t>(exView1.size()); ++j) {
        TEST_FLOATING_EQUALITY(nrmOrigC[exView2[j]], nrmOrigC_aft[exView2[j]], tol);
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, Describable, LO , GO , Scalar, Node )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

    out << "Test: MultiVector, Describable" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // create Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,3,comm);
    // test labeling
    const string lbl("mvecA");
    MV mvecA(map,2);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, BadMultiply, LO , GO , Scalar , Node )
  {
    // mfh 05 May 2016: Tpetra::MultiVector::multiply only checks
    // local dimensions in a debug build.
#ifdef HAVE_TPETRA_DEBUG
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

    out << "Test: MultiVector, BadMultiply" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const Scalar S1 = ScalarTraits<Scalar>::one(),
                 S0 = ScalarTraits<Scalar>::zero();
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    {
      // create local Maps
      RCP<const Map<LO,GO,Node> > map3l = createLocalMapWithNode<LO,GO,Node>(3,comm),
                                  map2l = createLocalMapWithNode<LO,GO,Node>(2,comm);
      MV mvecA(map3l,2),
         mvecB(map2l,3),
         mvecD(map2l,2);
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
      RCP<const Map<LO,GO,Node> > map3n = createContigMapWithNode<LO,GO,Node>(INVALID,3,comm);
      RCP<const Map<LO,GO,Node> > map2n = createContigMapWithNode<LO,GO,Node>(INVALID,2,comm);
      RCP<const Map<LO,GO,Node> > map2l = createLocalMapWithNode<LO,GO,Node>(2,comm),
                                  map3l = createLocalMapWithNode<LO,GO,Node>(3,comm);
      MV mv3nx2(map3n,2),
         mv2nx2(map2n,2),
         mv2lx2(map2l,2),
         mv2lx3(map2l,3),
         mv3lx2(map3l,2),
         mv3lx3(map3l,3);
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
      RCP<const Map<LO,GO,Node> > map3n = createContigMapWithNode<LO,GO,Node>(INVALID,3,comm);
      RCP<const Map<LO,GO,Node> > map2n = createContigMapWithNode<LO,GO,Node>(INVALID,2,comm);
      RCP<const Map<LO,GO,Node> > map2l = createLocalMapWithNode<LO,GO,Node>(2,comm),
                                  map3l = createLocalMapWithNode<LO,GO,Node>(3,comm);
      MV mv3nx2(map3n,2),
         mv2nx2(map2n,2),
         mv2x3(map2l,3),
         mv3x2(map3l,2);
      // non-matching input lengths
      TEST_THROW( mv3nx2.multiply(NO_TRANS,CONJ_TRANS,S1,mv3nx2,mv2x3,S0), std::runtime_error);   // (3n x 2) x (3 x 2) (trans) not compat
      TEST_THROW( mv3nx2.multiply(NO_TRANS,NO_TRANS  ,S1,mv3nx2,mv3x2,S0), std::runtime_error);   // (3n x 2) x (3 x 2) (nontrans) not compat
      // non-matching output sizes
      TEST_THROW( mv3nx2.multiply(NO_TRANS,CONJ_TRANS,S1,mv3nx2,mv3x2,S0), std::runtime_error);   // (3n x 2) x (2 x 3) doesn't fit 3nx2
      TEST_THROW( mv3nx2.multiply(NO_TRANS,NO_TRANS  ,S1,mv3nx2,mv2x3,S0), std::runtime_error);   // (3n x 2) x (2 x 3) doesn't fit 3nx2
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
#endif // HAVE_TPETRA_DEBUG
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, Multiply, LO , GO , Scalar , Node )
  {
    using Teuchos::View;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    int lclSuccess = 1;
    int gblSuccess = 0;

    out << "Test: MultiVector, Multiply" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    RCP<const Map<LO,GO,Node> > map3n = createContigMapWithNode<LO,GO,Node>(INVALID,3,comm),
                                map2n = createContigMapWithNode<LO,GO,Node>(INVALID,2,comm);
    RCP<const Map<LO,GO,Node> > lmap3 = createLocalMapWithNode<LO,GO,Node>(3,comm),
                                lmap2 = createLocalMapWithNode<LO,GO,Node>(2,comm);
    const Scalar S1 = ScalarTraits<Scalar>::one(),
                 S0 = ScalarTraits<Scalar>::zero();
    const Mag    M0 = ScalarTraits<Mag>::zero();
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    // deterministic input/output
    {
      MV mv3x2l(lmap3,2),
         mv2x3l(lmap2,3),
         mv2x2l(lmap2,2),
         mv3x3l(lmap3,3);
      // fill multivectors with ones
      mv3x2l.putScalar(ScalarTraits<Scalar>::one());
      mv2x3l.putScalar(ScalarTraits<Scalar>::one());
      // fill expected answers Array
      Teuchos::Array<Scalar> check2(4,3); // each entry (of four) is the product [1 1 1]*[1 1 1]' = 3
      Teuchos::Array<Scalar> check3(9,2); // each entry (of nine) is the product [1 1]*[1 1]' = 2
      // test
      mv3x3l.multiply(NO_TRANS  ,NO_TRANS  ,S1,mv3x2l,mv2x3l,S0);
      { auto tmpView = mv3x3l.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView(0,9),check3,M0); }
      mv2x2l.multiply(NO_TRANS  ,CONJ_TRANS,S1,mv2x3l,mv2x3l,S0);
      { auto tmpView = mv2x2l.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView(0,4),check2,M0); }
      mv2x2l.multiply(CONJ_TRANS,NO_TRANS  ,S1,mv3x2l,mv3x2l,S0);
      { auto tmpView = mv2x2l.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView(0,4),check2,M0); }
      mv3x3l.multiply(CONJ_TRANS,CONJ_TRANS,S1,mv2x3l,mv3x2l,S0);
      { auto tmpView = mv3x3l.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView(0,9),check3,M0); }
    }
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                         outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );

    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    // random input/output
    {
      Array<Scalar>     tmvCopy1(6), tmvCopy2(6);
      ArrayView<Scalar> sdmView(Teuchos::null);
      MV tmv3x2(lmap3,2),
         tmv2x3(lmap2,3),
         tmv2x2(lmap2,2),
         tmv3x3(lmap3,3);
      // fill multivectors with random, get copy of contents
      tmv3x2.randomize();  tmv3x2.get1dCopy(tmvCopy1(),3);
      tmv2x3.randomize();  tmv2x3.get1dCopy(tmvCopy2(),2);
      // point SerialDenseMatrices at copies
      SerialDenseMatrix<int,Scalar> sdm3x2(View,tmvCopy1.getRawPtr(),3,3,2);
      SerialDenseMatrix<int,Scalar> sdm2x3(View,tmvCopy2.getRawPtr(),2,2,3);
      // space for answers
      SerialDenseMatrix<int,Scalar> sdm2x2(2,2), sdm3x3(3,3);
      // test: perform local Tpetra::MultiVector multiply and Teuchos::SerialDenseMatrix multiply, then check that answers are equivalent
      {
        tmv3x3.multiply(NO_TRANS,NO_TRANS,S1,tmv3x2,tmv2x3,S0);
        sdm3x3.multiply(NO_TRANS,NO_TRANS,S1,sdm3x2,sdm2x3,S0);
        {
          auto tmpView = tmv3x3.get1dView(); sdmView = arrayView(sdm3x3.values(),sdm3x3.numRows()*sdm3x3.numCols());
          TEST_COMPARE_FLOATING_ARRAYS(tmpView,sdmView,testingTol<Scalar>() * errorTolSlack);
        }
      }
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                           outArg (gblSuccess));
      TEST_ASSERT( gblSuccess == 1 );

      {
        tmv2x2.multiply(NO_TRANS,CONJ_TRANS,S1,tmv2x3,tmv2x3,S0);
        sdm2x2.multiply(NO_TRANS,CONJ_TRANS,S1,sdm2x3,sdm2x3,S0);
        { 
          auto tmpView = tmv2x2.get1dView(); sdmView = arrayView(sdm2x2.values(),sdm2x2.numRows()*sdm2x2.numCols()); 
          TEST_COMPARE_FLOATING_ARRAYS(tmpView,sdmView,testingTol<Scalar>() * errorTolSlack);
        }
      }
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                           outArg (gblSuccess));
      TEST_ASSERT( gblSuccess == 1 );

      {
        tmv2x2.multiply(CONJ_TRANS,NO_TRANS,S1,tmv3x2,tmv3x2,S0);
        Kokkos::fence ();
        sdm2x2.multiply(CONJ_TRANS,NO_TRANS,S1,sdm3x2,sdm3x2,S0);
        { 
          auto tmpView = tmv2x2.get1dView(); sdmView = arrayView(sdm2x2.values(),sdm2x2.numRows()*sdm2x2.numCols()); 
          TEST_COMPARE_FLOATING_ARRAYS(tmpView,sdmView,testingTol<Scalar>() * errorTolSlack);
        }
      }
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                           outArg (gblSuccess));
      TEST_ASSERT( gblSuccess == 1 );

      {
        tmv3x3.multiply(CONJ_TRANS,CONJ_TRANS,S1,tmv2x3,tmv3x2,S0);
        Kokkos::fence ();
        sdm3x3.multiply(CONJ_TRANS,CONJ_TRANS,S1,sdm2x3,sdm3x2,S0);
        { 
          auto tmpView = tmv3x3.get1dView(); sdmView = arrayView(sdm3x3.values(),sdm3x3.numRows()*sdm3x3.numCols()); 
          TEST_COMPARE_FLOATING_ARRAYS(tmpView,sdmView,testingTol<Scalar>() * errorTolSlack);
        }
      }
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                           outArg (gblSuccess));
      TEST_ASSERT( gblSuccess == 1 );
    }

    // case 2: C(local) = A^T(distr) * B  (distr)  : one of these
    {
      MV mv3nx2(map3n,2),
         mv3nx3(map3n,3),
         // locals
         mv2x2(lmap2,2),
         mv2x3(lmap2,3),
         mv3x2(lmap3,2),
         mv3x3(lmap3,3);
      // fill multivectors with ones
      mv3nx3.putScalar(ScalarTraits<Scalar>::one());
      mv3nx2.putScalar(ScalarTraits<Scalar>::one());
      // fill expected answers Array

      Teuchos::Array<Scalar> check(9,3*numImages);
      // test
      mv2x2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx2,S0);
      { auto tmpView = mv2x2.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0); }
      mv2x3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx2,mv3nx3,S0);
      { auto tmpView = mv2x3.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0); }
      mv3x2.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx3,mv3nx2,S0);
      { auto tmpView = mv3x2.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);}
      mv3x3.multiply(CONJ_TRANS,NO_TRANS,S1,mv3nx3,mv3nx3,S0);
      { auto tmpView = mv3x3.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check(0,tmpView.size()),M0);}
    }
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                         outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );

    // case 3: C(distr) = A  (distr) * B^X(local)  : two of these
    {
      MV mv3nx2(map3n,2),
         mv3nx3(map3n,3),
         // locals
         mv2x3(lmap2,3);
      // fill multivectors with ones
      mv2x3.putScalar(S1);
      // fill expected answers Array
      Teuchos::Array<Scalar> check2(9,2), check3(6,3);
      // test
      mv3nx3.putScalar(S1); mv3nx2.putScalar(S1);
      mv3nx3.multiply(NO_TRANS,  NO_TRANS,S1,mv3nx2,mv2x3,S0);
      { auto tmpView = mv3nx3.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check2,M0);}
      mv3nx3.putScalar(S1); mv3nx2.putScalar(S1);
      mv3nx2.multiply(NO_TRANS,CONJ_TRANS,S1,mv3nx3,mv2x3,S0);
      { auto tmpView = mv3nx2.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check3,M0);}
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                         outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }

  // Test Tpetra::MultiVector::elementWiseMultiply.
  //
  // Be sure to exercise all combinations of the cases alpha =
  // {-1,0,1,other} and beta = {-1,0,1,other}, as these commonly have
  // special cases.
  //
  // Also be sure to exercise the common case (also often with a
  // special-case implementation) where all the MultiVectors have one
  // column.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, ElementWiseMultiply, LO , GO , ST , Node )
  {
    using Teuchos::View;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<ST> STS;
    typedef typename STS::magnitudeType MT;
    typedef Teuchos::ScalarTraits<MT> STM;
    typedef Tpetra::Map<LO,GO,Node> map_type;
    typedef Tpetra::MultiVector<ST,LO,GO,Node> MV;
    typedef Tpetra::Vector<ST,LO,GO,Node> V;
    typedef typename Kokkos::ArithTraits<ST>::val_type IST;

    out << "Tpetra::MultiVector::elementWiseMultiply test" << endl;
    Teuchos::OSTab tab0 (out);

    // Create a Map.
    RCP<const Comm<int> > comm = getDefaultComm ();
    const size_t lclNumRows = 3;
    const GST gblNumRows = comm->getSize () * lclNumRows;
    const GO indexBase = 0;
    RCP<const map_type> map3n =
      rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

    const MT M0 = STM::zero ();
    const ST S0 = STS::zero ();
    const ST S1 = STS::one ();

    // In what follows, '@' (without single quotes) denotes
    // element-wise multiplication -- that is, what
    // MultiVector::elementWiseMultiply implements.

    const size_t maxNumVecs = 3;

    // Test for various numbers of columns.
    for (size_t numVecs = 1; numVecs <= maxNumVecs; ++numVecs) {
      out << "Test numVecs = " << numVecs << endl;
      Teuchos::OSTab tab1 (out);

      // A (always) has 1 vector, and B and C have numVecs vectors.
      V A (map3n);
      MV B (map3n, numVecs);
      MV C (map3n, numVecs);
      MV C_exp (map3n, numVecs);
      Array<MT> C_norms (C.getNumVectors ());
      Array<MT> C_zeros (C.getNumVectors ());
      std::fill (C_zeros.begin (), C_zeros.end (), M0);

      int caseNum = 0;

      caseNum++;
      out << "Case " << caseNum << ": C = 0*C + 0*(A @ B)" << endl;
      // Fill A and B initially with nonzero values, just for
      // generality.  C should get filled with zeros afterwards.
      // Prefill C with NaN, to ensure that the method follows BLAS
      // update rules.
      {
        A.putScalar (S1);
        B.putScalar (S1);

        // Prefill C with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        C.putScalar (nan);

        C.elementWiseMultiply (S0, A, B, S0);

        C_exp.putScalar (S0);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 1*C + 0*(A @ B)" << endl;
      // Fill A and B with NaN to check that the method follows BLAS
      // update rules.
      {
        const ST S3 = S1 + S1 + S1;

        // Prefill A and B with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        A.putScalar (nan);
        B.putScalar (nan);
        C.putScalar (S3);

        C.elementWiseMultiply (S0, A, B, S1);

        C_exp.putScalar (S3);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = (-1)*C + 0*(A @ B)" << endl;
      // Fill A and B with NaN to check that the method follows BLAS
      // update rules.
      {
        // Prefill A and B with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        A.putScalar (nan);
        B.putScalar (nan);
        C.putScalar (S1);

        C.elementWiseMultiply (S0, A, B, -S1);

        C_exp.putScalar (-S1);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 2*C + 0*(A @ B)" << endl;
      // Fill A and B with NaN to check that the method follows BLAS
      // update rules.
      {
        const ST S2 = S1 + S1;

        // Prefill A and B with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        A.putScalar (nan);
        B.putScalar (nan);
        C.putScalar (S1);

        C.elementWiseMultiply (S0, A, B, S2);

        C_exp.putScalar (S2);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 0*C + 1*(A @ B)" << endl;
      // A and B will be filled with 1s, so C should get filled with 1s.
      // Prefill C with NaN, to ensure that the method follows BLAS
      // update rules.
      {
        A.putScalar (S1);
        B.putScalar (S1);

        // Prefill C with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        C.putScalar (nan);

        C.elementWiseMultiply (S1, A, B, S0);

        C_exp.putScalar (S1);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 0*C + (-1)*(A @ B)" << endl;
      // A and B will be filled with 1, so C should get filled with -1.
      // Prefill C with NaN, to ensure that the method follows BLAS
      // update rules.
      {
        A.putScalar (S1);
        B.putScalar (S1);

        // Prefill C with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        C.putScalar (nan);

        C.elementWiseMultiply (-S1, A, B, S0);

        C_exp.putScalar (-S1);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 1*C + 1*(A @ B)" << endl;
      // Fill A with 1, B with 2, and C with 3.  C should be 5 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S1 + S1 + S1;
        const ST S5 = S2 + S3;
        A.putScalar (S1);
        B.putScalar (S2);
        C.putScalar (S3);

        C.elementWiseMultiply (S1, A, B, S1);

        C_exp.putScalar (S5);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = (-1)*C + 1*(A @ B)" << endl;
      // Fill A with 1, B with 2, and C with 3.  C should be -1 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S1 + S1 + S1;
        A.putScalar (S1);
        B.putScalar (S2);
        C.putScalar (S3);

        C.elementWiseMultiply (S1, A, B, -S1);

        C_exp.putScalar (-S1);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 1*C + (-1)*(A @ B)" << endl;
      // Fill A with 2, B with 3, and C with 1.  C should be -5 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S2 + S1;
        const ST S5 = S2 + S3;

        A.putScalar (S2);
        B.putScalar (S3);
        C.putScalar (S1);

        C.elementWiseMultiply (-S1, A, B, S1);

        C_exp.putScalar (-S5);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = (-1)*C + (-1)*(A @ B)" << endl;
      // Fill A with 1, B with 2, and C with 3.  C should be -5 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S1 + S1 + S1;
        const ST S5 = S2 + S3;
        A.putScalar (S1);
        B.putScalar (S2);
        C.putScalar (S3);

        C.elementWiseMultiply (-S1, A, B, -S1);

        C_exp.putScalar (-S5);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 0*C + 2*(A @ B)" << endl;
      // Fill A with 3 and B with 4.  C should be 24 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S2 + S1;
        const ST S4 = S3 + S1;
        const ST S24 = S2 * S3 * S4;

        A.putScalar (S3);
        B.putScalar (S4);

        // Prefill C with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        C.putScalar (nan);

        C.elementWiseMultiply (S2, A, B, S0);

        C_exp.putScalar (S24);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = (-2)*C + 2*(A @ B)" << endl;
      // Fill A with 3, B with 4, and C with 5.  C should be 14 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S2 + S1;
        const ST S4 = S3 + S1;
        const ST S5 = S4 + S1;
        const ST S14 = S5 * S2 + S4;

        A.putScalar (S3);
        B.putScalar (S4);
        C.putScalar (S5);

        C.elementWiseMultiply (S2, A, B, -S2);

        C_exp.putScalar (S14);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }

  // Test Tpetra::MultiVector::elementWiseMultiply on a large multivector
  //
  // Be sure to exercise all combinations of the cases alpha =
  // {-1,0,1,other} and beta = {-1,0,1,other}, as these commonly have
  // special cases.
  //
  // Also be sure to exercise the common case (also often with a
  // special-case implementation) where all the MultiVectors have one
  // column.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, ElementWiseMultiplyLg, LO , GO , ST , Node )
  {
    using Teuchos::View;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<ST> STS;
    typedef typename STS::magnitudeType MT;
    typedef Teuchos::ScalarTraits<MT> STM;
    typedef Tpetra::Map<LO,GO,Node> map_type;
    typedef Tpetra::MultiVector<ST,LO,GO,Node> MV;
    typedef Tpetra::Vector<ST,LO,GO,Node> V;
    typedef typename Kokkos::ArithTraits<ST>::val_type IST;

    out << "Tpetra::MultiVector::elementWiseMultiplyLg test" << endl;
    Teuchos::OSTab tab0 (out);

    // Create a Map.
    RCP<const Comm<int> > comm = getDefaultComm ();
    const size_t lclNumRows = 15000;
    const GST gblNumRows = comm->getSize () * lclNumRows;
    const GO indexBase = 0;
    RCP<const map_type> map3n =
      rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

    const MT M0 = STM::zero ();
    const ST S0 = STS::zero ();
    const ST S1 = STS::one ();

    // In what follows, '@' (without single quotes) denotes
    // element-wise multiplication -- that is, what
    // MultiVector::elementWiseMultiply implements.

    const size_t maxNumVecs = 3;

    // Test for various numbers of columns.
    for (size_t numVecs = 1; numVecs <= maxNumVecs; ++numVecs) {
      out << "Test numVecs = " << numVecs << endl;
      Teuchos::OSTab tab1 (out);

      // A (always) has 1 vector, and B and C have numVecs vectors.
      V A (map3n);
      MV B (map3n, numVecs);
      MV C (map3n, numVecs);
      MV C_exp (map3n, numVecs);
      Array<MT> C_norms (C.getNumVectors ());
      Array<MT> C_zeros (C.getNumVectors ());
      std::fill (C_zeros.begin (), C_zeros.end (), M0);

      int caseNum = 0;

      caseNum++;
      out << "Case " << caseNum << ": C = 0*C + 0*(A @ B)" << endl;
      // Fill A and B initially with nonzero values, just for
      // generality.  C should get filled with zeros afterwards.
      // Prefill C with NaN, to ensure that the method follows BLAS
      // update rules.
      {
        A.putScalar (S1);
        B.putScalar (S1);

        // Prefill C with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        C.putScalar (nan);

        C.elementWiseMultiply (S0, A, B, S0);

        C_exp.putScalar (S0);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 1*C + 0*(A @ B)" << endl;
      // Fill A and B with NaN to check that the method follows BLAS
      // update rules.
      {
        const ST S3 = S1 + S1 + S1;

        // Prefill A and B with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        A.putScalar (nan);
        B.putScalar (nan);
        C.putScalar (S3);

        C.elementWiseMultiply (S0, A, B, S1);

        C_exp.putScalar (S3);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = (-1)*C + 0*(A @ B)" << endl;
      // Fill A and B with NaN to check that the method follows BLAS
      // update rules.
      {
        // Prefill A and B with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        A.putScalar (nan);
        B.putScalar (nan);
        C.putScalar (S1);

        C.elementWiseMultiply (S0, A, B, -S1);

        C_exp.putScalar (-S1);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 2*C + 0*(A @ B)" << endl;
      // Fill A and B with NaN to check that the method follows BLAS
      // update rules.
      {
        const ST S2 = S1 + S1;

        // Prefill A and B with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        A.putScalar (nan);
        B.putScalar (nan);
        C.putScalar (S1);

        C.elementWiseMultiply (S0, A, B, S2);

        C_exp.putScalar (S2);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 0*C + 1*(A @ B)" << endl;
      // A and B will be filled with 1s, so C should get filled with 1s.
      // Prefill C with NaN, to ensure that the method follows BLAS
      // update rules.
      {
        A.putScalar (S1);
        B.putScalar (S1);

        // Prefill C with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        C.putScalar (nan);

        C.elementWiseMultiply (S1, A, B, S0);

        C_exp.putScalar (S1);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 0*C + (-1)*(A @ B)" << endl;
      // A and B will be filled with 1, so C should get filled with -1.
      // Prefill C with NaN, to ensure that the method follows BLAS
      // update rules.
      {
        A.putScalar (S1);
        B.putScalar (S1);

        // Prefill C with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        C.putScalar (nan);

        C.elementWiseMultiply (-S1, A, B, S0);

        C_exp.putScalar (-S1);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 1*C + 1*(A @ B)" << endl;
      // Fill A with 1, B with 2, and C with 3.  C should be 5 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S1 + S1 + S1;
        const ST S5 = S2 + S3;
        A.putScalar (S1);
        B.putScalar (S2);
        C.putScalar (S3);

        C.elementWiseMultiply (S1, A, B, S1);

        C_exp.putScalar (S5);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = (-1)*C + 1*(A @ B)" << endl;
      // Fill A with 1, B with 2, and C with 3.  C should be -1 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S1 + S1 + S1;
        A.putScalar (S1);
        B.putScalar (S2);
        C.putScalar (S3);

        C.elementWiseMultiply (S1, A, B, -S1);

        C_exp.putScalar (-S1);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 1*C + (-1)*(A @ B)" << endl;
      // Fill A with 2, B with 3, and C with 1.  C should be -5 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S2 + S1;
        const ST S5 = S2 + S3;

        A.putScalar (S2);
        B.putScalar (S3);
        C.putScalar (S1);

        C.elementWiseMultiply (-S1, A, B, S1);

        C_exp.putScalar (-S5);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = (-1)*C + (-1)*(A @ B)" << endl;
      // Fill A with 1, B with 2, and C with 3.  C should be -5 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S1 + S1 + S1;
        const ST S5 = S2 + S3;
        A.putScalar (S1);
        B.putScalar (S2);
        C.putScalar (S3);

        C.elementWiseMultiply (-S1, A, B, -S1);

        C_exp.putScalar (-S5);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = 0*C + 2*(A @ B)" << endl;
      // Fill A with 3 and B with 4.  C should be 24 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S2 + S1;
        const ST S4 = S3 + S1;
        const ST S24 = S2 * S3 * S4;

        A.putScalar (S3);
        B.putScalar (S4);

        // Prefill C with NaN, if NaN exists for ST.
        const ST nan = static_cast<ST> (Kokkos::ArithTraits<IST>::nan ());
        C.putScalar (nan);

        C.elementWiseMultiply (S2, A, B, S0);

        C_exp.putScalar (S24);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }

      caseNum++;
      out << "Case " << caseNum << ": C = (-2)*C + 2*(A @ B)" << endl;
      // Fill A with 3, B with 4, and C with 5.  C should be 14 after.
      {
        const ST S2 = S1 + S1;
        const ST S3 = S2 + S1;
        const ST S4 = S3 + S1;
        const ST S5 = S4 + S1;
        const ST S14 = S5 * S2 + S4;

        A.putScalar (S3);
        B.putScalar (S4);
        C.putScalar (S5);

        C.elementWiseMultiply (S2, A, B, -S2);

        C_exp.putScalar (S14);
        C_exp.update (S1, C, -S1);
        C_exp.normInf (C_norms ());

        TEST_COMPARE_FLOATING_ARRAYS( C_norms, C_zeros, M0 );
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, BadConstAA, LO , GO , Scalar , Node )
  {
    // constructor takes ArrayView<ArrayView<Scalar> A, NumVectors
    // A.size() == NumVectors
    // A[i].size() >= MyLength
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

    out << "Test: MultiVector, BadConstAA" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    // multivector has two vectors, each proc having two values per vector
    RCP<const Map<LO,GO,Node> > map2 = createContigMapWithNode<LO,GO,Node>(INVALID,2,comm);
    RCP<const Map<LO,GO,Node> > map3 = createContigMapWithNode<LO,GO,Node>(INVALID,3,comm);
    // we need 4 scalars to specify values on each proc
    Array<Scalar> values(4);
    Array<ArrayView<const Scalar> > arrOfarr(2,ArrayView<const Scalar>(Teuchos::null));
    Array<ArrayView<const Scalar> > emptyArr;
    arrOfarr[0] = values(0,2);
    arrOfarr[1] = values(2,2);
    // arrOfarr.size() == 0
    TEST_THROW(MV mvec(map2,emptyArr(),0), std::runtime_error);
#ifdef HAVE_TPETRA_DEBUG
    // individual ArrayViews could be too small
    TEST_THROW(MV mvec(map3,arrOfarr(),2), std::runtime_error);
#endif

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, BadDot, LO , GO , Scalar , Node )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::Vector<Scalar,LO,GO,Node>       V;

    out << "Test: MultiVector, BadDot" << endl;
    Teuchos::OSTab tab0 (out);

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Map<LO,GO,Node> > map1 = createContigMapWithNode<LO,GO,Node>(INVALID,1,comm),
                                map2 = createContigMapWithNode<LO,GO,Node>(INVALID,2,comm);
    {
      MV mv12(map1,1),
         mv21(map2,1),
         mv22(map2,2);
      Array<Scalar> dots(2);
      // incompatible maps
      TEST_THROW(mv12.dot(mv21,dots()),std::runtime_error);
      // incompatible numvecs
      TEST_THROW(mv22.dot(mv21,dots()),std::runtime_error);
      // too small output array
#ifdef HAVE_TPETRA_DEBUG
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
#ifdef HAVE_TPETRA_DEBUG
      TEST_THROW(v1.dot(v2,dots()),std::runtime_error);
      TEST_THROW(v2.dot(v1,dots()),std::runtime_error);
#endif
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, OrthoDot, LO , GO , Scalar , Node )
  {
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

    out << "Test: MultiVector, OrthoDot" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar S0 = ScalarTraits<Scalar>::zero();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    const bool zeroOut = true;
    MV mvec1(map,numVectors,zeroOut),
       mvec2(map,numVectors,zeroOut);
    Array<Scalar> dots1(numVectors), dots2(numVectors), zeros(numVectors);
    Array<Mag>    norms1(numVectors), norms2(numVectors), ans(numVectors);
    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Scalar>::zero());
    // these should be numerically orthogonal even in finite arithmetic, because both are zero. 1-norms are zero.
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots2,zeros,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
    TEST_EQUALITY_CONST( static_cast<Scalar> (mvec1.getVector(0)->dot(*mvec2.getVector(0))), S0 );
    mvec1.norm1(norms1());
    mvec2.norm1(norms2());
    std::fill(ans.begin(), ans.end(), M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms2,ans,M0);
    // replace local entries s.t.
    // mvec1 = [1 1]  and  mvec2 = [0 0]
    //         [0 0]               [1 1]
    // still numerically orthogonal even in finite arithmetic. norms are numImages.
    for (size_t j=0; j < numVectors; ++j) {
      mvec1.replaceLocalValue(0,j,ScalarTraits<Scalar>::one());
      mvec2.replaceGlobalValue(map->getGlobalElement(1),j,ScalarTraits<Scalar>::one());
    }
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots2,zeros,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
    TEST_EQUALITY_CONST( static_cast<Scalar> (mvec1.getVector(0)->dot(*mvec2.getVector(0))), S0 );
    mvec1.norm1(norms1());
    mvec2.norm1(norms2());
    std::fill(ans.begin(), ans.end(), as<Mag>(numImages));
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms2,ans,M0);
    // sum into local entries s.t.
    // mvec1 = [1 1]  and  mvec2 = [-1 -1]
    //         [1 1]               [ 1  1]
    // still numerically orthogonal even in finite arithmetic. norms are 2*numImages.
    for (size_t j=0; j < numVectors; ++j) {
      mvec1.sumIntoLocalValue(1,j,ScalarTraits<Scalar>::one());
      mvec2.sumIntoGlobalValue(map->getGlobalElement(0),j,-ScalarTraits<Scalar>::one());
    }
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots2,zeros,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
    TEST_EQUALITY_CONST( static_cast<Scalar> (mvec1.getVector(0)->dot(*mvec2.getVector(0))), S0 );
    mvec1.norm1(norms1());
    mvec2.norm1(norms2());
    std::fill(ans.begin(), ans.end(), as<Mag>(2*numImages));
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms2,ans,M0);

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, CopyView, LO , GO , Scalar , Node )
  {
    using std::endl;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

    out << "Test: MultiVector, CopyView" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar S0 = ScalarTraits<Scalar>::zero();
    const Mag M0 = ScalarTraits<Mag>::zero();
    const Mag tol = errorTolSlack * testingTol<Scalar>();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 7;
    const size_t numVectors = 13;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    MV A(map,numVectors,false);
    {
      A.randomize();
      TEUCHOS_TEST_FOR_EXCEPT(numVectors != 13);
      Range1D inds1(8,12);
      // get a subview and a subcopy of certain vectors of A
      // check that the norms are the same
      // change the view, delete it, verify that the copy doesn't change but that A does
      A.randomize();
      Array<Mag> A_bef(numVectors),
                 A_aft (numVectors),
                 Av_bef(inds1.size()),
                 Av_aft (inds1.size()),
                 Ac_bef(inds1.size()),
                 Ac_aft (inds1.size());
      A.norm1(A_bef());
      // get view and its norms
      RCP<MV> Av = A.subViewNonConst(inds1);
      Av->norm1(Av_bef());
      // get copy and its norms
      RCP<MV> Ac = A.subCopy(inds1);
      Ac->norm1(Ac_bef());
      // set view to zero
      Av->putScalar(ScalarTraits<Scalar>::zero());
      // get norms of view
      Av->norm1(Av_aft());
      // free the view, copying data back to A
      Av = Teuchos::null;
      // get norms of A and copy
      Ac->norm1(Ac_aft());
      A.norm1(A_aft());
      // norms of copy and view before should match norms of A
      for (size_t i=0; i < as<size_t>(inds1.size()); ++i) {
        TEST_FLOATING_EQUALITY( A_bef[inds1.lbound()+i], Ac_bef[i], tol );
      }
      TEST_COMPARE_FLOATING_ARRAYS(Ac_bef,Av_bef,tol);
      // norms of copy (before and after) should match
      TEST_COMPARE_FLOATING_ARRAYS(Ac_bef,Ac_aft,tol);
      // norms of view after should be zero, as should corresponding A norms
      for (size_t i=0; i < as<size_t>(inds1.size()); ++i) {
        TEST_EQUALITY_CONST( Av_aft[i], M0 );
        TEST_EQUALITY_CONST( A_aft[inds1.lbound()+i], M0 );
      }
    }
    {
      A.randomize();
      TEUCHOS_TEST_FOR_EXCEPT(numVectors != 13);
      Tuple<size_t,5> inds = tuple<size_t>(0,5,6,7,12);
      // get a subview and a subcopy of certain vectors of A
      // check that the norms are the same
      // change the view, delete it, verify that the copy doesn't change but that A does
      Array<Mag> A_bef(numVectors),
                 A_aft (numVectors),
                 Av_bef(inds.size()),
                 Av_aft (inds.size()),
                 Ac_bef(inds.size()),
                 Ac_aft (inds.size());
      A.norm1(A_bef());
      // get view and its norms
      RCP<MV> Av = A.subViewNonConst(inds);
      Av->norm1(Av_bef());
      // get copy and its norms
      RCP<MV> Ac = A.subCopy(inds);
      Ac->norm1(Ac_bef());
      // set view to zero
      Av->putScalar(ScalarTraits<Scalar>::zero());
      // get norms of view
      Av->norm1(Av_aft());
      // free the view, copying data back to A
      Av = Teuchos::null;
      // get norms of A and copy
      Ac->norm1(Ac_aft());
      A.norm1(A_aft());
      // norms of copy and view before should match norms of A
      for (size_t i=0; i < as<size_t>(inds.size()); ++i) {
        TEST_FLOATING_EQUALITY( A_bef[inds[i]], Ac_bef[i], tol );
      }
      TEST_COMPARE_FLOATING_ARRAYS(Ac_bef,Av_bef,tol);
      // norms of copy (before and after) should match
      TEST_COMPARE_FLOATING_ARRAYS(Ac_bef,Ac_aft,tol);
      // norms of view after should be zero, as should corresponding A norms
      for (size_t i=0; i < as<size_t>(inds.size()); ++i) {
        TEST_EQUALITY_CONST( Av_aft[i], M0 );
        TEST_EQUALITY_CONST( A_aft[inds[i]], M0 );
      }
    }
    {
      A.randomize();
      Array<Mag> Anorms(numVectors);
      A.norm1(Anorms());
      TEUCHOS_TEST_FOR_EXCEPT(numVectors != 13);
      for (size_t vc=0; vc < 2; ++vc) {
        // vc == 0 -> view
        // vc == 1 -> copy
        for (size_t t=0; t < 4; ++t) {
          //  t |   outer   |   inner
          // ---|-----------|-----------
          //  0 | ArrayView | ArrayView
          //  1 |  Range1D  | ArrayView
          //  2 | ArrayView |  Range1D
          //  3 |  Range1D  |  Range1D
          //
          // outer grabs 5-9
          // inner grabs 1-3 of those, corresponding to 6-8
          RCP<const MV> sub1, sub2;
          if ((t & 1) == 0) {
            Tuple<size_t,5> inds = tuple<size_t>(5,6,7,8,9);
            if (vc == 0) sub1 = A.subView(inds);
            else         sub1 = A.subCopy(inds);
          }
          else {
            Range1D inds(5,9);
            if (vc == 0) sub1 = A.subView(inds);
            else         sub1 = A.subCopy(inds);
          }
          TEST_EQUALITY_CONST(sub1->getNumVectors(), 5);
          if ((t & 2) == 0) {
            Tuple<size_t,3> inds = tuple<size_t>(1,2,3);
            if (vc == 0) sub2 = sub1->subView(inds);
            else         sub2 = sub1->subCopy(inds);
          }
          else {
            Range1D inds(1,3);
            if (vc == 0) sub2 = sub1->subView(inds);
            else         sub2 = sub1->subCopy(inds);
          }
          TEST_EQUALITY_CONST(sub2->getNumVectors(), 3);
          Array<Mag> subnorms(3);
          sub2->norm1(subnorms());
          TEST_COMPARE_FLOATING_ARRAYS(Anorms(6,3),subnorms(),tol);
        }
      }
    }
    {
      A.randomize();

      out << "Check that get1dView and get1dCopy have the same values (type 1)" << endl;
      {
        ArrayRCP<const Scalar> view;
        Array<Scalar> copy(numLocal*numVectors);
        TEST_NOTHROW( view = A.get1dView() );
        TEST_NOTHROW( A.get1dCopy(copy,numLocal) );
        TEST_COMPARE_FLOATING_ARRAYS(view,copy,M0);
      }

      A.randomize();
      out << "Check that get1dView and get1dCopy have the same values (type 2)" << endl;
      {
        ArrayRCP<const Scalar> view;
        Array<Scalar> copy(numLocal*numVectors);
        TEST_NOTHROW( A.get1dCopy(copy,numLocal) );
        TEST_NOTHROW( view = A.get1dView() );
        TEST_COMPARE_FLOATING_ARRAYS(view,copy,M0);
      }

      out << "Check that get1dViewNonConst and get1dCopy have the same values" << endl;
      {
        ArrayRCP<Scalar> view;
        Array<Scalar> copy(numLocal*numVectors);
        TEST_NOTHROW( A.get1dCopy(copy,numLocal) );
        TEST_NOTHROW( view = A.get1dViewNonConst() );
        TEST_COMPARE_FLOATING_ARRAYS(view,copy,M0);
        // clear view, ensure that A is zero
        std::fill(view.begin(), view.end(), S0);
        view = Teuchos::null;
        Array<Mag> norms(numVectors), zeros(numVectors,M0);
        A.norm1(norms());
        TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
      }
    }
    {
      A.randomize();

      out << "Check that get2dView and get2dCopy have the same values" << endl;
      {
        ArrayRCP<ArrayRCP<const Scalar> > views;
        Array<Scalar> copyspace(numLocal*numVectors);
        Array<ArrayView<Scalar> > copies(numVectors);
        for (size_t j=0; j < numVectors; ++j) {
          copies[j] = copyspace(numLocal*j,numLocal);
        }
        TEST_NOTHROW( A.get2dCopy(copies) );
        TEST_NOTHROW( views = A.get2dView() );
        for (size_t j=0; j < numVectors; ++j) {
          TEST_COMPARE_FLOATING_ARRAYS(views[j],copies[j],M0);
        }
      }

      out << "Check that get2dViewNonConst and get2dCopy have the same values" << endl;
      {
        ArrayRCP<ArrayRCP<Scalar> > views;
        Array<Scalar> copyspace(numLocal*numVectors);
        Array<ArrayView<Scalar> > copies(numVectors);
        for (size_t j=0; j < numVectors; ++j) {
          copies[j] = copyspace(numLocal*j,numLocal);
        }
        TEST_NOTHROW( A.get2dCopy(copies()) );
        TEST_NOTHROW( views = A.get2dViewNonConst() );
        for (size_t j=0; j < numVectors; ++j) {
          TEST_COMPARE_FLOATING_ARRAYS(views[j],copies[j],M0);
        }
        // clear view, ensure that A is zero
        for (size_t j=0; j < numVectors; ++j) {
          std::fill(views[j].begin(), views[j].end(), S0);
        }
        views = Teuchos::null;
        Array<Mag> norms(numVectors), zeros(numVectors,M0);
        A.norm1(norms());
        TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, OffsetView, LO , GO , Scalar , Node )
  {
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

    out << "Test: MultiVector, OffsetView" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar S0 = ScalarTraits<Scalar>::zero();
    const Mag M0 = ScalarTraits<Mag>::zero();
    const Mag tol = errorTolSlack * testingTol<Scalar>();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal1 = 3;
    const size_t numLocal2 = 4;
    const size_t numLocal = numLocal1 + numLocal2;
    const size_t numVectors = 6;
    Array<size_t> even(tuple<size_t>(1,3,5));
    Array<size_t>  odd(tuple<size_t>(0,2,4));
    TEUCHOS_TEST_FOR_EXCEPTION( even.size() != odd.size(), std::logic_error, "Test setup assumption violated.");
    RCP<const Map<LO,GO,Node> > fullMap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    RCP<const Map<LO,GO,Node> >    map1 = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal1,comm);
    RCP<const Map<LO,GO,Node> >    map2 = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal2,comm);
    RCP<MV> A = rcp(new MV(fullMap,numVectors,false));
    {
      // contig source multivector
      RCP<MV> A1 = A->offsetViewNonConst(map1, 0);
      RCP<MV> A2 = A->offsetViewNonConst(map2, numLocal1);
      TEST_EQUALITY( A1->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A2->getLocalLength(), numLocal2 );
      TEST_EQUALITY( A1->getNumVectors(), numVectors );
      TEST_EQUALITY( A2->getNumVectors(), numVectors );
      Array<Mag>  A_befr(numVectors),
                 A1_befr(numVectors),
                 A2_befr(numVectors),
                  A_aft1(numVectors),
                 A1_aft1(numVectors),
                 A2_aft1(numVectors),
                  A_aft2(numVectors),
                 A1_aft2(numVectors),
                 A2_aft2(numVectors);
      // compute norms of A, A1 and A2
      A->randomize();
      A->norm2(A_befr());
      A1->norm2(A1_befr());
      A2->norm2(A2_befr());
      // set A1 = zeros, compute norms of A, A1 and A2
      A1->putScalar(S0);
      A->norm2(A_aft1());
      A1->norm2(A1_aft1());
      A2->norm2(A2_aft1());
      // set A2 = zeros, compute norms of A, A1 and A2
      A2->putScalar(S0);
      A->norm2(A_aft2());
      A1->norm2(A1_aft2());
      A2->norm2(A2_aft2());
      // change to A1 should not affect A2
      // change to A2 should not affect A1
      // change to A1 or A2 should change A
      // A should be zero after setting A1 to zero and A2 to zero
      for (size_t i=0; i<numVectors; ++i) {
        TEST_EQUALITY_CONST( A_aft1[i] < A_befr[i] + tol, true ); // shrunk as A1 = 0
        TEST_EQUALITY_CONST( A_aft2[i] < A_aft1[i] + tol, true ); // shrunk as A2 = 0
        TEST_EQUALITY_CONST( A_aft2[i] , M0 );                    // ... to zero
        TEST_EQUALITY_CONST( A1_aft1[i] , M0 );                   // was set to zero
        TEST_EQUALITY_CONST( A1_aft2[i] , M0 );                   // should not have been changed
        TEST_FLOATING_EQUALITY( A2_befr[i], A2_aft1[i], tol);     // should not have been changed
        TEST_EQUALITY_CONST( A2_aft2[i] , M0 );                   // was set to zero
      }
    }

    {
      // contig source multivector: repeat with "offset view"
      // constructor that takes RCP<const Map>.
      MV A1 (*A, map1, 0);
      MV A2 (*A, map2, numLocal1);
      TEST_EQUALITY( A1.getLocalLength(), numLocal1 );
      TEST_EQUALITY( A2.getLocalLength(), numLocal2 );
      TEST_EQUALITY( A1.getNumVectors(), numVectors );
      TEST_EQUALITY( A2.getNumVectors(), numVectors );
      Array<Mag>  A_befr(numVectors),
                 A1_befr(numVectors),
                 A2_befr(numVectors),
                  A_aft1(numVectors),
                 A1_aft1(numVectors),
                 A2_aft1(numVectors),
                  A_aft2(numVectors),
                 A1_aft2(numVectors),
                 A2_aft2(numVectors);
      // compute norms of A, A1 and A2
      A->randomize();
      A->norm2(A_befr());
      A1.norm2(A1_befr());
      A2.norm2(A2_befr());
      // set A1 = zeros, compute norms of A, A1 and A2
      A1.putScalar(S0);
      A->norm2(A_aft1());
      A1.norm2(A1_aft1());
      A2.norm2(A2_aft1());
      // set A2 = zeros, compute norms of A, A1 and A2
      A2.putScalar(S0);
      A->norm2(A_aft2());
      A1.norm2(A1_aft2());
      A2.norm2(A2_aft2());
      // change to A1 should not affect A2
      // change to A2 should not affect A1
      // change to A1 or A2 should change A
      // A should be zero after setting A1 to zero and A2 to zero
      for (size_t i=0; i<numVectors; ++i) {
        TEST_EQUALITY_CONST( A_aft1[i] < A_befr[i] + tol, true ); // shrunk as A1 = 0
        TEST_EQUALITY_CONST( A_aft2[i] < A_aft1[i] + tol, true ); // shrunk as A2 = 0
        TEST_EQUALITY_CONST( A_aft2[i] , M0 );                    // ... to zero
        TEST_EQUALITY_CONST( A1_aft1[i] , M0 );                   // was set to zero
        TEST_EQUALITY_CONST( A1_aft2[i] , M0 );                   // should not have been changed
        TEST_FLOATING_EQUALITY( A2_befr[i], A2_aft1[i], tol);     // should not have been changed
        TEST_EQUALITY_CONST( A2_aft2[i] , M0 );                   // was set to zero
      }
    }

    {
      // contig source multivector: repeat with "offset view"
      // constructor that takes const Map.
      MV A1 (*A, *map1, 0);
      MV A2 (*A, *map2, numLocal1);
      TEST_EQUALITY( A1.getLocalLength(), numLocal1 );
      TEST_EQUALITY( A2.getLocalLength(), numLocal2 );
      TEST_EQUALITY( A1.getNumVectors(), numVectors );
      TEST_EQUALITY( A2.getNumVectors(), numVectors );
      Array<Mag>  A_befr(numVectors),
                 A1_befr(numVectors),
                 A2_befr(numVectors),
                  A_aft1(numVectors),
                 A1_aft1(numVectors),
                 A2_aft1(numVectors),
                  A_aft2(numVectors),
                 A1_aft2(numVectors),
                 A2_aft2(numVectors);
      // compute norms of A, A1 and A2
      A->randomize();
      A->norm2(A_befr());
      A1.norm2(A1_befr());
      A2.norm2(A2_befr());
      // set A1 = zeros, compute norms of A, A1 and A2
      A1.putScalar(S0);
      A->norm2(A_aft1());
      A1.norm2(A1_aft1());
      A2.norm2(A2_aft1());
      // set A2 = zeros, compute norms of A, A1 and A2
      A2.putScalar(S0);
      A->norm2(A_aft2());
      A1.norm2(A1_aft2());
      A2.norm2(A2_aft2());
      // change to A1 should not affect A2
      // change to A2 should not affect A1
      // change to A1 or A2 should change A
      // A should be zero after setting A1 to zero and A2 to zero
      for (size_t i=0; i<numVectors; ++i) {
        TEST_EQUALITY_CONST( A_aft1[i] < A_befr[i] + tol, true ); // shrunk as A1 = 0
        TEST_EQUALITY_CONST( A_aft2[i] < A_aft1[i] + tol, true ); // shrunk as A2 = 0
        TEST_EQUALITY_CONST( A_aft2[i] , M0 );                    // ... to zero
        TEST_EQUALITY_CONST( A1_aft1[i] , M0 );                   // was set to zero
        TEST_EQUALITY_CONST( A1_aft2[i] , M0 );                   // should not have been changed
        TEST_FLOATING_EQUALITY( A2_befr[i], A2_aft1[i], tol);     // should not have been changed
        TEST_EQUALITY_CONST( A2_aft2[i] , M0 );                   // was set to zero
      }
    }

    {
      // non-contig source multivector
      RCP<MV> A1e = A->subViewNonConst(even)->offsetViewNonConst(map1, 0);
      RCP<MV> A2e = A->subViewNonConst(even)->offsetViewNonConst(map2, numLocal1);
      RCP<MV> A1o = A->subViewNonConst(odd)->offsetViewNonConst(map1, 0);
      RCP<MV> A2o = A->subViewNonConst(odd)->offsetViewNonConst(map2, numLocal1);
      TEST_EQUALITY( A1e->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A1o->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A2e->getLocalLength(), numLocal2 );
      TEST_EQUALITY( A2o->getLocalLength(), numLocal2 );
      const size_t numSubVecs = (size_t)even.size();
      TEST_EQUALITY( A1e->getNumVectors(), numSubVecs );
      TEST_EQUALITY( A2e->getNumVectors(), numSubVecs );
      TEST_EQUALITY( A1o->getNumVectors(), numSubVecs );
      TEST_EQUALITY( A2o->getNumVectors(), numSubVecs );
      A->randomize();
      Array<Mag> b1(numSubVecs), b2(numSubVecs), b3(numSubVecs), bw(numVectors); // before putScalar(): unchanged 1, 2, 3; whole
      Array<Mag> a1(numSubVecs), a2(numSubVecs), a3(numSubVecs), aw(numVectors); // after putScalar(): ...
      Array<Mag> changed(numSubVecs), zeros(numSubVecs,M0);
      for (int i=0; i<4; ++i) {
        ArrayView<RCP<MV> > allMVs; // (changed,three unchanged)
        switch (i) {
        case 0:
          allMVs = tuple<RCP<MV> >(A1e,A2e,A1o,A2o); break;
        case 1:
          allMVs = tuple<RCP<MV> >(A2e,A1o,A2o,A1e); break;
        case 2:
          allMVs = tuple<RCP<MV> >(A1o,A2o,A1e,A2e); break;
        case 3:
          allMVs = tuple<RCP<MV> >(A2o,A1e,A2e,A1o); break;
        }
        allMVs[1]->norm2(b1()); allMVs[2]->norm2(b2()); allMVs[3]->norm2(b3());
        A->norm2(bw());
        allMVs[0]->putScalar(S0);
        allMVs[0]->norm2(changed());
        allMVs[1]->norm2(a1()); allMVs[2]->norm2(a2()); allMVs[3]->norm2(a3());
        A->norm2(aw());
        TEST_COMPARE_FLOATING_ARRAYS(b1,a1,tol);
        TEST_COMPARE_FLOATING_ARRAYS(b2,a2,tol);
        TEST_COMPARE_FLOATING_ARRAYS(b3,a3,tol);
        TEST_COMPARE_ARRAYS(changed(), zeros());
        for (size_t ii = 0; ii < numVectors; ++ii) {
          TEST_EQUALITY_CONST( aw[ii] < bw[ii] + tol, true ); // shrunk
        }
      }
    }

    {
      // non-contig source multivector: repeat with "offset view"
      // constructor that takes RCP<const Map>.
      RCP<MV> A1e (new MV (* (A->subViewNonConst (even)), map1, 0));
      RCP<MV> A2e (new MV (* (A->subViewNonConst (even)), map2, numLocal1));
      RCP<MV> A1o (new MV (* (A->subViewNonConst (odd)), map1, 0));
      RCP<MV> A2o (new MV (* (A->subViewNonConst (odd)), map2, numLocal1));

      TEST_EQUALITY( A1e->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A1o->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A2e->getLocalLength(), numLocal2 );
      TEST_EQUALITY( A2o->getLocalLength(), numLocal2 );
      const size_t numSubVecs = (size_t)even.size();
      TEST_EQUALITY( A1e->getNumVectors(), numSubVecs );
      TEST_EQUALITY( A2e->getNumVectors(), numSubVecs );
      TEST_EQUALITY( A1o->getNumVectors(), numSubVecs );
      TEST_EQUALITY( A2o->getNumVectors(), numSubVecs );
      A->randomize();
      Array<Mag> b1(numSubVecs), b2(numSubVecs), b3(numSubVecs), bw(numVectors); // before putScalar(): unchanged 1, 2, 3; whole
      Array<Mag> a1(numSubVecs), a2(numSubVecs), a3(numSubVecs), aw(numVectors); // after putScalar(): ...
      Array<Mag> changed(numSubVecs), zeros(numSubVecs,M0);
      for (int i=0; i<4; ++i) {
        std::vector<RCP<MV> > allMVs; // (changed,three unchanged)
        switch (i) {
        case 0:
          allMVs = {A1e, A2e, A1o, A2o}; break;
        case 1:
          allMVs = {A2e, A1o, A2o, A1e}; break;
        case 2:
          allMVs = {A1o, A2o, A1e, A2e}; break;
        case 3:
          allMVs = {A2o, A1e, A2e, A1o}; break;
        }
        allMVs[1]->norm2(b1()); allMVs[2]->norm2(b2()); allMVs[3]->norm2(b3());
        A->norm2(bw());
        allMVs[0]->putScalar(S0);
        allMVs[0]->norm2(changed());
        allMVs[1]->norm2(a1()); allMVs[2]->norm2(a2()); allMVs[3]->norm2(a3());
        A->norm2(aw());
        TEST_COMPARE_FLOATING_ARRAYS(b1,a1,tol);
        TEST_COMPARE_FLOATING_ARRAYS(b2,a2,tol);
        TEST_COMPARE_FLOATING_ARRAYS(b3,a3,tol);
        TEST_COMPARE_ARRAYS(changed(), zeros());
        for (size_t ii = 0; ii < numVectors; ++ii) {
          TEST_EQUALITY_CONST( aw[ii] < bw[ii] + tol, true ); // shrunk
        }
      }
    }

    {
      RCP<const MV> A1 = A->offsetView(map1, 0);
      RCP<const MV> A2 = A->offsetView(map2, numLocal1);
      TEST_EQUALITY( A1->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A2->getLocalLength(), numLocal2 );
      TEST_EQUALITY( A1->getNumVectors(), numVectors );
      TEST_EQUALITY( A2->getNumVectors(), numVectors );
      Array<Mag>  A_bef(numVectors),
                 A1_bef(numVectors),
                 A2_bef(numVectors),
                  A_aft(numVectors),
                 A1_aft(numVectors),
                 A2_aft(numVectors);
      // compute norms of A, A1 and A2
      A->randomize();
      A->norm2(A_bef());
      A1->norm2(A1_bef());
      A2->norm2(A2_bef());
      A->putScalar(S0);
      A->norm2(A_aft());
      A1->norm2(A1_aft());
      A2->norm2(A2_aft());
      for (size_t i=0; i<numVectors; ++i) {
        TEST_EQUALITY_CONST( A_bef[i] < A1_bef[i] + A2_bef[i] + tol, true );
        TEST_EQUALITY_CONST( A_aft[i], S0 );
        TEST_EQUALITY_CONST( A1_aft[i], S0 );
        TEST_EQUALITY_CONST( A2_aft[i], S0 );
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  //This test is for issue #9160. Row subview of a col subview (e.g. A->getVector(i)->getOffsetView(...))
  //would have numVectors() equal to A->numVectors(), not the col subview's numVectors() as it should.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, VectorOffsetView, LO , GO , Scalar , Node )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

    out << "Test: MultiVector, VectorOffsetView" << endl;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal1 = 3;
    const size_t numLocal2 = 4;
    const size_t numLocal = numLocal1 + numLocal2;
    const size_t numVectors = 3;
    RCP<const Map<LO,GO,Node> > fullMap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    RCP<const Map<LO,GO,Node> >    map1 = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal1,comm);
    RCP<const Map<LO,GO,Node> >    map2 = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal2,comm);
    RCP<MV> A = rcp(new MV(fullMap,numVectors,false));
    RCP<MV> Acol0 = A->getVectorNonConst(0);
    //MultiVector has 3 interfaces for the same thing:
    // 1.) MV(MV&, RCP<Map>, size_t): offset view of 1st arg, using given map.
    //     Other 2 interfaces implemented in terms of this.
    // 2.) MV(MV&, Map&, LO): same as 1.) but with different arg types, and it copy-ctors the map (shallow copy)
    // 3.) Methods offsetView/offsetViewNonConst (RCP<Map>, size_t)

    //Test first with a zero offset (corresponding to map1)
    //1
    {
      RCP<MV> Acol0sub = rcp(new MV(*Acol0, map1, 0));
      TEST_EQUALITY( Acol0sub->getNumVectors(), Acol0->getNumVectors() );
    }
    //2
    {
      RCP<MV> Acol0sub = rcp(new MV(*Acol0, *map1, 0));
      TEST_EQUALITY( Acol0sub->getNumVectors(), Acol0->getNumVectors() );
    }
    //3
    {
      {
        RCP<MV> Acol0sub = Acol0->offsetViewNonConst(map1, 0);
        TEST_EQUALITY( Acol0sub->getNumVectors(), Acol0->getNumVectors() );
      }
      {
        RCP<const MV> Acol0sub = Acol0->offsetView(map1, 0);
        TEST_EQUALITY( Acol0sub->getNumVectors(), Acol0->getNumVectors() );
      }
    }
    //Now test with offset numLocal1, corresponding to map2's
    //1
    {
      RCP<MV> Acol0sub = rcp(new MV(*Acol0, map2, numLocal1));
      TEST_EQUALITY( Acol0sub->getNumVectors(), Acol0->getNumVectors() );
    }
    //2
    {
      RCP<MV> Acol0sub = rcp(new MV(*Acol0, *map2, numLocal1));
      TEST_EQUALITY( Acol0sub->getNumVectors(), Acol0->getNumVectors() );
    }
    //3
    {
      {
        RCP<MV> Acol0sub = Acol0->offsetViewNonConst(map2, numLocal1);
        TEST_EQUALITY( Acol0sub->getNumVectors(), Acol0->getNumVectors() );
      }
      {
        RCP<const MV> Acol0sub = Acol0->offsetView(map2, numLocal1);
        TEST_EQUALITY( Acol0sub->getNumVectors(), Acol0->getNumVectors() );
      }
    }
  }


  // This unit test exercises the following situation: Given a
  // Tpetra::MultiVector X, partition it into row blocks [X1; X2]
  // (Matlab notation) using offsetView (or offsetViewNonConst).  The
  // sum of the local number of rows in X1 and X2 equals the local
  // number of rows in X, but either X1 or X2 might have a zero number
  // of local rows.  We exercise each of the latter cases, in two
  // separate tests.  Repeat both cases for offsetView (const X1 and
  // X2) and offsetViewNonConst (nonconst X1 and X2).
  //
  // The most interesting thing this test exercises is that neither of
  // the above cases should throw exceptions.  This was not originally
  // true for the case where X2 has zero local rows.  Thanks to
  // Deaglan Halligan for pointing this out (on 23 Oct 2013).
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, OffsetViewZeroLength, LO , GO , Scalar , Node )
  {
    typedef Tpetra::global_size_t GST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    out << "Test: MultiVector, OffsetViewZeroLength" << endl;
    Teuchos::OSTab tab0 (out);

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    // Get a communicator and Kokkos node instance.
    RCP<const Comm<int> > comm = getDefaultComm ();

    // Create a Map with a nonzero number of entries on each process.
    const size_t numLocalEntries = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocalEntries, indexBase, comm));

    // Create a MultiVector X using that Map.  Give it some number of
    // columns (vectors) other than 1, just to exercise the most
    // general case.
    const size_t numVecs = 3;
    MV X (map, numVecs);

    // Make sure that X has the right (local) dimensions.
    TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
    TEST_EQUALITY( X->getNumVectors (), numVecs );

    // Create a Map with zero entries on every process.
    RCP<const map_type> mapZero =
      rcp (new map_type (INVALID, 0, indexBase, comm));

    // Case 1: X1 has the same local number of rows as X, and X2 has
    // zero local rows.  Thus, X2 will be a zero-length view of X,
    // starting at the end of the local part of X (so the offset is
    // numLocalEntries).
    {
      RCP<const MV> X1;
      RCP<const MV> X2;
      try {
        X1 = X.offsetView (map, 0);
        X2 = X.offsetView (mapZero, numLocalEntries);
      } catch (...) {
        out << "The following case failed: X = [X1; X2] where X2 has zero "
          "local rows." << std::endl;
        throw;
      }
      // Make sure that offsetView() didn't change X's dimensions.
      TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X->getNumVectors (), numVecs );

      // Make sure that X1 and X2 have the right (local) dimensions.
      TEST_EQUALITY( X1->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X1->getNumVectors (), numVecs );
      TEST_EQUALITY_CONST( X2->getLocalLength (), static_cast<size_t> (0) );
      TEST_EQUALITY( X2->getNumVectors (), numVecs );

      // Make sure the pointers are the same, by extracting the
      // Kokkos::DualView objects.  Get the host pointer, just in case
      // MV allocation favors host space for initial allocations and
      // defers device allocations.

      auto X_local = X->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X1_local = X1->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X2_local = X2->getLocalViewHost(Tpetra::Access::ReadOnly);

      // Make sure the pointers match.  It doesn't really matter to
      // what X2_local points, as long as it has zero rows.
      TEST_EQUALITY( X1_local.data (), X_local.data () );

      // Make sure the local dimensions of X1 are correct.
      TEST_EQUALITY( X1_local.extent (0), X_local.extent (0) );
      TEST_EQUALITY( X1_local.extent (1), X_local.extent (1) );

      // Make sure the local dimensions of X2 are correct.
      TEST_EQUALITY_CONST( X2_local.extent (0), static_cast<size_t> (0) );
      TEST_EQUALITY( X2_local.extent (1), X_local.extent (1) );

      // Make sure that nothing bad happens on deallocation.
      try {
        X1 = Teuchos::null;
        X2 = Teuchos::null;
      } catch (...) {
        out << "Failed to deallocate X1 or X2." << std::endl;
        throw;
      }
    }

    // Nonconst version of Case 1.
    {
      RCP<MV> X1_nonconst;
      RCP<MV> X2_nonconst;
      try {
        X1_nonconst = X.offsetViewNonConst (map, 0);
        X2_nonconst = X.offsetViewNonConst (mapZero, numLocalEntries);
      } catch (...) {
        out << "The following case failed: X = [X1; X2] where X2 has zero "
          "local rows, and X1 and X2 are nonconst." << std::endl;
        throw;
      }
      // Make sure that offsetView() didn't change X's dimensions.
      TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X->getNumVectors (), numVecs );

      // Make sure that X1 and X2 have the right (local) dimensions.
      TEST_EQUALITY( X1_nonconst->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X1_nonconst->getNumVectors (), numVecs );
      TEST_EQUALITY_CONST( X2_nonconst->getLocalLength (), static_cast<size_t> (0) );
      TEST_EQUALITY( X2_nonconst->getNumVectors (), numVecs );

      // Make sure the pointers are the same, by extracting the
      // Kokkos::DualView objects.  Get the host pointer, just in case
      // MV allocation favors host space for initial allocations and
      // defers device allocations.

      auto X_local = X->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X1_local = X1_nonconst->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X2_local = X2_nonconst->getLocalViewHost(Tpetra::Access::ReadOnly);

      // Make sure the pointers match.  It doesn't really matter to
      // what X2_local points, as long as it has zero rows.
      TEST_EQUALITY( X1_local.data (), X_local.data () );

      // Make sure the local dimensions of X1 are correct.
      TEST_EQUALITY( X1_local.extent (0), X_local.extent (0) );
      TEST_EQUALITY( X1_local.extent (1), X_local.extent (1) );

      // Make sure the local dimensions of X2 are correct.
      TEST_EQUALITY_CONST( X2_local.extent (0), static_cast<size_t> (0) );
      TEST_EQUALITY( X2_local.extent (1), X_local.extent (1) );

      // Make sure that nothing bad happens on deallocation.
      try {
        X1_nonconst = Teuchos::null;
        X2_nonconst = Teuchos::null;
      } catch (...) {
        out << "Failed to deallocate X1 or X2." << std::endl;
        throw;
      }
    }

    // Case 2: X1 has zero rows, and X2 has the same local number of
    // rows as X.  Thus, X1 will be a zero-length view of X, starting
    // at the beginning of the local part of X (so the offset is 0).
    {
      RCP<const MV> X1;
      RCP<const MV> X2;
      try {
        X1 = X.offsetView (mapZero, 0);
        X2 = X.offsetView (map, 0);
      } catch (...) {
        out << "The following case failed: X = [X1; X2] where X1 has zero "
          "local rows." << std::endl;
        throw;
      }
      // Make sure that offsetView() didn't change X's dimensions.
      TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X->getNumVectors (), numVecs );

      // Make sure that X1 and X2 have the right (local) dimensions.
      TEST_EQUALITY_CONST( X1->getLocalLength (), static_cast<size_t> (0) );
      TEST_EQUALITY( X1->getNumVectors (), numVecs );
      TEST_EQUALITY( X2->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X2->getNumVectors (), numVecs );

      // Make sure the pointers are the same, by extracting the
      // Kokkos::DualView objects.  Get the host pointer, just in case
      // MV allocation favors host space for initial allocations and
      // defers device allocations.

      auto X_local = X->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X1_local = X1->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X2_local = X2->getLocalViewHost(Tpetra::Access::ReadOnly);
      // Make sure the pointers match.  It doesn't really matter to
      // what X1_local points, as long as it has zero rows.
      TEST_EQUALITY( X2_local.data (), X_local.data () );

      // Make sure the local dimensions of X1 are correct.
      TEST_EQUALITY_CONST( X1_local.extent (0), static_cast<size_t> (0) );
      TEST_EQUALITY( X1_local.extent (1), X_local.extent (1) );

      // Make sure the local dimensions of X2 are correct.
      TEST_EQUALITY( X2_local.extent (0), X_local.extent (0) );
      TEST_EQUALITY( X2_local.extent (1), X_local.extent (1) );

      // Make sure that nothing bad happens on deallocation.
      try {
        X1 = Teuchos::null;
        X2 = Teuchos::null;
      } catch (...) {
        out << "Failed to deallocate X1 or X2." << std::endl;
        throw;
      }
    }

    // Nonconst version of Case 2.
    {
      RCP<MV> X1_nonconst;
      RCP<MV> X2_nonconst;
      try {
        X1_nonconst = X.offsetViewNonConst (mapZero, 0);
        X2_nonconst = X.offsetViewNonConst (map, 0);
      } catch (...) {
        out << "The following case failed: X = [X1; X2] where X1 has zero "
          "local rows, and X1 and X2 are nonconst." << std::endl;
        throw;
      }
      // Make sure that offsetView() didn't change X's dimensions.
      TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X->getNumVectors (), numVecs );

      // Make sure that X1 and X2 have the right (local) dimensions.
      TEST_EQUALITY_CONST( X1_nonconst->getLocalLength (), static_cast<size_t> (0) );
      TEST_EQUALITY( X1_nonconst->getNumVectors (), numVecs );
      TEST_EQUALITY( X2_nonconst->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X2_nonconst->getNumVectors (), numVecs );

      // Make sure the pointers are the same, by extracting the
      // Kokkos::DualView objects.  Get the host pointer, just in case
      // MV allocation favors host space for initial allocations and
      // defers device allocations.

      auto X_local = X->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X1_local = X1_nonconst->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X2_local = X2_nonconst->getLocalViewHost(Tpetra::Access::ReadOnly);

      // Make sure the pointers match.  It doesn't really matter to
      // what X1_local points, as long as it has zero rows.
      TEST_EQUALITY( X2_local.data (), X_local.data () );

      // Make sure the local dimensions of X1 are correct.
      TEST_EQUALITY_CONST( X1_local.extent (0), static_cast<size_t> (0) );
      TEST_EQUALITY( X1_local.extent (1), X_local.extent (1) );

      // Make sure the local dimensions of X2 are correct.
      TEST_EQUALITY( X2_local.extent (0), X_local.extent (0) );
      TEST_EQUALITY( X2_local.extent (1), X_local.extent (1) );

      // Make sure that nothing bad happens on deallocation.
      try {
        X1_nonconst = Teuchos::null;
        X2_nonconst = Teuchos::null;
      } catch (...) {
        out << "Failed to deallocate X1 or X2." << std::endl;
        throw;
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, ZeroScaleUpdate, LO , GO , Scalar , Node )
  {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::global_size_t GST;

    out << "Test: MultiVector, ZeroScaleUpdate" << endl;
    Teuchos::OSTab tab0 (out);

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const Mag M0 = ScalarTraits<Mag>::zero ();

    const Scalar zero = STS::zero ();
    const Scalar one = STS::one ();
    const Scalar two = one + one;
    const Scalar four = two + two;

    RCP<const Comm<int> > comm = getDefaultComm();

    // create a Map
    const size_t numLocal = 2;
    const size_t numVectors = 2;
    const size_t LDA = 2;
    RCP<const Map<LO, GO, Node> > map =
      createContigMapWithNode<LO, GO, Node> (INVALID, numLocal, comm);
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
    values[0] = one;
    values[1] = one;
    values[2] = two;
    values[3] = two;
    values[4] = four;
    values[5] = four;
    MV A (map, values (0,4), LDA, numVectors);
    MV B (map, values (2,4), LDA, numVectors);
    Array<Mag> norms (numVectors);
    Array<Mag> zeros (numVectors);
    std::fill (zeros.begin (), zeros.end (), M0);
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
      MV A2 (A, Teuchos::Copy);
      A2.scale (two);
      A2.update (-one, B, one);
      A2.norm1 (norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }
    //   set A2 = A
    //   check that it equals B: scale,subtraction in situ
    {
      MV A2 (A, Teuchos::Copy);
      A2.update (-one, B, two);
      A2.norm1 (norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }
    //   set C random
    //   set it to zero by combination with A,B
    {
      MV C (map, numVectors);
      C.randomize ();
      C.update (-one, B, two, A, zero);
      C.norm1 (norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }
    //   set C random
    //   scale it ex-situ
    //   check that it equals B: subtraction in situ
    {
      MV C (map, numVectors);
      C.scale (two, A);
      C.update (one, B, -one);
      C.norm1 (norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, ScaleAndAssign, LO , GO , Scalar , Node )
  {
    std::cerr << std::endl;

    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::Vector<Scalar,LO,GO,Node>       V;

    out << "Test: MultiVector, ScaleAndAssign" << endl;
    Teuchos::OSTab tab0 (out);

    int lclSuccess = 1;
    int gblSuccess = 0; // to be set below

    Teuchos::ScalarTraits<Scalar>::seedrandom(0);   // consistent seed
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Mag tol = errorTolSlack * testingTol<Scalar>();
    const Mag M0 = ScalarTraits<Mag>::zero();

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 23;
    const size_t numVectors = 11;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    // Use random multivector A
    // Set B = A * 2 manually.
    // Therefore, if C = 2*A, then C == B
    // If C = A and C *= 2, then C == B
    // This test operator= and all of our scale ops
    // We'll do Vector and MultiVector variations
    // Also, ensure that other vectors aren't changed

    out << "Create A, and fill with random numbers" << endl;
    MV A(map, numVectors, false);
    A.randomize();

    out << "Stash away norms of columns of A" << endl;
    Array<Mag> Anrms(numVectors);
    A.norm1(Anrms());

    out << "Test B := A*2, using different methods" << endl;
    // set B = A * 2, using different techniques
    // * deep_copy(B,A) and scale B in place
    // * get 1-vector subview(Range1D), MultiVector::operator=
    // * get 1-vector subview(ArrayView), MultiVector::operator=
    // * get data view, assign
    TEUCHOS_TEST_FOR_EXCEPT(numVectors < 4);

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (! gblSuccess) {
      return; // no point in continuing
    }

    MV B (map, numVectors, false);
    for (size_t j = 0; j < numVectors; ++j) {
      Teuchos::OSTab tab1 (out);

      // assign j-th vector of B to 2 * j-th vector of A
      switch (j % 4) {
        case 0:
          {
            std::ostringstream os;
            os << ">>> Proc " << comm->getSize ();
            os << ": A.modified_host: " << (A.need_sync_device ()?1:0);
            os << ", A.modified_device: " << (A.need_sync_host ()?1:0);
            os << ", B.modified_host: " << (B.need_sync_device ()?1:0);
            os << ", B.modified_device: " << (B.need_sync_host ()?1:0);
            os << std::endl;
            std::cerr << os.str ();
          }
          {
            out << "Method 0" << endl;

            RCP<V> bj = B.getVectorNonConst(j);
            RCP<const V> aj = A.getVector(j);
            deep_copy((*bj),(*aj));

            ArrayRCP<Scalar> bjview = bj->get1dViewNonConst();
            for (size_t i=0; i < numLocal; ++i) {
              bjview[i] *= as<Scalar>(2);
            }
          }
          
          lclSuccess = success ? 1 : 0;
          gblSuccess = 0; // output argument
          reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
          TEST_EQUALITY_CONST( gblSuccess, 1 );
          if (! gblSuccess) {
            return; // no point in continuing
          }
          break;
        case 1:
          {
            out << "Method 1" << endl;

            RCP<MV>       bj = B.subViewNonConst(Range1D(j,j));
            RCP<const MV> aj = A.subView(Range1D(j,j));
            ////(*bj) = (*aj);
            deep_copy((*bj),(*aj));
            ArrayRCP<Scalar> bjview = bj->get1dViewNonConst();
            for (size_t i=0; i < numLocal; ++i) {
              bjview[i] *= as<Scalar>(2);
            }
          }
          lclSuccess = success ? 1 : 0;
          gblSuccess = 0; // output argument
          reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
          TEST_EQUALITY_CONST( gblSuccess, 1 );
          if (! gblSuccess) {
            return; // no point in continuing
          }
          break;
        case 2:
          {
            out << "Method 2" << endl;

            RCP<MV> bj = B.subViewNonConst(tuple<size_t>(j));
            RCP<const MV> aj = A.subView(tuple<size_t>(j));
            //RCP<MV>       bj = B.subViewNonConst(Range1D(j,j));
            //RCP<const MV> aj = A.subView(Range1D(j,j));
            //(*bj) = (*aj);
            deep_copy((*bj),(*aj));
            ArrayRCP<const Scalar> ajview = aj->get1dView();
            ArrayRCP<Scalar> bjview = bj->get1dViewNonConst();
            for (size_t i=0; i < numLocal; ++i) {
              bjview[i] *= as<Scalar>(2);
            }
          }
          lclSuccess = success ? 1 : 0;
          gblSuccess = 0; // output argument
          reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
          TEST_EQUALITY_CONST( gblSuccess, 1 );
          if (! gblSuccess) {
            return; // no point in continuing
          }
          break;
        case 3:
          {
            out << "Method 3" << endl;

            ArrayRCP<Scalar>       bjview = B.getDataNonConst(j);
            ArrayRCP<const Scalar> ajview = A.getData(j);
            for (size_t i=0; i < numLocal; ++i) {
              bjview[i] = as<Scalar>(2) * ajview[i];
            }
          }
          lclSuccess = success ? 1 : 0;
          gblSuccess = 0; // output argument
          reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
          TEST_EQUALITY_CONST( gblSuccess, 1 );
          if (! gblSuccess) {
            return; // no point in continuing
          }
          break;
      }
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (! gblSuccess) {
      return; // no point in continuing
    }

    out << "Check that A wasn't modified" << endl;
    // check that A wasn't modified
    {
      Array<Mag> Anrms_aft(numVectors);
      A.norm1(Anrms_aft());
      TEST_COMPARE_FLOATING_ARRAYS(Anrms(),Anrms_aft(),tol);
    }

    out << "Check that C.scale(2,A) == B" << endl;
    // check that C.Scale(A,2.0) == B
    {
      MV C (map, numVectors, false);
      C.scale (as<Scalar> (2), A);
      C.update (-1.0,B,1.0);

      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm1(Cnorms());
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms(),zeros,tol);
    }

    out << "Check that C := A, C.scale(2) == B" << endl;
    // check that C=A, C.Scale(2.0) == B
    {
      MV C (createCopy(A));
      C.scale(as<Scalar>(2));
      C.update(-1.0,B,1.0);
      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm1(Cnorms);
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms,zeros,tol);
    }

    out << "Check that C := A, C.scale(tuple(2)) == B" << endl;
    // check that C=A, C.Scale(tuple(2)) == B
    {
      MV C(createCopy(A));
      Array<Scalar> twos(numVectors,as<Scalar>(2));
      C.scale(twos());
      C.update(-1.0,B,1.0);
      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm1(Cnorms());
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms(),zeros,tol);
    }

    // Make sure that we succeeded on all processes.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Vector, ZeroScaleUpdate, LO , GO , Scalar , Node )
  {
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::Vector<Scalar,LO,GO,Node>       V;

    out << "Test: Vector, ZeroScaleUpdate" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,2,comm);
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
      V A2(createCopy(A));
      A2.scale(as<Scalar>(2));
      A2.update(as<Scalar>(-1),B,as<Scalar>(1));
      norm = A2.norm1(); A2.norm1(norms());
      TEST_EQUALITY(norm,M0);
      TEST_EQUALITY(norm,norms[0]);
    }
    //   set A2 = A
    //   check that it equals B: scale,subtraction in situ
    {
      V A2(createCopy(A));
      A2.update(as<Scalar>(-1),B,as<Scalar>(2));
      norm = A2.norm1(); A2.norm1(norms());
      TEST_EQUALITY(norm,M0);
      TEST_EQUALITY(norm,norms[0]);
    }
    //   set C random
    //   set it to zero by combination with A,B
    {
      V C(map);
      C.randomize();
      C.update(as<Scalar>(-1),B,as<Scalar>(2),A,as<Scalar>(0));
      norm = C.norm1(); C.norm1(norms());
      TEST_EQUALITY(norm,M0);
      TEST_EQUALITY(norm,norms[0]);
    }
    //   set C random
    //   scale it ex-situ
    //   check that it equals B: subtraction in situ
    {
      V C(map);
      C.randomize();
      C.scale(as<Scalar>(2),A);
      C.update(as<Scalar>(1),B,as<Scalar>(-1));
      norm = C.norm1(); C.norm1(norms());
      TEST_EQUALITY(norm,M0);
      TEST_EQUALITY(norm,norms[0]);
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, CopyConst, LO , GO , Scalar , Node )
  {
    using std::endl;
    using Teuchos::toString;
    typedef Tpetra::global_size_t GST;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename MV::mag_type Mag;
    typedef Teuchos::ScalarTraits<Mag> STM;

    out << "Test: MultiVector, CopyConst" << endl;
    Teuchos::OSTab tab0 (out);

    const Mag M0 = STM::zero ();
    // This test should even pass in the field of the integers mod 2.
    // In that case, TWO == ZERO, THREE == ONE, and FOUR == TWO,
    // though, so the test won't be very conclusive.
    const Scalar ZERO = STS::zero ();
    const Scalar ONE = STS::one ();
    const Scalar TWO = ONE + ONE;
    const Scalar THREE = TWO + ONE;
    const Scalar FOUR = THREE + ONE;
    const Scalar FIVE = FOUR + ONE;

    // Create a Map
    const size_t numLocal = 13;
    const size_t numVectors = 7;
    const GO indexBase = 0;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    auto comm = getDefaultComm ();
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));

    out << "Part 1:" << endl;
    {
      Teuchos::OSTab tab1 (out);

      // Create original MultiVector.  Fill it with nonzero initial
      // data.  Thus, if we see zeros in the inf-norms below, we know
      // that something is wrong.
      out << "Create original MultiVector (mvorig)" << endl;
      MV mvorig (map, numVectors);
      out << "Fill all entries of mvorig with " << FIVE << endl;
      mvorig.putScalar (FIVE);

      // Remember the inf-norms of all of the columns of mvorig.
      Array<Mag> norig (numVectors);
      mvorig.normInf (norig ());
      out << "Inf-norms of mvorig's columns: " << toString (norig ()) << endl;

      // Create nonconst strided (noncontiguous) subview of mvorig.
      Tuple<size_t,3> inds = tuple<size_t> (1,3,5);
      out << "mvview = mvorig.subViewNonConst(" << toString (inds) << ")" << endl;
      RCP<MV> mvview = mvorig.subViewNonConst (inds);

      Array<Mag> nsub (inds.size ());
      for (size_t j = 0; j < static_cast<size_t> (inds.size ()); ++j) {
        nsub[j] = norig[inds[j]];
      }

      out << "Make three different deep copies of mvview" << endl;

      // Make a deep copy of mvview, in three different ways.  All
      // three copies should be independent of mvview and of each
      // other.

      out << "  MV mvcopy (*mvview, Teuchos::Copy);" << endl;
      MV mvcopy (*mvview, Teuchos::Copy);

      out << "  MV mvcopy2 = Tpetra::createCopy (*mvview);" << endl;
      MV mvcopy2 = Tpetra::createCopy (*mvview);

      out << "  MV mvcopy3 (map, numVecs, false); mvcopy3.randomize (); "
        "Tpetra::deep_copy (mvcopy3, *mvview);" << endl;
      MV mvcopy3 (mvview->getMap (), mvview->getNumVectors (), false);
      mvcopy3.randomize ();
      Tpetra::deep_copy (mvcopy3, *mvview);

      out << "Deep copy must preserve norms of the copied columns" << endl;

      // Deep copy must preserve norms of the copied columns.
      Array<Mag> ncopy (inds.size ());
      mvcopy.normInf (ncopy ());
      TEST_COMPARE_FLOATING_ARRAYS(ncopy,nsub,M0);
      out << "  Inf-norms of mvcopy's columns: " << toString (ncopy ()) << endl;

      Array<Mag> ncopy2 (inds.size ());
      mvcopy2.normInf (ncopy2 ());
      TEST_COMPARE_FLOATING_ARRAYS(ncopy2,nsub,M0);
      out << "  Inf-norms of mvcopy2's columns: " << toString (ncopy2 ()) << endl;

      Array<Mag> ncopy3 (inds.size ());
      mvcopy3.normInf (ncopy3 ());
      TEST_COMPARE_FLOATING_ARRAYS(ncopy3,nsub,M0);
      out << "  Inf-norms of mvcopy3's columns: " << toString (ncopy3 ()) << endl;

      out << "Test whether the copies are independent of their source, "
        "and of each other" << endl;

      // Change all the entries in both the view, and in all copies of
      // the view.  We will test below whether all of these
      // MultiVectors are independent.
      mvview->putScalar (ONE);
      mvcopy.putScalar (TWO);
      mvcopy2.putScalar (THREE);
      mvcopy3.putScalar (FOUR);

      Array<Mag> nsub_aft (inds.size ());
      mvview->normInf (nsub_aft ());
      out << "  Inf-norms of mvview's columns: " << toString (nsub_aft ()) << endl;
      Array<Mag> ones (inds.size (), STM::one ());
      TEST_COMPARE_FLOATING_ARRAYS(nsub_aft, ones, M0);

      Array<Mag> ncopy_aft (inds.size ());
      mvcopy.normInf (ncopy_aft ());
      out << "  Inf-norms of mvcopy's columns: " << toString (ncopy_aft ()) << endl;
      Array<Mag> twos (inds.size (), STM::one () + STM::one ());
      TEST_COMPARE_FLOATING_ARRAYS(ncopy_aft, twos, M0);

      Array<Mag> ncopy2_aft (inds.size ());
      mvcopy2.normInf (ncopy2_aft ());
      out << "  Inf-norms of mvcopy2's columns: " << toString (ncopy2_aft ()) << endl;
      Array<Mag> threes (inds.size (), STM::one () + STM::one () + STM::one ());
      TEST_COMPARE_FLOATING_ARRAYS(ncopy2_aft,threes,M0);

      Array<Mag> ncopy3_aft (inds.size ());
      mvcopy3.normInf (ncopy3_aft ());
      out << "  Inf-norms of mvcopy3's columns: " << toString (ncopy3_aft ()) << endl;
      Array<Mag> fours (inds.size (), STM::one () + STM::one () + STM::one () + STM::one ());
      TEST_COMPARE_FLOATING_ARRAYS(ncopy3_aft,fours,M0);
    }

    out << "Part 2:" << endl;
    {
      // create random MV
      MV morig (map,numVectors);
      morig.randomize ();
      // Test the copy constructor in its deep copy mode.
      MV mcopy1 (morig, Teuchos::Copy);
      // Test createCopy.  It should do the same as above, except that
      // mcopy2 always has view semantics, even in the (old)
      // KokkosClassic version of Tpetra.
      MV mcopy2 = createCopy (morig);

      // verify that all three have identical values
      Array<Mag> norig (numVectors);
      Array<Mag> ncopy1 (numVectors);
      Array<Mag> ncopy2 (numVectors);
      morig.normInf (norig);
      mcopy1.normInf (ncopy1);
      mcopy2.normInf (ncopy2);
      TEST_COMPARE_FLOATING_ARRAYS(norig, ncopy1, M0);
      TEST_COMPARE_FLOATING_ARRAYS(norig, ncopy2, M0);

      // Change all three MultiVectors to have different values.  The
      // three MultiVectors should be independent of each other.  That
      // is, none of them should be a view of any of the others.
      morig.putScalar (ZERO);
      mcopy1.putScalar (ONE);
      mcopy2.putScalar (TWO);

      // Check the above assertion about independence.
      Array<Mag> zeros (numVectors, STM::zero ());
      Array<Mag> ones (numVectors, STM::one ());
      Array<Mag> twos (numVectors, STM::one () + STM::one ());
      morig.normInf (norig);
      mcopy1.normInf (ncopy1);
      mcopy2.normInf (ncopy2);
      TEST_COMPARE_FLOATING_ARRAYS(norig, zeros, M0);
      TEST_COMPARE_FLOATING_ARRAYS(ncopy1, ones, M0);
      TEST_COMPARE_FLOATING_ARRAYS(ncopy2, twos, M0);
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Vector, CopyConst, LO , GO , Scalar , Node )
  {
    typedef Tpetra::Vector<Scalar,LO,GO,Node>       V;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;

    out << "Test: Vector, CopyConst" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,2,comm);
    // create random MV
    V morig(map);
    morig.randomize();
    // copy it
    V mcopy1(createCopy(morig)), mcopy2(createCopy(morig));
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

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Vector, Indexing, LO , GO , Scalar , Node )
  {
    typedef Tpetra::Vector<Scalar,LO,GO,Node>       V;
    typedef ScalarTraits<Scalar>              SCT;
    typedef typename SCT::magnitudeType Magnitude;

    out << "Test: Vector, Indexing" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,100,comm);
    // create two random Vector objects
    V v1(map), v2(map);
    v1.randomize();
    v2.randomize();
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
    err = v1.norm1();
    TEST_EQUALITY_CONST(err,SCT::zero());

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, SingleVecNormalize, LO , GO , Scalar , Node )
  {
    typedef Map<LO, GO, Node> map_type;

    out << "Test: MultiVector, SingleVecNormalize" << endl;
    Teuchos::OSTab tab0 (out);

    // this documents a usage case in Anasazi::SVQBOrthoManager, which was failing
    // error turned out to be a neglected return in both implementations of update(),
    // after passing the buck to scale() in the case of alpha==0 or beta==0 or gamma=0
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Magnitude M1  = ScalarTraits<Magnitude>::one();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVectors = 6;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));
    // create random MV
    MV mv(map,numVectors);
    mv.randomize();
    // compute the norms
    Array<Magnitude> norms(numVectors);
    mv.norm2(norms());
    for (size_t j=0; j<numVectors; ++j) {
      // get a view of column j, normalize it using update()
      RCP<MV> mvj = mv.subViewNonConst(tuple<size_t>(j));
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
    TEST_COMPARE_FLOATING_ARRAYS(norms,ones,testingTol<Scalar>()*errorTolSlack);

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, CountDot, LO , GO , Scalar , Node )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;

    out << "Test: MultiVector, CountDot" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
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
    MV mvec1(map,values(),2,numVectors),
       mvec2(map,values(),2,numVectors);

    // Make sure that MultiVector construction succeeded on all processes.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess,
                                  Teuchos::outArg (gblSuccess));
    if (gblSuccess) {
      out << "Successfully constructed MultiVector on all processes" << endl;
    } else {
      out << "FAILED to construct MultiVector on one or more processes" << endl;
      success = false;
      return;
    }

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

    // Make sure that the test passed on all processes, not just Proc 0.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, CountDotNonTrivLDA, LO , GO , Scalar , Node )
  {
    out << "Test dot products of MultiVectors created from an input "
      "Teuchos::ArrayView with nontrivial LDA." << endl;
    Teuchos::OSTab tab1 (out);

    // same as CountDot, but the A,LDA has a non-trivial LDA (i.e., LDA != myLen)
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    const size_t LDA = 3;
    RCP<const Map<LO,GO,Node> > map =
      createContigMapWithNode<LO,GO,Node> (INVALID, numLocal, comm);
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

    out << "Create the two MultiVectors" << endl;
    MV mvec1(map,values(),LDA,numVectors);
    MV mvec2(map,values(),LDA,numVectors);

    out << "Do the dot products" << endl;
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

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, CountNorm1, LO , GO , Scalar , Node )
  {
    typedef Tpetra::global_size_t GST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;

    out << "Test: MultiVector, CountNorm1" << endl;
    Teuchos::OSTab tab0 (out);

    const MT M0 = Teuchos::ScalarTraits<MT>::zero ();

    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm ();
    const int numImages = comm->getSize();

    // create a Map
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    RCP<const Map<LO, GO, Node> > map =
      createContigMapWithNode<LO, GO, Node> (INVALID, numLocal, comm);

    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    //                               [0 1 2]
    // norm1(values) = [0 2 4]
    // over all procs, this is [0 2*nprocs 4*nprocs]
    // mean is [0 1 2]
    Array<Scalar> values (6);
    values[0] = as<Scalar>(0);
    values[1] = as<Scalar>(0);
    values[2] = as<Scalar>(1);
    values[3] = as<Scalar>(1);
    values[4] = as<Scalar>(2);
    values[5] = as<Scalar>(2);
    MV mvec (map, values (), 2, numVectors);

    // Make sure that all processes successfully constructed the MultiVector.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess,
                                  Teuchos::outArg (gblSuccess));
    if (gblSuccess) {
      out << "Successfully constructed the MultiVector on all processes" << endl;
    } else {
      out << "FAILED to construct the MultiVector on some process" << endl;
      success = false;
      return;
    }

    // compute, check norms
    {
      Array<MT> norms(numVectors), answer(numVectors);
      answer[0] = as<MT>(0);
      answer[1] = as<MT>(2*numImages);
      answer[2] = as<MT>(4*numImages);
      mvec.norm1(norms());
      TEST_COMPARE_FLOATING_ARRAYS(norms,answer,M0);
    }
    {
      // compute, check means
      Array<Scalar> means(numVectors), answer(numVectors);
      mvec.meanValue(means());
      answer[0] = as<Scalar>(0);
      answer[1] = as<Scalar>(1);
      answer[2] = as<Scalar>(2);
      TEST_COMPARE_FLOATING_ARRAYS(means,answer,M0);
      for (size_t j=0; j < numVectors; ++j) {
        TEST_EQUALITY( mvec.getVector(j)->meanValue(), answer[j] );
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, CountNormInf, LO , GO , Scalar , Node )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;

    out << "Test: MultiVector, CountNormInf" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
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
    MV mvec(map,values(),2,numVectors);
    Array<MT> norms(numVectors), answer(numVectors);
    answer[0] = as<MT>(0);
    answer[1] = as<MT>(1);
    answer[2] = as<MT>(2);
    // do the dots
    mvec.normInf(norms());
    // check the answers
    TEST_COMPARE_FLOATING_ARRAYS(norms,answer,M0);

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, Norm2, LO , GO , Scalar , Node )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;

    out << "Test: MultiVector, Norm2" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 13;
    const size_t numVectors = 7;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    MV mvec(map,numVectors);
    // randomize the multivector
    mvec.randomize();
    // take norms; they should not be zero
    Array<MT> normsRand(numVectors), normsZero(numVectors);
    mvec.norm2(normsRand());
    // zero the vector
    mvec.putScalar(ScalarTraits<Scalar>::zero());
    // take norms; they should be zero
    mvec.norm2(normsZero());
    // check the answers
    bool local_success = true;
    for (size_t i=0; i<numVectors; ++i) {
      TEST_ARRAY_ELE_INEQUALITY(normsRand,i,M0);
      TEST_ARRAY_ELE_EQUALITY(normsZero,i,M0);
    }
    success &= local_success;

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, BadCombinations, LO , GO , Scalar , Node )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;

    out << "Test: MultiVector, BadCombinations" << endl;
    Teuchos::OSTab tab0 (out);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // create a Map
    const Scalar rnd = ScalarTraits<Scalar>::random();
    // two maps: one has two entires per node, the other disagrees on node 0
    RCP<const Map<LO,GO,Node> > map1 = createContigMapWithNode<LO,GO,Node>(INVALID,2,comm),
                                map2 = createContigMapWithNode<LO,GO,Node>(INVALID,myImageID == 0 ? 1 : 2,comm);
    // multivectors from different maps are incompatible for all ops
    // multivectors from the same map are compatible only if they have the same number of
    //    columns
    MV m1n1(map1,1), m1n2(map1,2), m2n2(map2,2), m1n2_2(map1,2);
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
    TEST_THROW(m1n2.reciprocal(m1n1), std::runtime_error);                  // reciprocal
    TEST_THROW(m1n2.reciprocal(m2n2), std::runtime_error);

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, Typedefs,        LO , GO , Scalar , Node )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;

    out << "Test: MultiVector, Typedefs" << endl;
    Teuchos::OSTab tab0 (out);

    TEST_EQUALITY_CONST( (std::is_same< scalar_type         , Scalar  >::value) == true, true );
    TEST_EQUALITY_CONST( (std::is_same< local_ordinal_type  , LO >::value) == true, true );
    TEST_EQUALITY_CONST( (std::is_same< global_ordinal_type , GO >::value) == true, true );
    TEST_EQUALITY_CONST( (std::is_same< node_type           , Node    >::value) == true, true );
  }

#if defined(HAVE_TEUCHOS_COMPLEX) && (defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) || defined(HAVE_TPETRA_INST_COMPLEX_FLOAT))
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, ComplexDotOneColumn, RealType, LO, GO, Node )
  {
    using Teuchos::rcp_implicit_cast;

    typedef RealType magnitude_type;
    typedef std::complex<RealType> scalar_type;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    typedef Tpetra::MultiVector<scalar_type, LO, GO, Node> MV;

    typedef Teuchos::SerialComm<int> comm_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    constexpr bool debug = true;

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test: MultiVector, ComplexDotOneColumn" << endl;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Create Comm and Map" << endl;

    // We use a SerialComm so that we don't complicate the test by
    // introducing interprocess communication.  The point of this test
    // is to test conjugation, that is, to make sure that the dot
    // product of x and y is the conjugate transpose of x times y.
    RCP<const Comm<int> > serialComm =
      rcp_implicit_cast<const comm_type > (rcp (new comm_type));

    const size_t numLocalElts = 1;
    const Tpetra::global_size_t numGlobalElts = serialComm->getSize () * numLocalElts;
    const GO indexBase = 0;
    RCP<const map_type> map (new map_type (numGlobalElts, indexBase, serialComm,
                                           Tpetra::GloballyDistributed));

    myOut << "Create MultiVectors x and y" << endl;
    MV x (map, 1);
    MV y (map, 1);

    std::vector<scalar_type> results (1, STS::zero ()); // dot product result

    myOut << "Modify entries of x and y" << endl;

    // dot([i], [i]) should be 1, not -1.
    x.replaceLocalValue (LO (0), 0, scalar_type (STM::zero (), STM::one ()));
    y.replaceLocalValue (LO (0), 0, scalar_type (STM::zero (), STM::one ()));

    myOut << "Compute dot product of x and y" << endl;
    x.dot (y, results);
    TEST_EQUALITY( results[0], STS::one() );

    myOut << "Modify entries of x and y" << endl;

    // dot([-i], [i]) should be -1, not +1.
    x.replaceLocalValue (LO (0), 0, scalar_type (STM::zero (), -STM::one ()));
    y.replaceLocalValue (LO (0), 0, scalar_type (STM::zero (), STM::one ()));

    myOut << "Compute dot product of x and y" << endl;
    x.dot (y, results);
    TEST_EQUALITY( results[0], -STS::one() );

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*getDefaultComm (), REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }
#endif // defined(HAVE_TEUCHOS_COMPLEX) && (defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) || defined(HAVE_TPETRA_INST_COMPLEX_FLOAT))

  // Test that MultiVector can be declared with no template
  // parameters, so that every template parameter has its default
  // value.
  TEUCHOS_UNIT_TEST( MultiVector, AllDefaultTemplateParameters )
  {
    // If you are letting all template parameters take their default
    // values, you must follow the class name MultiVector with <>.
    typedef MultiVector<> mv_type;
    typedef mv_type::scalar_type scalar_type;
    typedef mv_type::local_ordinal_type local_ordinal_type;
    typedef mv_type::global_ordinal_type global_ordinal_type;

    out << "Test: MultiVector, AllDefaultTemplateParameters" << std::endl;
    Teuchos::OSTab tab0 (out);

    // Verify that the default Scalar type is double.  We can't put
    // the std::is_same expression in the macro, since it has a comma
    // (commas separate arguments in a macro).
    const bool defaultScalarMatchesTpetra =
      std::is_same<scalar_type,
                   Tpetra::Details::DefaultTypes::scalar_type>::value;
    TEST_ASSERT( defaultScalarMatchesTpetra );

    // Verify that the default LocalOrdinal type is the same as Map's
    // default LocalOrdinal type.  This assumes that all of Map's
    // template parameters have default values.
    //
    // We can't put the std::is_same expression in the macro, since it has
    // a comma (commas separate arguments in a macro).
    typedef Tpetra::Map<>::local_ordinal_type map_local_ordinal_type;
    const bool defaultLocalOrdinalIsInt =
      std::is_same<local_ordinal_type, map_local_ordinal_type>::value;
    TEST_ASSERT( defaultLocalOrdinalIsInt );

    // Verify that the default GlobalOrdinal type has size no less
    // than the default LocalOrdinal type.  Currently (as of 17 Jun
    // 2014), the default GlobalOrdinal type is the same as the
    // default LocalOrdinal type, but at some point we may want to
    // change it to default to a 64-bit integer type.
    TEST_ASSERT( sizeof (global_ordinal_type) >= sizeof (local_ordinal_type) );

    // Make sure that the test passed on all processes, not just Proc 0.
    RCP<const Comm<int> > comm = getDefaultComm ();
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  // Test MultiVector::replaceMap.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, ReplaceMap, LO, GO, Scalar, Node )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using std::endl;
    typedef Tpetra::global_size_t GST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::Map<LO,GO,Node> map_type;

    out << "Test: MultiVector, ReplaceMap" << endl;
    Teuchos::OSTab tab0 (out);
    //
    // Create a Map, on which every process in the communicator has nonzero rows.
    //
    RCP<const Comm<int> > comm = getDefaultComm ();
    const size_t lclNumRows = 5;
    //const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GST gblNumRows = lclNumRows * static_cast<GST> (comm->getSize ());
    const GO indexBase = 0;
    RCP<const map_type> map = rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));
    TEST_ASSERT( ! map.is_null () );
    TEST_EQUALITY( map->getGlobalNumElements (), gblNumRows );
    TEST_EQUALITY( map->getLocalNumElements (), lclNumRows );
    TEST_EQUALITY( map->getIndexBase (), indexBase );

    //
    // Create a MultiVector with that Map.
    //
    const size_t numCols = 3;
    MV X (map, numCols);
    TEST_EQUALITY( X.getNumVectors (), numCols );
    TEST_ASSERT( ! X.getMap ().is_null () );
    if (! X.getMap ().is_null ()) {
      TEST_ASSERT( X.getMap ()->isSameAs (*map) );
    }

    //
    // Split the Comm, so that Proc 0 is in its own Comm, and the
    // other processes are in the other Comm.
    //
    const int color = (comm->getRank () == 0) ? 0 : 1;
    // Make the key the same on all processes with the same color, so
    // that MPI_Comm_split will order the processes in the split
    // communicators by their current rank.
    const int key = 0;
    RCP<const Comm<int> > subsetComm = comm->split (color, key);
    TEST_ASSERT( ! subsetComm.is_null () );
    if (! subsetComm.is_null ()) {
      TEST_ASSERT( (comm->getRank () == 0 && subsetComm->getSize () == 1) ||
                   (comm->getRank () != 0 && subsetComm->getSize () == comm->getSize () - 1) );
    }

    // Make a Map which exists on Processes 1 .. P-1, not on Proc 0.
    // On Proc 0, we pass in a null input Comm, to tell the Map that
    // we want to exclude that process.
    RCP<const map_type> subsetMap =
      map->replaceCommWithSubset (comm->getRank () == 0 ? Teuchos::null : subsetComm);
    TEST_ASSERT( (comm->getRank () == 0 && subsetMap.is_null ()) ||
                 (comm->getRank () != 0 && ! subsetMap.is_null ()) );

    //
    // Replace the MultiVector's original Map with a subset Map.
    //
    X.replaceMap (subsetMap);
    TEST_ASSERT( (comm->getRank () == 0 && X.getMap ().is_null ()) ||
                 (comm->getRank () != 0 && ! X.getMap ().is_null ()) );

    // The number of columns must not change, even on excluded processes.
    TEST_EQUALITY( X.getNumVectors (), numCols );

    if (comm->getRank () == 0) {
      TEST_EQUALITY( X.getLocalLength (), static_cast<size_t> (0) );
    }
    else { // my rank is not zero
      TEST_EQUALITY( X.getLocalLength (), lclNumRows );
      if (! subsetMap.is_null () && ! X.getMap ().is_null ()) {
        // This is a collective on the subset communicator.
        TEST_ASSERT( X.getMap ()->isSameAs (*subsetMap) );
      }
    }

    //
    // Replace the MultiVector's subset Map with its original Map.
    //
    X.replaceMap (map);
    TEST_ASSERT( ! X.getMap ().is_null () );
    if (! X.getMap ().is_null ()) {
      TEST_ASSERT( ! X.getMap ()->getComm ().is_null () );
    }
    TEST_EQUALITY( X.getNumVectors (), numCols );
    TEST_EQUALITY( X.getLocalLength (), lclNumRows );
    TEST_EQUALITY( X.getGlobalLength (), gblNumRows );


    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess) {
      out << "Test PASSED on all processes" << endl;
    } else {
      out << "Test FAILED on one or more processes" << endl;
      success = false;
    }
  }

  // Make sure that deep_copy compiles, and actually does a deep copy.
  //
  // NOTE: This test only exercises deep_copy for MVs of the same type.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, DeepCopy, LO, GO, Scalar, Node )
  {
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename MV::mag_type mag_type;

    out << "Test: MultiVector, DeepCopy" << endl;
    Teuchos::OSTab tab0 (out);

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLocal = 13;
    const size_t numVecs  = 7;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, getDefaultComm ()));

    Teuchos::Array<mag_type> norms (numVecs);

    MV X (map, numVecs);
    X.putScalar (static_cast<Scalar> (42.0));
    X.normInf (norms ());
    for (size_t j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms[j], static_cast<mag_type> (42.0) );
      norms[j] = Teuchos::ScalarTraits<mag_type>::zero ();
    }

    MV Y (X, Teuchos::Copy);
    Y.normInf (norms ());
    for (size_t j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms[j], static_cast<mag_type> (42.0) );
      norms[j] = Teuchos::ScalarTraits<mag_type>::zero ();
    }

    MV Z (map, numVecs);
    Tpetra::deep_copy (Z, X);
    Z.normInf (norms ());
    for (size_t j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms[j], static_cast<mag_type> (42.0) );
      norms[j] = Teuchos::ScalarTraits<mag_type>::zero ();
    }

    MV W (map, numVecs);
    W.update (STS::one (), Y, -STS::one (), Z, STS::zero ());
    W.normInf (norms ());
    for (size_t j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms[j], Teuchos::ScalarTraits<mag_type>::zero () );
    }
  }

  // Test Tpetra::MultiVector dual view semantics.
  //
  // Create a Tpetra::MultiVector X, and fill it with some
  // characteristic number.  Use getLocalView() to modify its data in
  // either memory space (we can test both sides here -- for devices
  // with a single memory space, that will just be a redundant test),
  // and sync.  Then, use an independent mechanism (that doesn't
  // involve Kokkos::DualView or Kokkos::View) to see whether the
  // Tpetra::MultiVector saw the change.
  //
  // This tests ensures that getLocalView() actually returns a view of
  // the data, NOT a deep copy.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, DualViewSemantics, LO, GO, Scalar, Node )
  {
    typedef Tpetra::global_size_t GST;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename MV::device_type device_type;

    out << "Test: MultiVector's DualView semantics" << endl;
    Teuchos::OSTab tab0 (out);

    int lclSuccess = 1;
    int gblSuccess = 1;
    std::ostringstream errStrm;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    RCP<const Comm<int> > comm = getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    const size_t numLclRows = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLclRows, indexBase, comm));

    const size_t numVecs = 3;
    RCP<MV> X;
    try {
      X = rcp (new MV (map, numVecs, false));
      lclSuccess = 1;
    }
    catch (std::exception& e) {
      errStrm << "Process " << myRank << ": MV constructor threw exception: "
              << e.what () << endl;
      lclSuccess = 0;
    }
    gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
    if (gblSuccess != 1) {
      out << "MV constructor threw an exception on one or more processes!" << endl;
      for (int r = 0; r < numProcs; ++r) {
        if (r == myRank) {
          std::cerr << errStrm.str ();
        }
        comm->barrier ();
        comm->barrier ();
        comm->barrier ();
      }
      return; // no sense in continuing.
    }

    // Don't use a negative number, in case Scalar is an unsigned integer.
    const Scalar ONE = STS::one ();
    const Scalar TWO = STS::one () + STS::one ();
    try {
      X->putScalar (TWO);
    }
    catch (std::exception& e) {
      errStrm << "Process " << myRank << ": MV::putScalar threw exception: "
              << e.what () << endl;
      lclSuccess = 0;
    }
    gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
    if (gblSuccess != 1) {
      out << "MV::putScalar threw an exception on one or more processes!" << endl;
      for (int r = 0; r < numProcs; ++r) {
        if (r == myRank) {
          std::cerr << errStrm.str ();
        }
        comm->barrier ();
        comm->barrier ();
        comm->barrier ();
      }
      return; // no sense in continuing.
    }

    // Modify the data through the host View, by setting all of its
    // entries to a different number than before.  (ONE and TWO differ
    // even in the finite field Z_2.)
    {
      auto X_lcl_h = X->getLocalViewHost(Tpetra::Access::OverwriteAll);
      Kokkos::deep_copy (X_lcl_h, ONE);
    }
    // Now compute the inf-norms of the columns of X.  (We want a
    // separate mechanism from methods that return Kokkos::DualView or
    // Kokkos::View.)  All inf-norms should be ONE, not TWO.
    {
      typedef typename MV::mag_type mag_type;
      Kokkos::DualView<mag_type*, device_type> norms ("norms", numVecs);
      norms.template modify<device_type> ();
      X->normInf (norms.template view<device_type> ());
      norms.sync_host ();
      for (size_t k = 0; k < numVecs; ++k) {
        TEST_EQUALITY_CONST( norms.h_view(k), ONE );
      }
    }
  }


  // Test constructor that takes a Kokkos::DualView.
  //
  // Create a Kokkos::DualView X_lcl, and fill it with some
  // characteristic number.  Create a Tpetra::MultiVector X_gbl that
  // views X_lcl, modify X_lcl's data, and sync.  Then, test whether
  // X_gbl saw the change.
  //
  // This tests whether the Tpetra::MultiVector constructor that takes
  // a Kokkos::DualView actually views the DualView.  (It must NOT
  // make a deep copy.)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, DualViewCtor, LO, GO, Scalar, Node )
  {
    typedef Tpetra::global_size_t GST;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
    typedef typename MV::device_type device_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    out << "Tpetra::MultiVector DualView constructor test" << endl;
    Teuchos::OSTab tab0 (out);

    const Scalar ONE = STS::one ();
    const Scalar TWO = ONE + ONE;
    //int lclSuccess = 1;
    //int gblSuccess = 1;

    // This typedef (a 2-D Kokkos::DualView specialization) must exist.
    typedef typename MV::dual_view_type dual_view_type;

    // We'll need this for error checking before we need it in Tpetra.
    RCP<const Comm<int> > comm = getDefaultComm ();

    // Create the Kokkos::DualView.
    const size_t numLclRows = 10;
    const size_t numVecs = 3;
    dual_view_type X_lcl ("X_lcl", numLclRows, numVecs);

    // Modify the Kokkos::DualView's data on the host.
    {
      auto X_lcl_h = X_lcl.view_host ();
      X_lcl.modify_host ();
      Kokkos::deep_copy (X_lcl_h, ONE);
      X_lcl.template sync<device_type> ();
    }

    // Hand off the Kokkos::DualView to a Tpetra::MultiVector.
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLclRows, indexBase, comm));
    MV X_gbl (map, X_lcl);

    // Make sure (using an independent mechanism, in this case the
    // inf-norm) that X_gbl's constructor didn't change the values in
    // X_lcl.
    typedef typename MV::mag_type mag_type;
    Kokkos::DualView<mag_type*, device_type> norms ("norms", numVecs);
    norms.template modify<device_type> ();
    X_gbl.normInf (norms.template view<device_type> ());
    norms.sync_host ();
    for (size_t k = 0; k < numVecs; ++k) {
      TEST_EQUALITY_CONST( norms.h_view(k), ONE );
    }

    // Now change the values in X_lcl.  X_gbl should see them.  Just
    // for variety, we do this on the device, not on the host.
    {
      auto X_lcl_d = X_lcl.template view<device_type> ();
      X_lcl.template modify<device_type> ();
      Kokkos::deep_copy (X_lcl_d, TWO);
      X_lcl.sync_host ();
    }

    // Make sure that X_gbl saw the changes made to X_lcl's data.
    norms.template modify<device_type> ();
    X_gbl.normInf (norms.template view<device_type> ());
    norms.sync_host ();
    for (size_t k = 0; k < numVecs; ++k) {
      TEST_EQUALITY_CONST( norms.h_view(k), TWO );
    }
  }


  // Test constructor that takes a Kokkos::View (on device).
  //
  // Create a Kokkos::View X_lcl, and fill it with some characteristic
  // number.  Create a Tpetra::MultiVector X_gbl that views X_lcl, and
  // modify X_lcl's data.  Then, test whether X_gbl saw the change.
  //
  // This tests whether the Tpetra::MultiVector constructor that takes
  // a Kokkos::View actually views the Kokkos::View.  (It must NOT
  // make a deep copy.)
  //
  // NOTE: It is undefined for users to modify the device View without
  // respecting Tpetra::MultiVector's DualView semantics.  That is,
  // they need to use modify() and sync() correctly for the
  // Tpetra::MultiVector (or the underlying Kokkos::DualView).
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, ViewCtor, LO, GO, Scalar, Node )
  {
    typedef Tpetra::global_size_t GST;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename MV::device_type device_type;
    // This typedef (a 2-D Kokkos::DualView specialization) must exist.
    typedef typename MV::dual_view_type dual_view_type;

    out << "Tpetra::MultiVector View constructor test" << endl;
    Teuchos::OSTab tab0 (out);

    const Scalar ONE = STS::one ();
    const Scalar TWO = ONE + ONE;
#define TPETRA_MULTIVECTOR_VIEWCTOR_DO_NOT_TEST
#if !defined(TPETRA_MULTIVECTOR_VIEWCTOR_DO_NOT_TEST)
    const Scalar THREE = TWO + ONE;
#endif
    int lclSuccess = 1;
    int gblSuccess = 0; // to be set below

    // We'll need this for error checking before we need it in Tpetra.
    RCP<const Comm<int> > comm = getDefaultComm ();

    // Create the Kokkos::View X_lcl.
    const size_t numLclRows = 10;
    const size_t numVecs = 3;

    /// KJ : release local object, this workflow is problematic. 
    ///      a user create a device view and hand it to tpetra. 
    ///      tpetra now has unmatched referecne count for host and device view
    ///      as the local device view is alive. this is the case that we do not want 
    ///      to encourage users.
    typename dual_view_type::t_dev X_lcl ("X_lcl", numLclRows, numVecs);

    // Modify the Kokkos::View's data.
    Kokkos::deep_copy (X_lcl, ONE);

    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST(gblSuccess, 1);
      if (gblSuccess != 1) {
        return;
      }
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": checkpoint 1" << std::endl;
      std::cerr << os.str ();
    }

    // Hand off the Kokkos::View to a Tpetra::MultiVector.
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLclRows, indexBase, comm));

    /// KJ : release local object, this workflow is problematic. 
    ///      a user create a device view and hand it to tpetra. 
    ///      tpetra now has unmatched referecne count for host and device view
    ///      as the local device view is alive. this is the case that we do not want 
    ///      to encourage users.
    MV X_gbl (map, X_lcl);
    
    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST(gblSuccess, 1);
      if (gblSuccess != 1) {
        return;
      }
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": checkpoint 2" << std::endl;
      std::cerr << os.str ();
    }

    // Make sure (using an independent mechanism, in this case the
    // inf-norm) that X_gbl's constructor didn't change the values in
    // X_lcl.
    typedef typename MV::mag_type mag_type;
    Kokkos::DualView<mag_type*, device_type> norms ("norms", numVecs);
    norms.template modify<device_type> ();
    X_gbl.normInf (norms.template view<device_type> ());
    norms.sync_host ();
    for (size_t k = 0; k < numVecs; ++k) {
      TEST_EQUALITY_CONST( norms.h_view(k), ONE );
    }

    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST(gblSuccess, 1);
      if (gblSuccess != 1) {
        return;
      }
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": checkpoint 3" << std::endl;
      std::cerr << os.str ();
    }

    // Now change the values in X_lcl.  X_gbl should see them.  Be
    // sure to tell X_gbl that we want to modify its data on device.
    Kokkos::deep_copy (X_lcl, TWO);

    // Tpetra::MultiVector::normInf _should_ either read from the most
    // recently modified memory space, or do a sync to device first.
    // In either case, we don't need to do an explicit sync of X_gbl
    // before calling normInf.

    // Make sure that X_gbl saw the changes made to X_lcl's data.
    norms.template modify<device_type> ();
    X_gbl.normInf (norms.template view<device_type> ());
    norms.sync_host ();
    for (size_t k = 0; k < numVecs; ++k) {
      TEST_EQUALITY_CONST( norms.h_view(k), TWO );
    }

    {
      std::ostringstream os;
      os << ">>> Proc " << comm->getSize ();
      os << ": X_gbl.modified_host: " << (X_gbl.need_sync_device()?1:0)
         << ", X_gbl.modified_device: " << (X_gbl.need_sync_host()?1:0);
      os << std::endl;
      std::cerr << os.str ();
    }

    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST(gblSuccess, 1);
      if (gblSuccess != 1) {
        return;
      }
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": checkpoint 4" << std::endl;
      std::cerr << os.str ();
    }

#if !defined(TPETRA_MULTIVECTOR_VIEWCTOR_DO_NOT_TEST)
    // Just as X_gbl views X_lcl, X_lcl should also view X_gbl.  Thus,
    // if we modify X_gbl in host memory, and sync to device memory,
    // X_lcl should also be changed.

    // We modified on device above, and we're about to modify on host
    // now, so we need to sync to host first.
    auto X_host = X_gbl.getLocalViewHost(Tpetra::Access::ReadWrite);
    
    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST(gblSuccess, 1);
      if (gblSuccess != 1) {
        return;
      }
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": checkpoint 5" << std::endl;
      std::cerr << os.str ();
    }
    
    Kokkos::deep_copy (X_host, THREE);
    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST(gblSuccess, 1);
      if (gblSuccess != 1) {
        return;
      }
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": checkpoint 6" << std::endl;
      std::cerr << os.str ();
    }

    // FIXME (mfh 01 Mar 2015) We avoid writing a separate functor to
    // check the contents of X_lcl, by copying to host and checking
    // there.  Once Tpetra can use C++11, we should instead check on
    // device, by using a parallel_reduce functor.
    typename dual_view_type::t_dev::HostMirror X_lcl_host =
      Kokkos::create_mirror_view (X_lcl);
    Kokkos::deep_copy(X_lcl_host,X_lcl);

    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST(gblSuccess, 1);
      if (gblSuccess != 1) {
        return;
      }
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": checkpoint 7" << std::endl;
      std::cerr << os.str ();
    }

    bool same = true;
    for (size_t j = 0; j < numVecs; ++j) {
      for (size_t i = 0; i < numLclRows; ++i) {
        if (X_lcl_host(i,j) != X_host(i,j)) {
          same = false;
          break;
        }
      }
    }
    lclSuccess = same ? 1 : 0;
    gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
    if (gblSuccess != 1) {
      out << "We modified X_gbl in host memory, and sync'd to device memory, "
        "but X_lcl did not change!" << endl;
    }

    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST(gblSuccess, 1);
      if (gblSuccess != 1) {
        return;
      }
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": DONE" << std::endl;
      std::cerr << os.str ();
    }
#endif
#undef TPETRA_MULTIVECTOR_VIEWCTOR_DO_NOT_TEST
  }

// Macro used inside the SubViewSomeZeroRows test below.  It tests for
// global error, and if so, prints each process' error message and
// quits the test early.
//
// 'out' only prints on Process 0.  It's really not OK for other
// processes to print to stdout, but it usually works and we need to
// do it for debugging.
#define SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( WHAT_STRING ) do { \
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess)); \
  TEST_EQUALITY_CONST( gblSuccess, 1 ); \
  if (gblSuccess != 1) { \
    out << WHAT_STRING << " FAILED on one or more processes!" << endl; \
    for (int p = 0; p < numProcs; ++p) { \
      if (myRank == p && lclSuccess != 1) { \
        std::cout << errStrm.str () << std::flush; \
      } \
      comm->barrier (); \
      comm->barrier (); \
      comm->barrier (); \
    } \
    return; \
  } \
} while (false)

  // Exercise getVector, subView(Range1D) and subCopy(Range1D) where
  // some processes have zero rows.  Contributed by Andrew Bradley.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, SubViewSomeZeroRows, LO, GO, ST, Node )
  {
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::Vector<ST, LO, GO, Node> V;
    typedef Tpetra::MultiVector<ST, LO, GO, Node> MV;

    out << "Tpetra::MultiVector: Test subView and subCopy when some processes "
      "have zero rows" << endl;

    int lclSuccess = 1;
    int gblSuccess = 1;
    std::ostringstream errStrm; // for error collection

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    // Create a Map that puts everything on Process 0 and nothing on
    // the other processes.
    const Tpetra::global_size_t gblNumInds = 10;
    const size_t lclNumInds = (myRank == 0) ? 10 : 0;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (gblNumInds, lclNumInds, indexBase, comm));
    // Create a MultiVector with this Map.  Give it at least three
    // columns, so that when we try to take a subview with two
    // columns, it's actually a nontrivial subview.  (subView might
    // have an optimization when the input column range is exactly the
    // original set of columns.)
    const size_t origNumVecs = 5;
    MV mv (map, origNumVecs);

    // Make sure that getVector(NonConst) works on this MultiVector.
    // While doing that, fill each column with data that distinguish
    // it from the other columns.
    RCP<V> v_j;
    try {
      for (size_t j = 0; j < mv.getNumVectors (); ++j) {
        v_j = mv.getVectorNonConst (j);
        // Use j+1, so that no column gets filled with zeros.  That
        // will distinguish the zeroth column from a zero-filled
        // (Multi)Vector.
        v_j->putScalar (static_cast<ST> (j+1));
      }
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.getVector(j), "
        "mv.getVectorNonConst(j), or v_j->putScalar() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.getVector(j), mv.getVectorNonConst(j), or v_j->putScalar()" );

    // Make sure that every column got the right data.  Don't use
    // getVectorNonConst(j) to test this; we need independent
    // confirmation.
    try {
      Array<typename MV::mag_type> norms (mv.getNumVectors ());
      mv.normInf (norms ());
      for (size_t j = 0; j < mv.getNumVectors (); ++j) {
        TEST_EQUALITY( norms[j], static_cast<ST> (j+1) );
      }
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.normInf() threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.normInf()" );

    // Test getting a subView of the first two columns.
    out << "Test subView(Range1D(0,1))" << endl;
    Teuchos::Range1D r (0, 1);
    RCP<const MV> mv_sv;
    try {
      mv_sv = mv.subView (r);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.subView(Range1D(0,1)) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.subView(Range1D(0,1))" );

    // Make sure that the two columns being viewed actually have the
    // right data in them.
    try {
      Array<typename MV::mag_type> norms (mv_sv->getNumVectors ());
      mv_sv->normInf (norms ());
      TEST_EQUALITY( norms[0], static_cast<ST> (1) );
      TEST_EQUALITY( norms[1], static_cast<ST> (2) );
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv_sv->normInf() "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv_sv->normInf()" );

    // Make sure that the supposed "view" is actually a view.  That
    // is, if I change the original MultiVector's data, the view
    // should see the change right away.
    try {
      // Even if ST is bool, at least one of the columns will see the
      // change (mod it by 2 to see why).
      mv.putScalar (static_cast<ST> (mv.getNumVectors ()));
      Array<typename MV::mag_type> norms (mv_sv->getNumVectors ());
      mv_sv->normInf (norms ());
      TEST_EQUALITY_CONST( norms[0], static_cast<ST> (mv.getNumVectors ()) );
      TEST_EQUALITY_CONST( norms[1], static_cast<ST> (mv.getNumVectors ()) );
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.putScalar() or mv_sv->normInf() "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.putScalar() or mv_sv->normInf()" );

    // Restore the MultiVector's original data.
    try {
      for (size_t j = 0; j < mv.getNumVectors (); ++j) {
        v_j = mv.getVectorNonConst (j);
        v_j->putScalar (static_cast<ST> (j+1));
      }
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.getVectorNonConst(j) or "
        "v_j->putScalar() threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.getVectorNonConst(j) or v_j->putScalar()" );

    // Test subCopy (which reportedly has the same issue as subView).
    out << "Test subCopy(Range1D(0,1))" << endl;
    RCP<const MV> mv_sc;
    try {
      mv_sc = mv.subCopy (r);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.subCopy(Range1D(0,1)) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.subCopy(Range1D(0,1))" );

    // Make sure that the two copied columns actually have the right
    // data in them.
    try {
      Array<typename MV::mag_type> norms (mv_sv->getNumVectors ());
      mv_sc->normInf (norms ());
      TEST_EQUALITY( norms[0], static_cast<ST> (1) );
      TEST_EQUALITY( norms[1], static_cast<ST> (2) );
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv_sc->normInf() "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv_sc->normInf()" );

    // Make sure that the supposed "copy" is actually a (deep) copy.
    // That is, if I change the original MultiVector's data, the copy
    // should NOT see the change.
    try {
      // Even if ST is bool, at least one of the columns will see the
      // change (mod it by 2 to see why).
      mv.putScalar (static_cast<ST> (mv.getNumVectors ()));
      Array<typename MV::mag_type> norms (mv_sc->getNumVectors ());
      mv_sc->normInf (norms ());
      TEST_EQUALITY_CONST( norms[0], static_cast<ST> (1) );
      TEST_EQUALITY_CONST( norms[1], static_cast<ST> (2) );
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.putScalar() or mv_sc->normInf() "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.putScalar() or mv_sc->normInf()" );

    // Restore the MultiVector's original data.
    try {
      for (size_t j = 0; j < mv.getNumVectors (); ++j) {
        v_j = mv.getVectorNonConst (j);
        v_j->putScalar (static_cast<ST> (j+1));
      }
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.getVectorNonConst(j) or "
        "v_j->putScalar() threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.getVectorNonConst(j) or v_j->putScalar()" );

    // Test getting a subView of just the first column.
    out << "Test subView(Range1D(1,1))" << endl;
    Teuchos::Range1D r11 (1, 1);
    try {
      mv_sv = mv.subView (r11);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.subView(Range1D(1,1)) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.subView(Range1D(1,1))" );

    // Single-column subviews must always have constant stride.
    {
      bool constStride = false;
      try {
        constStride = ! mv_sv.is_null () && mv_sv->isConstantStride ();
      } catch (std::exception& e) {
        lclSuccess = 0;
        errStrm << "Process " << myRank << ": mv_sv->isConstantStride() "
          "threw exception: " << e.what () << endl;
      }
      SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv_sv->isConstantStride()" );
      lclSuccess = constStride ? 1 : 0;
      SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv_sv having constant stride" );
    }

    // Test subCopy of just the first column.
    out << "Test subCopy(Range1D(1,1))" << endl;
    try {
      mv_sc = mv.subCopy (r11);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.subCopy(Range1D(1,1)) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.subCopy(Range1D(1,1))" );

    // Single-column subcopies must always have constant stride.
    {
      bool constStride = false;
      try {
        constStride = ! mv_sc.is_null () && mv_sc->isConstantStride ();
      } catch (std::exception& e) {
        lclSuccess = 0;
        errStrm << "Process " << myRank << ": mv_sc->isConstantStride() "
          "threw exception: " << e.what () << endl;
      }
      SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv_sc->isConstantStride()" );
      lclSuccess = constStride ? 1 : 0;
      SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv_sc having constant stride" );
    }

    //
    // Make a noncontiguous subview of the original MultiVector.  Test
    // that both multiple-column and single-column subviews and
    // subcopies work.
    //

    // Start with a noncontiguous subview of the original MV.
    out << "Test subView([0, 2, 4])" << endl;
    RCP<const MV> X_noncontig;
    try {
      Array<size_t> colsToView (3);
      colsToView[0] = 0;
      colsToView[1] = 2;
      colsToView[2] = 4;
      X_noncontig = mv.subView (colsToView);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": mv.subView([0, 2, 4]) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "mv.subView([0, 2, 4])" );

    // Test getting a multiple-column noncontiguous subview of the
    // noncontiguous subview.
    out << "Test multi-column noncontig subview of noncontig subview" << endl;
    try {
      // View columns 0 and 2 of X_noncontig, which should be columns
      // 0 and 4 of the original MV.
      Array<size_t> colsToView (2);
      colsToView[0] = 0;
      colsToView[1] = 2;
      mv_sv = X_noncontig->subView (colsToView);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X_noncontig->subView([0, 2]) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X_noncontig->subView([0, 2])" );

    // Test getting a multiple-column noncontiguous subcopy of the
    // noncontiguous subview.
    out << "Test multi-column noncontig subcopy of noncontig subview" << endl;
    try {
      // Copy columns 0 and 2 of X_noncontig, which should be columns
      // 0 and 4 of the original MV.
      Array<size_t> colsToCopy (2);
      colsToCopy[0] = 0;
      colsToCopy[1] = 2;
      mv_sc = X_noncontig->subCopy (colsToCopy);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X_noncontig->subCopy([0, 2]) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X_noncontig->subCopy([0, 2])" );

    // Test getting a single-column subview of the noncontiguous
    // subview, using subView(Teuchos::ArrayView<const size_t>).
    out << "Test single-column noncontig subview of noncontig subview, "
      "using Teuchos::ArrayView<const size_t>" << endl;
    try {
      // View column 2 of X_noncontig, which should be column 4 of the
      // original MV.
      Array<size_t> colsToView (1);
      colsToView[0] = 2;
      mv_sv = X_noncontig->subView (colsToView);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X_noncontig->subView([2]) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X_noncontig->subView([2])" );

    // Test getting a single-column subview of the noncontiguous
    // subview, using subView(Teuchos::Range1D).
    out << "Test single-column noncontig subview of noncontig subview, "
      "using Teuchos::Range1D" << endl;
    try {
      // View column 2 of X_noncontig, which should be column 4 of the
      // original MV.
      mv_sv = X_noncontig->subView (Teuchos::Range1D (2, 2));
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X_noncontig->subView(Range1D(2,2)) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X_noncontig->subView(Range1D(2,2))" );

    // Test getting a single-column subcopy of the noncontiguous
    // subview, using subCopy(Teuchos::ArrayView<const size_t>).
    out << "Test single-column noncontig subcopy of noncontig subview, "
      "using Teuchos::ArrayView<const size_t>" << endl;
    try {
      // Copy column 2 of X_noncontig, which should be column 4 of the
      // original MV.
      Array<size_t> colsToCopy (1);
      colsToCopy[0] = 2;
      mv_sv = X_noncontig->subCopy (colsToCopy);
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X_noncontig->subCopy([2]) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X_noncontig->subCopy([2])" );

    // Test getting a single-column subview of the noncontiguous
    // subview, using subCopy(Teuchos::Range1D).
    out << "Test single-column noncontig subview of noncontig subview, "
      "using Teuchos::Range1D" << endl;
    try {
      // Copy column 2 of X_noncontig, which should be column 4 of the
      // original MV.
      mv_sc = X_noncontig->subCopy (Teuchos::Range1D (2, 2));
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X_noncontig->subCopy(Range1D(2,2)) "
        "threw exception: " << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X_noncontig->subCopy(Range1D(2,2))" );
  }


  // Create a MultiVector with zero rows on some processes, but a
  // nonzero number of columns.  Make sure that getLocalLength(),
  // getGlobalLength(), and getNumVectors() return the correct values.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, DimsWithSomeZeroRows, LO, GO, ST, Node )
  {
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<ST, LO, GO, Node> MV;
    typedef Tpetra::global_size_t GST;

    out << "Tpetra::MultiVector: Test MultiVector dimensions when some "
      "processes have zero rows" << endl;

    int lclSuccess = 1;
    int gblSuccess = 1;
    std::ostringstream errStrm; // for error collection

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    // Create a Map that puts nothing on Process 0 and something on
    // the other processes.
    const size_t lclNumRowsWhenPopulated = 1;
    const size_t lclNumRows = (myRank == 0) ?
      size_t (0) :
      lclNumRowsWhenPopulated;
    const GST gblNumRows = (numProcs == 1) ?
      GST (0) :
      GST ((numProcs - 1) * lclNumRowsWhenPopulated);
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

    const size_t numCols = 3;
    MV X (map, numCols);

    size_t reportedNumCols = 0;
    try {
      reportedNumCols = X.getNumVectors ();
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X.getNumVectors() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X.getNumVectors() threw exception" );
    if (reportedNumCols != numCols) {
      lclSuccess = 0;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "reportedNumCols != numCols" );

    size_t reportedLclNumRows = 0;
    try {
      reportedLclNumRows = X.getLocalLength ();
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X.getNumVectors() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X.getNumVectors() threw exception" );
    if (reportedLclNumRows != lclNumRows) {
      lclSuccess = 0;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "reportedLclNumRows != lclNumRows" );

    size_t reportedGblNumRows = 0;
    try {
      reportedGblNumRows = X.getGlobalLength ();
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X.getNumVectors() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X.getNumVectors() threw exception" );
    if (reportedGblNumRows != gblNumRows) {
      lclSuccess = 0;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "reportedGblNumRows != gblNumRows" );
  }


  // Create a MultiVector with zero rows on ALL processes, but a
  // nonzero number of columns.  Make sure that getLocalLength(),
  // getGlobalLength(), and getNumVectors() return the correct values.
  // Then, do the same thing with a globally 0 x 0 multivector.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, DimsWithAllZeroRows, LO, GO, ST, Node )
  {
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<ST, LO, GO, Node> MV;
    typedef Tpetra::global_size_t GST;
    typedef Tpetra::MultiVector<ST, LO, GO, Node> MV;

    out << "Tpetra::MultiVector: Test MultiVector dimensions when ALL "
      "processes have zero rows" << endl;

    int lclSuccess = 1;
    int gblSuccess = 1;
    std::ostringstream errStrm; // for error collection

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    // Create a Map that puts nothing on Process 0 and something on
    // the other processes.
    const size_t lclNumRows = 0;
    const GST gblNumRows = GST (numProcs * lclNumRows);
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

    size_t numCols = 3;
    MV X (map, numCols);

    size_t reportedNumCols = 0;
    try {
      reportedNumCols = X.getNumVectors ();
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X.getNumVectors() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X.getNumVectors() threw exception" );
    if (reportedNumCols != numCols) {
      lclSuccess = 0;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "reportedNumCols != numCols" );

    size_t reportedLclNumRows = 0;
    try {
      reportedLclNumRows = X.getLocalLength ();
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X.getNumVectors() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X.getNumVectors() threw exception" );
    if (reportedLclNumRows != lclNumRows) {
      lclSuccess = 0;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "reportedLclNumRows != lclNumRows" );

    size_t reportedGblNumRows = 0;
    try {
      reportedGblNumRows = X.getGlobalLength ();
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": X.getNumVectors() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "X.getNumVectors() threw exception" );
    if (reportedGblNumRows != gblNumRows) {
      lclSuccess = 0;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "reportedGblNumRows != gblNumRows" );

    numCols = 0;
    MV Y (map, 0);

    reportedNumCols = 0;
    try {
      reportedNumCols = Y.getNumVectors ();
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": Y.getNumVectors() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "Y.getNumVectors() threw exception" );
    if (reportedNumCols != numCols) {
      lclSuccess = 0;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "reportedNumCols != numCols" );

    reportedLclNumRows = 0;
    try {
      reportedLclNumRows = Y.getLocalLength ();
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": Y.getNumVectors() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "Y.getNumVectors() threw exception" );
    if (reportedLclNumRows != lclNumRows) {
      lclSuccess = 0;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "reportedLclNumRows != lclNumRows" );

    reportedGblNumRows = 0;
    try {
      reportedGblNumRows = Y.getGlobalLength ();
    } catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << myRank << ": Y.getNumVectors() threw exception: "
              << e.what () << endl;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "Y.getNumVectors() threw exception" );
    if (reportedGblNumRows != gblNumRows) {
      lclSuccess = 0;
    }
    SUBVIEWSOMEZEROROWS_REPORT_GLOBAL_ERR( "reportedGblNumRows != gblNumRows" );
  }

  // Swap test
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, Swap, LO , GO , Scalar , Node ) {
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar,LO, GO, Node> MV;
    typedef Tpetra::global_size_t GST;

    Scalar ONE  = Teuchos::ScalarTraits<Scalar>::one();
    Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int numProcs = comm->getSize ();

    // Create a Map that puts nothing on Process 0 and something on
    // the other processes.
    const size_t lclNumRows = 4;
    const GST gblNumRows = GST (numProcs * lclNumRows);
    const GO indexBase = 0;
    RCP<const map_type> map = rcp (new map_type (gblNumRows, indexBase, comm));

    size_t numCols = 3;
    MV Xo (map, numCols), Yo (map, numCols), Xn (map, numCols), Yn (map, numCols);

    // Comparison vectors (unswapped)
    Xo.putScalar(ZERO); Yo.putScalar(ONE);

    // Swapping vectors (swapped)
    Yn.putScalar(ZERO); Xn.putScalar(ONE);
    Xn.swap(Yn);

    // Compare
    TEST_COMPARE_FLOATING_ARRAYS(Xo.get1dView(),Xn.get1dView(),testingTol<Scalar>() * errorTolSlack);
    TEST_COMPARE_FLOATING_ARRAYS(Yo.get1dView(),Yn.get1dView(),testingTol<Scalar>() * errorTolSlack);

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, DualViewRefcountCheck, LO , GO , Scalar , Node ) {
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar,LO, GO, Node> MV;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    RCP<const map_type> map = rcp (new map_type (100, 0, comm));
    MV x(map, 1);

    const void* devicePtr = x.getLocalViewDevice(Tpetra::Access::ReadOnly).data();
    const void* hostPtr = x.getLocalViewHost(Tpetra::Access::ReadOnly).data();

    if(devicePtr != hostPtr)
    {
      //Host and device views are not the same.
      //Make sure that (assuming the 'device' space is not host accessible) checking out
      //host and device views at the same time is not allowed.
      bool threw = false;
      try
      {
        auto xDevice = x.getLocalViewDevice(Tpetra::Access::ReadOnly);
        //this shouldn't be allowed, since xDevice holds a reference to device view.
        auto xHost = x.getLocalViewHost(Tpetra::Access::ReadOnly);
      }
      catch(...)
      {
        threw = true;
      }
      TEST_EQUALITY(threw, true);
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, CopyCounterCheck, LO , GO , Scalar , Node ) {
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar,LO, GO, Node> MV;
    using device_view = typename MV::dual_view_type::t_dev;
    using host_view   = typename MV::dual_view_type::t_host;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    RCP<const map_type> map = rcp (new map_type (100, 0, comm));
    MV x(map, 1);
    x.putScalar(Teuchos::ScalarTraits<Scalar>::one());

    const device_view x_d = x.getLocalViewDevice(Tpetra::Access::ReadWrite);

    host_view y_h = create_mirror_view(x_d);  

    size_t correct_count;
    // Check to see if we'll be deep_copy-ing between memory spaces
    if(std::is_same<typename device_view::memory_space,typename host_view::memory_space>::value) {
      correct_count = 0;
    }
    else {
      correct_count = 1;
    }


    // Stop / Start  (reset first to clear counts from previous unit test calls)   
    Tpetra::Details::DeepCopyCounter::reset();   
    Tpetra::Details::DeepCopyCounter::start();
    Kokkos::deep_copy(y_h,x_d);
    Tpetra::Details::DeepCopyCounter::stop();   
    size_t count = Tpetra::Details::DeepCopyCounter::get_count_different_space();   
    TEST_EQUALITY(count,correct_count);


    // Reset / get_count (should be zero now)
    Tpetra::Details::DeepCopyCounter::reset();   
    count = Tpetra::Details::DeepCopyCounter::get_count_different_space();   
    TEST_EQUALITY(count,0);


    // Second  Stop / Start (should have the original count)
    Tpetra::Details::DeepCopyCounter::start();
    Kokkos::deep_copy(y_h,x_d);
    Tpetra::Details::DeepCopyCounter::stop();   
    count = Tpetra::Details::DeepCopyCounter::get_count_different_space();   
    TEST_EQUALITY(count,correct_count);


    // This guy should not get counted, since the counter is stopped
    Kokkos::deep_copy(y_h,x_d);
    count = Tpetra::Details::DeepCopyCounter::get_count_different_space();   
    TEST_EQUALITY(count,correct_count);


    // Third Second  Stop / Start (should have double the original count)
    Tpetra::Details::DeepCopyCounter::start();
    Kokkos::deep_copy(y_h,x_d);
    Tpetra::Details::DeepCopyCounter::stop();   
    count = Tpetra::Details::DeepCopyCounter::get_count_different_space();   
    TEST_EQUALITY(count,2*correct_count);
          
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, FenceCounterCheck, LO , GO , Scalar , Node ) {
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    auto exec_space = typename Node::execution_space();
    const std::string space = exec_space.name();

    /***********************************************************************/
    // Global fences
    size_t global_correct_count=1;

    // Stop / Start  (reset first to clear counts from previous unit test calls)   
    Tpetra::Details::FenceCounter::reset();   
    Tpetra::Details::FenceCounter::start();
    Kokkos::fence();
    Tpetra::Details::FenceCounter::stop();   
    size_t global_count = Tpetra::Details::FenceCounter::get_count_global(space);   
    size_t instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,global_correct_count);
    TEST_EQUALITY(instance_count,0);

    // Reset / get_count (should be zero now)
    Tpetra::Details::FenceCounter::reset();   
    global_count =Tpetra::Details::FenceCounter::get_count_global(space);   
    instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,0);
    TEST_EQUALITY(instance_count,0);

    // Second  Stop / Start (should have the original count)
    Tpetra::Details::FenceCounter::start();
    Kokkos::fence();
    Tpetra::Details::FenceCounter::stop();   
    global_count =Tpetra::Details::FenceCounter::get_count_global(space);   
    instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,global_correct_count);
    TEST_EQUALITY(instance_count,0);

    // This guy should not get counted, since the counter is stopped
    Kokkos::fence();
    global_count =Tpetra::Details::FenceCounter::get_count_global(space);   
    instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,global_correct_count);
    TEST_EQUALITY(instance_count,0);

    // Third Second  Stop / Start (should have double the original count)
    Tpetra::Details::FenceCounter::start();
    Kokkos::fence();
    Tpetra::Details::FenceCounter::stop();   
    global_count =Tpetra::Details::FenceCounter::get_count_global(space);   
    instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,2*global_correct_count);
    TEST_EQUALITY(instance_count,0);

    /***********************************************************************/
    // Instance Fences
    size_t instance_correct_count = 1;

    // Stop / Start  (reset first to clear counts from previous unit test calls)   
    Tpetra::Details::FenceCounter::reset();   
    Tpetra::Details::FenceCounter::start();
    exec_space.fence();
    Tpetra::Details::FenceCounter::stop();   
    global_count =Tpetra::Details::FenceCounter::get_count_global(space);   
    instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,0);
    TEST_EQUALITY(instance_count,instance_correct_count);

    // Reset / get_count (should be zero now)
    Tpetra::Details::FenceCounter::reset();   
    global_count =Tpetra::Details::FenceCounter::get_count_global(space);   
    instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,0);
    TEST_EQUALITY(instance_count,0);

    // Second  Stop / Start (should have the original count)
    Tpetra::Details::FenceCounter::start();
    exec_space.fence();    
    Tpetra::Details::FenceCounter::stop();   
    global_count =Tpetra::Details::FenceCounter::get_count_global(space);   
    instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,0);
    TEST_EQUALITY(instance_count,instance_correct_count);

    // This guy should not get counted, since the counter is stopped
    exec_space.fence();        
    global_count =Tpetra::Details::FenceCounter::get_count_global(space);   
    instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,0);
    TEST_EQUALITY(instance_count,instance_correct_count);

    // Third Second  Stop / Start (should have double the original count)
    Tpetra::Details::FenceCounter::start();
    exec_space.fence();        
    Tpetra::Details::FenceCounter::stop();   
    global_count =Tpetra::Details::FenceCounter::get_count_global(space);   
    instance_count = Tpetra::Details::FenceCounter::get_count_instance(space);   
    TEST_EQUALITY(global_count,0);
    TEST_EQUALITY(instance_count,2*instance_correct_count);      
  }






#ifdef KOKKOS_ENABLE_OPENMP
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, OpenMP_ThreadedSum, LO , GO , Scalar , Node ) {
    // Restrict to OpenMPNode and disable in debug mode (weird things happen w/ GCC 8.3.0 since RCP's
    // are not necessarily thread-safe
    if(typeid(Node)!=typeid(Tpetra::KokkosCompat::KokkosDeviceWrapperNode<Kokkos::OpenMP, Kokkos::HostSpace>) ||
       ::Tpetra::Details::Behavior::debug())
       return;

    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar,LO, GO, Node> MV;
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    RCP<const map_type> map = rcp (new map_type (1000, 0, comm));
    MV x(map, 1);
    x.putScalar(0.0);
    LO N= (LO) map->getLocalNumElements();

    // Vector parallel fill
    // If we
#pragma omp parallel for
    for(LO i=0; i<N; i++) {
      GO global_idx = map->getGlobalElement(i);
      double val = 1.0/global_idx;
      x.sumIntoGlobalValue(global_idx, 0, val, true);
    }
  }
#endif


//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP_BASE( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, basic             , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, large             , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, NonMemberConstructors, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, BadConstLDA       , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, BadConstAA        , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, CopyConst         , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(      Vector, CopyConst         , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(      Vector, Indexing          , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, OrthoDot          , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, CountDot          , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, CountDotNonTrivLDA, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, BadDot            , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, CountNorm1        , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, CountNormInf      , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, Norm2             , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, CopyView          , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, OffsetView        , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, VectorOffsetView  , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, ZeroScaleUpdate   , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(      Vector, ZeroScaleUpdate   , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, ScaleAndAssign    , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, BadMultiply       , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, SingleVecNormalize, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, Multiply          , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, ElementWiseMultiply,LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, ElementWiseMultiplyLg,LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, NonContigView     , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, Describable       , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, Typedefs          , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, ReplaceMap        , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, DeepCopy          , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, DualViewSemantics , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, DualViewCtor      , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, ViewCtor          , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, SubViewSomeZeroRows, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, DimsWithSomeZeroRows, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, DimsWithAllZeroRows, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, Swap, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, DualViewRefcountCheck, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, CopyCounterCheck, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, FenceCounterCheck, LO, GO, SCALAR, NODE )

#ifdef KOKKOS_ENABLE_OPENMP
  // Add special test for OpenMP
  #define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
    UNIT_TEST_GROUP_BASE( SCALAR, LO, GO, NODE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, OpenMP_ThreadedSum, LO, GO, SCALAR, NODE )
#else
  #define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
    UNIT_TEST_GROUP_BASE( SCALAR, LO, GO, NODE )
#endif






  typedef Tpetra::Map<>::local_ordinal_type default_local_ordinal_type;
  typedef Tpetra::Map<>::global_ordinal_type default_global_ordinal_type;

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
#  define TPETRA_MULTIVECTOR_COMPLEX_FLOAT_DOT_TEST( NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, ComplexDotOneColumn, float, default_local_ordinal_type, default_global_ordinal_type, NODE )
#else
#  define TPETRA_MULTIVECTOR_COMPLEX_FLOAT_DOT_TEST( NODE )
#endif // defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
#  define TPETRA_MULTIVECTOR_COMPLEX_DOUBLE_DOT_TEST( NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, ComplexDotOneColumn, double, default_local_ordinal_type, default_global_ordinal_type, NODE )
#else
#  define TPETRA_MULTIVECTOR_COMPLEX_DOUBLE_DOT_TEST( NODE )
#endif // defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)



#define VIEWMODETEST(NODE) \

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_N( TPETRA_MULTIVECTOR_COMPLEX_FLOAT_DOT_TEST )

  TPETRA_INSTANTIATE_N( TPETRA_MULTIVECTOR_COMPLEX_DOUBLE_DOT_TEST )

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

}
