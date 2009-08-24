#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

// FINISH: add test for MultiVector with a node containing zero local entries
// FINISH: add tests for local MultiVectors 

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

  using std::vector;
  using std::sort;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::SerialDenseMatrix;
  using Teuchos::Range1D;
  using Teuchos::Tuple;
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
  using Tpetra::LocallyReplicated;
  using Tpetra::GloballyDistributed;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::MultiVector;
  using Tpetra::global_size_t;
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
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, basic, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const size_t numLocal = 13;
    const size_t numVecs  = 7;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,0,comm) );
    MV mvec(map,numVecs,true);
    TEST_EQUALITY( mvec.getNumVectors(), numVecs );
    TEST_EQUALITY( mvec.getLocalLength(), numLocal );
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
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 13;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,0,comm) );
    TEST_THROW(MV mvec(map,0),  std::invalid_argument);
    if (std::numeric_limits<size_t>::is_signed) {
      TEST_THROW(MV mvec(map,-1), std::invalid_argument);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstLDA, Ordinal, Scalar )
  {
    // numlocal > LDA
    // ergo, the arrayview doesn't contain enough data to specify the entries
    // also, if bounds checking is enabled, check that bad bounds are caught
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const Ordinal indexBase = 0;
    const size_t numLocal = 2;
    const size_t numVecs = 2;
    // multivector has two vectors, each proc having two values per vector
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
    // we need 4 scalars to specify values on each proc
    Array<Scalar> values(4);
#ifdef HAVE_TPETRA_DEBUG
    // too small an ArrayView (less than 4 values) is met with an exception, if debugging is on
    TEST_THROW(MV mvec(map,values(0,3),2,numVecs), std::runtime_error);
    // it could also be too small for the given LDA: 
    TEST_THROW(MV mvec(map,values(),2+1,numVecs), std::runtime_error);
    // too small for number of entries in a Vector
    TEST_THROW(V   vec(map,values(0,1)), std::runtime_error);
#endif
    // LDA < numLocal throws an exception anytime
    TEST_THROW(MV mvec(map,values(0,4),1,numVecs), std::runtime_error);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, NonContigView, Ordinal, Scalar )
  {
    if (ScalarTraits<Scalar>::isOrdinal) return;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node> V;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Mag tol = errorTolSlack * ScalarTraits<Mag>::eps();
    const Mag M0  = ScalarTraits<Mag>::zero();
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 13;
    const size_t numVecs = 7;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
    //
    // we will create a non-contig subview of the vector; un-viewed vectors should not be changed
    Tuple<size_t,4> inView1 = tuple<size_t>(1,4,3,2);
    Tuple<size_t,3> exView1 = tuple<size_t>(0,5,6);
    Tuple<size_t,4> inView2 = tuple<size_t>(6,0,4,3);
    Tuple<size_t,4> exView2 = tuple<size_t>(1,2,5,7);
    const size_t numView = 4;
    TEST_FOR_EXCEPTION(numView != as<size_t>(inView1.size()), std::logic_error, "Someone ruined a test invariant.");
    TEST_FOR_EXCEPTION(numView != as<size_t>(inView1.size()), std::logic_error, "Someone ruined a test invariant.");
    TEST_FOR_EXCEPTION(numView != as<size_t>(inView2.size()), std::logic_error, "Someone ruined a test invariant.");
    {
      // test dot, all norms, randomize
      MV mvOrig1(map,numVecs), mvOrig2(map,numVecs+1), mvWeights(map,numVecs);
      mvWeights.randomize();
      RCP<const MV> mvW1 = mvWeights.subView(tuple<size_t>(0));
      RCP<const MV> mvSubWeights = mvWeights.subView(inView1);
      mvOrig1.randomize();
      mvOrig2.randomize();
      // 
      Array<Mag> nOrig2(numVecs), nOrig1(numVecs), nOrigI(numVecs), nOrigW(numVecs), nOrigW1(numVecs);
      Array<Scalar> meansOrig(numVecs), dotsOrig(numView);
      mvOrig1.norm1(nOrig1());
      mvOrig1.norm2(nOrig2());
      mvOrig1.normInf(nOrigI());
      mvOrig1.normWeighted(mvWeights,nOrigW());
      mvOrig1.normWeighted(*mvW1,nOrigW1());
      mvOrig1.meanValue(meansOrig());
      for (size_t j=0; j < numView; ++j) {
        RCP<const V> v1 = mvOrig1.getVector(inView1[j]),
          v2 = mvOrig2.getVector(inView2[j]);
        dotsOrig[j] = v1->dot(*v2);
      }
      // create the views, compute and test
      RCP<      MV> mvView1 = mvOrig1.subViewNonConst(inView1);
      RCP<const MV> mvView2 = mvOrig2.subView(inView2);
      Array<Mag> nView2(numView), nView1(numView), nViewI(numView), nViewW(numView), nViewW1(numView);
      Array<Scalar> meansView(numView), dotsView(numView);
      mvView1->norm1(nView1());
      mvView1->norm2(nView2());
      mvView1->normInf(nViewI());
      mvView1->normWeighted(*mvSubWeights,nViewW());
      mvView1->normWeighted(*mvW1,nViewW1());
      mvView1->meanValue(meansView());
      mvView1->dot( *mvView2, dotsView() );
      for (size_t j=0; j < numView; ++j) {
        TEST_FLOATING_EQUALITY(nOrig1[inView1[j]],  nView1[j],  tol);
        TEST_FLOATING_EQUALITY(nOrig2[inView1[j]],  nView2[j],  tol);
        TEST_FLOATING_EQUALITY(nOrigI[inView1[j]],  nViewI[j],  tol);
        TEST_FLOATING_EQUALITY(nOrigW[inView1[j]],  nViewW[j],  tol);
        TEST_FLOATING_EQUALITY(nOrigW1[inView1[j]], nViewW1[j], tol);
        TEST_FLOATING_EQUALITY(meansOrig[inView1[j]], meansView[j], tol);
        TEST_FLOATING_EQUALITY(dotsOrig[j], dotsView[j], tol);
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
      mvOrigA.randomize();
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
        (*doubleViewA) = (*doubleViewB) = (*doubleViewC);
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
        TEST_FOR_EXCEPTION(inView1.size() != 4, std::logic_error, "Someone ruined a test invariant.");
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
        //   B.recip(C)           ->   B == 1, A == 0,            reciprocal(mv)
        dvB->reciprocal(*dvC);
        //   C = 2*B              ->   A == 0, B == 1, C == 2
        dvC->scale(as<Mag>(2),*dvB);
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
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, LabeledObject, LO, GO, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // create Map
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,as<size_t>(3),0,comm) );
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadMultiply, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const Scalar S1 = ScalarTraits<Scalar>::one(),
                 S0 = ScalarTraits<Scalar>::zero();
    // case 1: C(local) = A^X(local) * B^X(local)  : four of these
    {
      // create local Maps
      RCP<Map<Ordinal,Ordinal,Node> > map3l = rcp( new Map<Ordinal,Ordinal,Node>(as<global_size_t>(3),as<Ordinal>(0),comm,LocallyReplicated) ),
                                      map2l = rcp( new Map<Ordinal,Ordinal,Node>(as<global_size_t>(2),as<Ordinal>(0),comm,LocallyReplicated) );
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
      RCP<Map<Ordinal,Ordinal,Node> > map3n = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(3),0,comm) ),
                                      map2n = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(2),0,comm) ),
                                      map2l = rcp( new Map<Ordinal,Ordinal,Node>(as<global_size_t>(2),as<Ordinal>(0),comm,LocallyReplicated) ),
                                      map3l = rcp( new Map<Ordinal,Ordinal,Node>(as<global_size_t>(3),as<Ordinal>(0),comm,LocallyReplicated) );
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
      RCP<Map<Ordinal,Ordinal,Node> > map3n = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(3),0,comm) ),
                                      map2n = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(2),0,comm) ),
                                      map2l = rcp( new Map<Ordinal,Ordinal,Node>(as<global_size_t>(2),as<Ordinal>(0),comm,LocallyReplicated) ),
                                      map3l = rcp( new Map<Ordinal,Ordinal,Node>(as<global_size_t>(3),as<Ordinal>(0),comm,LocallyReplicated) );
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
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, Multiply, Ordinal, Scalar )
  {
    using Teuchos::View;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const Ordinal indexBase = 0;
    // create a Map
    RCP<Map<Ordinal,Ordinal,Node> > map3n = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(3),indexBase,comm) ),
                                    map2n = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(2),indexBase,comm) ),
                                    lmap3 = rcp( new Map<Ordinal,Ordinal,Node>(as<global_size_t>(3),indexBase,comm,LocallyReplicated) ),
                                    lmap2 = rcp( new Map<Ordinal,Ordinal,Node>(as<global_size_t>(2),indexBase,comm,LocallyReplicated) );
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
      MV mv3nx2(map3n,2),
         mv3nx3(map3n,3),
         // locals
         mv2x3(lmap2,3);
      // fill multivectors with ones
      mv2x3.putScalar(S1);
      // fill expected answers Array
      ArrayRCP<const Scalar> tmpView;
      Teuchos::Array<Scalar> check2(9,2), check3(6,3);
      // test
      mv3nx3.putScalar(S1); mv3nx2.putScalar(S1);
      mv3nx3.multiply(NO_TRANS,  NO_TRANS,S1,mv3nx2,mv2x3,S0);
      tmpView = mv3nx3.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check2,M0);
      mv3nx3.putScalar(S1); mv3nx2.putScalar(S1);
      mv3nx2.multiply(NO_TRANS,CONJ_TRANS,S1,mv3nx3,mv2x3,S0);
      tmpView = mv3nx2.get1dView(); TEST_COMPARE_FLOATING_ARRAYS(tmpView,check3,M0);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadConstAA, Ordinal, Scalar )
  {
    // constructor takes ArrayView<ArrayView<Scalar> A, NumVectors
    // A.size() == NumVectors
    // A[i].size() >= MyLength
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 2;
    // multivector has two vectors, each proc having two values per vector
    RCP<Map<Ordinal,Ordinal,Node> > map2 = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(numLocal  ),indexBase,comm) ),
                                    map3 = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(numLocal+1),indexBase,comm) );
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
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadDot, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<Map<Ordinal,Ordinal,Node> > map1 = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(1),indexBase,comm) ),
                                    map2 = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(2),indexBase,comm) );
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
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar S0 = ScalarTraits<Scalar>::zero();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
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
    TEST_EQUALITY_CONST( mvec1.getVector(0)->dot(*mvec2.getVector(0)), S0);
    mvec1.norm1(norms1());
    mvec2.norm1(norms2());
    std::fill(ans.begin(), ans.end(), M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    // replace local entries s.t. 
    // mvec1 = [1 1]  and  mvec2 = [0 0]
    //         [0 0]               [1 1]
    // still numerically orthogonal even in finite arithmetic. norms are numImages.
    for (size_t j=0; j < numVectors; ++j) {
      mvec1.replaceLocalValue(0,j,ScalarTraits<Scalar>::one());
      mvec2.replaceGlobalValue(map->getGlobalIndex(1),j,ScalarTraits<Scalar>::one());
    }
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots2,zeros,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
    TEST_EQUALITY_CONST( mvec1.getVector(0)->dot(*mvec2.getVector(0)), S0);
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
      mvec2.sumIntoGlobalValue(map->getGlobalIndex(0),j,-ScalarTraits<Scalar>::one());
    }
    mvec1.dot(mvec2,dots1());
    mvec2.dot(mvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots2,zeros,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
    TEST_EQUALITY_CONST( mvec1.getVector(0)->dot(*mvec2.getVector(0)), S0);
    mvec1.norm1(norms1());
    mvec2.norm1(norms2());
    std::fill(ans.begin(), ans.end(), as<Mag>(2*numImages));
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms2,ans,M0);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CopyView, Ordinal, Scalar )
  {
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar S0 = ScalarTraits<Scalar>::zero();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 7;
    const size_t numVectors = 13;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
    MV A(map,numVectors,false);
    {
      A.randomize();
      TEST_FOR_EXCEPT(numVectors != 13);
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
      A.norm2(A_bef());
      // get view and its norms
      RCP<MV> Av = A.subViewNonConst(inds1);
      Av->norm2(Av_bef());
      // get copy and its norms
      RCP<MV> Ac = A.subCopy(inds1);
      Ac->norm2(Ac_bef());
      // set view to zero
      Av->putScalar(ScalarTraits<Scalar>::zero());
      // get norms of view
      Av->norm2(Av_aft());
      // free the view, copying data back to A
      Av = Teuchos::null;
      // get norms of A and copy
      Ac->norm2(Ac_aft());
      A.norm2(A_aft());
      // norms of copy and view before should match norms of A
      for (size_t i=0; i < as<size_t>(inds1.size()); ++i) {
        TEST_EQUALITY( A_bef[inds1.lbound()+i], Ac_bef[i] );
      }
      TEST_COMPARE_FLOATING_ARRAYS(Ac_bef,Av_bef,M0);
      // norms of copy (before and after) should match
      TEST_COMPARE_FLOATING_ARRAYS(Ac_bef,Ac_aft,M0);
      // norms of view after should be zero, as should corresponding A norms
      for (size_t i=0; i < as<size_t>(inds1.size()); ++i) {
        TEST_EQUALITY_CONST( Av_aft[i], M0 );
        TEST_EQUALITY_CONST( A_aft[inds1.lbound()+i], M0 );
      }
    }
    {
      A.randomize();
      TEST_FOR_EXCEPT(numVectors != 13);
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
      A.norm2(A_bef());
      // get view and its norms
      RCP<MV> Av = A.subViewNonConst(inds);
      Av->norm2(Av_bef());
      // get copy and its norms
      RCP<MV> Ac = A.subCopy(inds);
      Ac->norm2(Ac_bef());
      // set view to zero
      Av->putScalar(ScalarTraits<Scalar>::zero());
      // get norms of view
      Av->norm2(Av_aft());
      // free the view, copying data back to A
      Av = Teuchos::null;
      // get norms of A and copy
      Ac->norm2(Ac_aft());
      A.norm2(A_aft());
      // norms of copy and view before should match norms of A
      for (size_t i=0; i < as<size_t>(inds.size()); ++i) {
        TEST_EQUALITY( A_bef[inds[i]], Ac_bef[i] );
      }
      TEST_COMPARE_FLOATING_ARRAYS(Ac_bef,Av_bef,M0);
      // norms of copy (before and after) should match
      TEST_COMPARE_FLOATING_ARRAYS(Ac_bef,Ac_aft,M0);
      // norms of view after should be zero, as should corresponding A norms
      for (size_t i=0; i < as<size_t>(inds.size()); ++i) {
        TEST_EQUALITY_CONST( Av_aft[i], M0 );
        TEST_EQUALITY_CONST( A_aft[inds[i]], M0 );
      }
    }
    {
      A.randomize();
      Array<Mag> Anorms(numVectors);
      A.norm2(Anorms());
      TEST_FOR_EXCEPT(numVectors != 13);
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
          sub2->norm2(subnorms());
          TEST_COMPARE_FLOATING_ARRAYS(Anorms(6,3),subnorms(),M0);
        }
      }
    }
    {
      A.randomize();
      {
        // check that 1dView and 1dCopy have the same values
        ArrayRCP<const Scalar> view;
        Array<Scalar> copy(numLocal*numVectors);
        view = A.get1dView();
        A.get1dCopy(copy(),numLocal);
        TEST_COMPARE_FLOATING_ARRAYS(view,copy,M0);
      }
      {
        // check that 1dView and 1dCopy have the same values
        ArrayRCP<Scalar> view;
        Array<Scalar> copy(numLocal*numVectors);
        view = A.get1dViewNonConst();
        A.get1dCopy(copy(),numLocal);
        TEST_COMPARE_FLOATING_ARRAYS(view,copy,M0);
        // clear view, ensure that A is zero
        std::fill(view.begin(), view.end(), S0);
        view = Teuchos::null;
        Array<Mag> norms(numVectors), zeros(numVectors,M0);
        A.norm2(norms());
        TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
      }
      A.randomize();
      {
        // check that 1dView and 1dCopy have the same values
        ArrayRCP<ArrayRCP<const Scalar> > views;
        Array<Scalar> copyspace(numLocal*numVectors);
        Array<ArrayView<Scalar> > copies(numVectors);
        for (size_t j=0; j < numVectors; ++j) {
          copies[j] = copyspace(numLocal*j,numLocal);
        }
        views = A.get2dView();
        A.get2dCopy(copies());
        for (size_t j=0; j < numVectors; ++j) {
          TEST_COMPARE_FLOATING_ARRAYS(views[j],copies[j],M0);
        }
      }
      {
        // check that 1dView and 1dCopy have the same values
        ArrayRCP<ArrayRCP<Scalar> > views;
        Array<Scalar> copyspace(numLocal*numVectors);
        Array<ArrayView<Scalar> > copies(numVectors);
        for (size_t j=0; j < numVectors; ++j) {
          copies[j] = copyspace(numLocal*j,numLocal);
        }
        views = A.get2dViewNonConst();
        A.get2dCopy(copies());
        for (size_t j=0; j < numVectors; ++j) {
          TEST_COMPARE_FLOATING_ARRAYS(views[j],copies[j],M0);
        }
        // clear view, ensure that A is zero
        for (size_t j=0; j < numVectors; ++j) {
          std::fill(views[j].begin(), views[j].end(), S0);
        }
        views = Teuchos::null;
        Array<Mag> norms(numVectors), zeros(numVectors,M0);
        A.norm2(norms());
        TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
      }
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, ZeroScaleUpdate, Ordinal, Scalar )
  {
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 2;
    const size_t numVectors = 2;
    const size_t LDA = 2;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
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
      MV C(map,numVectors);
      C.randomize();
      C.update(as<Scalar>(-1),B,as<Scalar>(2),A,as<Scalar>(0));
      C.norm2(norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }
    //   set C random
    //   scale it ex-situ
    //   check that it equals B: subtraction in situ
    {
      MV C(map,numVectors);
      C.scale(as<Scalar>(2),A);
      C.update(as<Scalar>(1),B,as<Scalar>(-1));
      C.norm2(norms);
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, ScaleAndAssign, Ordinal, Scalar )
  {
    if (ScalarTraits<Scalar>::isOrdinal) return;
    Teuchos::ScalarTraits<Scalar>::seedrandom(0);   // consistent seed
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 23;
    const size_t numVectors = 11;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
    // Use random multivector A
    // Set B = A * 2 manually.
    // Therefore, if C = 2*A, then C == B
    // If C = A and C *= 2, then C == B
    // This test operator= and all of our scale ops
    // We'll do Vector and MultiVector variations
    // Also, ensure that other vectors aren't changed
    MV A(map,numVectors,false),
       B(map,numVectors,false);
    A.randomize();
    Array<Mag> Anrms(numVectors);
    A.norm2(Anrms());
    // set B = A * 2, using different techniques
    // * get vector, Vector::operator=
    // * get 1-vector subview(Range1D), MultiVector::operator=
    // * get 1-vector subview(ArrayView), MultiVector::operator=
    // * get data view, assign
    TEST_FOR_EXCEPT(numVectors < 4);
    for (size_t j = 0; j < numVectors; ++j) {
      // assign j-th vector of B to 2 * j-th vector of A
      switch (j % 4) {
        case 0:
          {
            RCP<V> bj = B.getVectorNonConst(j);
            RCP<const V> aj = A.getVector(j);
            (*bj) = (*aj);
            ArrayRCP<Scalar> bjview = bj->get1dViewNonConst();
            for (size_t i=0; i < numLocal; ++i) {
              bjview[i] *= as<Scalar>(2);
            }
          }
          break;
        case 1:
          {
            RCP<MV>       bj = B.subViewNonConst(Range1D(j,j));
            RCP<const MV> aj = A.subView(Range1D(j,j));
            (*bj) = (*aj);
            ArrayRCP<Scalar> bjview = bj->get1dViewNonConst();
            for (size_t i=0; i < numLocal; ++i) {
              bjview[i] *= as<Scalar>(2);
            }
          }
          break;
        case 2:
          {
            RCP<MV> bj = B.subViewNonConst(tuple<size_t>(j));
            RCP<const MV> aj = A.subView(tuple<size_t>(j));
            (*bj) = (*aj);
            ArrayRCP<Scalar> bjview = bj->get1dViewNonConst();
            for (size_t i=0; i < numLocal; ++i) {
              bjview[i] *= as<Scalar>(2);
            }
          }
          break;
        case 3:
          {
            ArrayRCP<Scalar>       bjview = B.getDataNonConst(j);
            ArrayRCP<const Scalar> ajview = A.getData(j);
            for (size_t i=0; i < numLocal; ++i) {
              bjview[i] = as<Scalar>(2) * ajview[i];
            }
          }
          break;
      }
    }
    // check that A wasn't modified
    {
      Array<Mag> Anrms_aft(numVectors);
      A.norm2(Anrms_aft());
      TEST_COMPARE_FLOATING_ARRAYS(Anrms(),Anrms_aft(),M0);
    }
    // check that C.Scale(A,2.0) == B
    {
      MV C(map,numVectors,false);
      C.scale(as<Scalar>(2), A);
      C.update(-1.0,B,1.0);
      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm2(Cnorms());
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms(),zeros,M0);
    }
    // check that C=A, C.Scale(2.0) == B
    {
      MV C(map,numVectors,false);
      C = A;
      C.scale(as<Scalar>(2));
      C.update(-1.0,B,1.0);
      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm2(Cnorms());
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms(),zeros,M0);
    }
    // check that C=A, C.Scale(tuple(2)) == B
    {
      MV C(map,numVectors,false);
      C = A;
      Array<Scalar> twos(numVectors,as<Scalar>(2));
      C.scale(twos());
      C.update(-1.0,B,1.0);
      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm2(Cnorms());
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms(),zeros,M0);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, ZeroScaleUpdate, Ordinal, Scalar )
  {
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(2),0,comm) );
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
      V C(map);
      C.randomize();
      C.update(as<Scalar>(-1),B,as<Scalar>(2),A,as<Scalar>(0));
      norm = C.norm2(); C.norm2(norms());
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
      norm = C.norm2(); C.norm2(norms());
      TEST_EQUALITY(norm,M0);
      TEST_EQUALITY(norm,norms[0]);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CopyConst, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 13;
    const size_t numVectors = 7;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
    {
      // create random MV
      MV mvorig(map,numVectors);
      mvorig.randomize();
      // create non-const subview, test copy constructor
      TEST_FOR_EXCEPT(numVectors != 7);
      Tuple<size_t,3> inds = tuple<size_t>(1,3,5);
      RCP<MV> mvview = mvorig.subViewNonConst(inds);
      Array<Mag> norig(numVectors), nsub(inds.size()), ncopy(inds.size());
      mvorig.normInf(norig());
      for (size_t j=0; j < as<size_t>(inds.size()); ++j) {
        nsub[j] = norig[inds[j]];
      }
      MV mvcopy(*mvview);
      mvcopy.normInf(ncopy());
      TEST_COMPARE_FLOATING_ARRAYS(ncopy,nsub,M0);
      // reset both the view and the copy of the view, ensure that they are independent
      Teuchos::Array<Mag> nsub_aft(inds.size()), ones(inds.size(),as<Mag>(1));
      Teuchos::Array<Mag> ncopy_aft(inds.size()), twos(inds.size(),as<Mag>(2));
      mvview->putScalar(as<Scalar>(1));
      mvcopy.putScalar(as<Scalar>(2));
      mvview->normInf(nsub_aft());
      mvcopy.normInf(ncopy_aft());
      TEST_COMPARE_FLOATING_ARRAYS(nsub_aft,ones,M0);
      TEST_COMPARE_FLOATING_ARRAYS(ncopy_aft,twos,M0);
    }
    {
      // create random MV
      MV morig(map,numVectors);
      morig.randomize();
      // test copy constructor with 
      // copy it
      MV mcopy1(morig), mcopy2(morig);
      // verify that all three have identical values
      Array<Mag> norig(numVectors), ncopy1(numVectors), ncopy2(numVectors);
      morig.normInf(norig);
      mcopy1.normInf(ncopy1);
      mcopy2.normInf(ncopy2);
      TEST_COMPARE_FLOATING_ARRAYS(norig,ncopy1,M0);
      TEST_COMPARE_FLOATING_ARRAYS(norig,ncopy2,M0);
      // modify all three
      morig.putScalar(as<Scalar>(0));
      mcopy1.putScalar(as<Scalar>(1));
      mcopy2.putScalar(as<Scalar>(2));
      // compute norms, check
      Array<Mag> zeros(numVectors,as<Mag>(0)), ones(numVectors,as<Mag>(1)), twos(numVectors,as<Mag>(2));
      morig.normInf(norig);
      mcopy1.normInf(ncopy1);
      mcopy2.normInf(ncopy2);
      TEST_COMPARE_FLOATING_ARRAYS(norig,zeros,M0);
      TEST_COMPARE_FLOATING_ARRAYS(ncopy1,ones,M0);
      TEST_COMPARE_FLOATING_ARRAYS(ncopy2,twos,M0);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, CopyConst, Ordinal, Scalar )
  {
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node>       V;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(2),0,comm) );
    // create random MV
    V morig(map);
    morig.randomize();
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
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(100),0,comm) );
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
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Magnitude M1  = ScalarTraits<Magnitude>::one();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 10;
    const size_t numVectors = 6;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
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
    TEST_COMPARE_FLOATING_ARRAYS(norms,ones,ScalarTraits<Magnitude>::eps()*as<Magnitude>(10.));
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountDot, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
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
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Magnitude M0 = ScalarTraits<Magnitude>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    const size_t LDA = 3;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
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
    TEST_COMPARE_FLOATING_ARRAYS(dots1,dots2,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,answer,M0);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountNorm1, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
    Array<Scalar> values(6);
    // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    //                               [0 1 2]
    // norm1(values) = [0 2 4]
    // over all procs, this is [0 2*nprocs 4*nprocs]
    // mean is [0 1 2]
    values[0] = as<Scalar>(0);
    values[1] = as<Scalar>(0);
    values[2] = as<Scalar>(1);
    values[3] = as<Scalar>(1);
    values[4] = as<Scalar>(2);
    values[5] = as<Scalar>(2);
    MV mvec(map,values(),2,numVectors);
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
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CountNormInf, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 2;
    const size_t numVectors = 3;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
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
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, Norm2, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const MT M0 = ScalarTraits<MT>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 13;
    const size_t numVectors = 7;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
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
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, NormWeighted, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Mag tol = errorTolSlack * ScalarTraits<Mag>::eps();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const Ordinal indexBase = 0;
    const size_t numLocal = 13;
    const size_t numVectors = 7;
    RCP<Map<Ordinal,Ordinal,Node> > map = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,numLocal,indexBase,comm) );
    MV    mvec(map,numVectors),
       weights(map,numVectors),
       weight1(map,1);
    // randomize the multivector
    mvec.randomize();
    // set the weights
    Array<Scalar> wvec(numVectors);
    Scalar w1 = ScalarTraits<Scalar>::random();
    for (size_t j=0; j < numVectors; ++j) {
      wvec[j] = ScalarTraits<Scalar>::random();
    }
    weights.putScalar(ScalarTraits<Scalar>::one());
    weights.scale(wvec());
    weight1.putScalar(w1);
    // take norms
    Array<Mag> normsW(numVectors), normsW1(numVectors);
    Array<Scalar> dots(numVectors);
    mvec.dot(mvec,dots());
    mvec.normWeighted(weights,normsW());
    mvec.normWeighted(weight1,normsW1());
    for (size_t j=0; j < numVectors; ++j) {
      Mag ww = ScalarTraits<Scalar>::real( ScalarTraits<Scalar>::conjugate(wvec[j]) * wvec[j] );
      Mag expnorm = ScalarTraits<Mag>::squareroot( 
                      ScalarTraits<Scalar>::real(dots[j]) / (as<Mag>(numImages * numLocal) * ww)
                    );
      Mag ww1 = ScalarTraits<Scalar>::real( ScalarTraits<Scalar>::conjugate(w1) * w1 );
      Mag expnorm1 = ScalarTraits<Mag>::squareroot(
                       ScalarTraits<Scalar>::real(dots[j]) / (as<Mag>(numImages * numLocal) * ww1)
                     );
      TEST_FLOATING_EQUALITY( expnorm, normsW[j], tol );
      TEST_FLOATING_EQUALITY( expnorm1, normsW1[j], tol );
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, BadCombinations, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // create a Map
    const Ordinal indexBase = 0;
    const Scalar rnd = ScalarTraits<Scalar>::random();
    // two maps: one has two entires per node, the other disagrees on node 0
    RCP<Map<Ordinal,Ordinal,Node> > map1 = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(2),indexBase,comm) ),
                                    map2 = rcp( new Map<Ordinal,Ordinal,Node>(INVALID,as<size_t>(myImageID == 0 ? 1 : 2),indexBase,comm) );
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
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

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
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, NormWeighted      , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CopyView          , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, ZeroScaleUpdate   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(      Vector, ZeroScaleUpdate   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, ScaleAndAssign    , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadCombinations   , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, BadMultiply       , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, SingleVecNormalize, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, Multiply          , ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, NonContigView     , ORDINAL, SCALAR ) \
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
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)
     UNIT_TEST_GROUP_ORDINAL(int)
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt; UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#    endif

#endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
