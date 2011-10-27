#include "MueLu_TestHelpers2.hpp"

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"

// Belos / Xpetra-MueLu adapters
#include "BelosXpetraAdapter.hpp"
#include "BelosMueLuAdapter.hpp"

namespace MueLuTests {
   
  //
  // Helpers function to build tests
  //

  // Test Belos adapters for the couple <MV,OP> 
  // TODO: add a bunch of 'const' on prototype
  template <class Scalar, class MV, class OP>
  int BelosAdaptersTest(RCP<OP> & belosOp, RCP<OP> & belosPrec, RCP<MV> & X, RCP<MV> & B, Teuchos::FancyOStream & out, bool & success) {
    RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, B));
    belosProblem->setLeftPrec(belosPrec);
    
    bool set = belosProblem->setProblem();
    TEST_EQUALITY(set, true);
    
    // Belos parameter list
    Teuchos::ParameterList belosList;
    belosList.set("Maximum Iterations",    10);   // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", 1e-7); // Relative convergence tolerance requested
    
    // Create an iterative solver manager.
    RCP<Belos::SolverManager<Scalar, MV, OP> > belosSolver = rcp(new Belos::BlockCGSolMgr<double,MV,OP>(belosProblem, rcp(&belosList,false)));
    
    // Perform solve
    Belos::ReturnType ret = belosSolver->solve();
    TEST_EQUALITY(ret, Belos::Converged);
    
    // Return number of iterations
    return belosSolver->getNumIters();
  }

  //
  // Helpers function to verify results
  //

  // Singleton for norm comparisons across tests
  template <class Scalar>
  bool BelosAdaptersTestResultsNorm(typename Teuchos::ScalarTraits<Scalar>::magnitudeType r) {
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType ref = -1;
    if (ref == -1) {
      //std::cout << "BelosAdaptersTestResults(): Set reference results" << std::endl;
      ref = r;
      return true;
    }
    //std::cout << "BelosAdaptersTestResults(): Compare" << std::endl;

    if (r != ref)
      std::cout << "ref  norm = " << ref << std::endl
                << "curr norm = " << r   << std::endl;

    return (r == ref);
  }

  // Test results
  template <class Scalar, class MV>
  bool BelosAdaptersTestResults(int numIters, RCP<MV> & X, Teuchos::FancyOStream & out, bool & success) {

    // Check numIters
    switch (TestHelpers::Parameters::getDefaultComm()->getSize()) { 
    case 0: TEST_EQUALITY(numIters, 5); break;
    case 4: 
      // Epetra TEST_EQUALITY(numIters, 6);
      // Tpetra TEST_EQUALITY(numIters, 7); 
      break;
    default:;
    }
    
    // Compute norm of X (using MV traits)
    typedef Belos::MultiVecTraits<Scalar, MV> MVT;
    std::vector<Scalar> norms(1);
    MVT::MvNorm(*X, norms);
    
    // Test norm equality across the unit tests
    return MueLuTests::BelosAdaptersTestResultsNorm<Scalar>(norms[0]);
  }

  //
  // Tests
  //

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

  // TEST:
  // - OP: Xpetra::Operator
  // - MV: Xpetra::MultiVector
  TEUCHOS_UNIT_TEST(BelosAdapters, XpetraOp_XpetraMV) {
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    RCP<TestProblem<SC,LO,GO,NO,LMO> > p = TestHelpers::getTestProblem<SC,LO,GO,NO,LMO>(lib);

    typedef Xpetra::MultiVector<SC> MV;
    typedef Belos::OperatorT<MV>    OP;
    
    // Construct a Belos LinearProblem object
    RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>    (p->GetA()));
    RCP<OP> belosPrec = rcp(new Belos::MueLuPrecOp<SC, LO, GO, NO, LMO>(p->GetH()));
    
    // Run Belos
    RCP<MultiVector> X = p->GetNewX0();
    int numIters = MueLuTests::BelosAdaptersTest<SC, MV, OP>(belosOp, belosPrec, X, p->GetRHS(), out, success);

    // Tests
    TEST_EQUALITY(MueLuTests::BelosAdaptersTestResults<Scalar>(numIters, X, out, success), true);
  }

#ifdef HAVE_MUELU_EPETRA
  // TEST:
  // - OP: Xpetra::Operator
  // - MV: Epetra::MultiVector
  TEUCHOS_UNIT_TEST(BelosAdapters, XpetraOp_EpetraMV) {
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();
    if (lib == Xpetra::UseEpetra) {  // Epetra specific test: run only once.

      RCP<TestProblem<SC,LO,GO,NO,LMO> > p = TestHelpers::getTestProblem<SC,LO,GO,NO,LMO>(lib);
      
      typedef Epetra_MultiVector   MV;
      typedef Belos::OperatorT<MV> OP;
      
      // Construct a Belos LinearProblem object
      RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>    (p->GetA()));
      RCP<OP> belosPrec = rcp(new Belos::MueLuPrecOp<SC, LO, GO, NO, LMO>(p->GetH()));
      
      // X, B
      RCP<MV> X = Utils::MV2NonConstEpetraMV(p->GetNewX0());
      RCP<MV> B = Utils::MV2NonConstEpetraMV(p->GetRHS());

      // Run Belos      
      int numIters = MueLuTests::BelosAdaptersTest<SC, MV, OP>(belosOp, belosPrec, X, B, out, success);
      
      // Tests
      TEST_EQUALITY(MueLuTests::BelosAdaptersTestResults<Scalar>(numIters, X, out, success), true);
    }
  }

  // TEST:
  // - OP: Belos::MultiVec<double>
  // - MV: Belos::Operator<double>
  TEUCHOS_UNIT_TEST(BelosAdapters, BelosMultiVec_BelosOperator) {
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();
    if (lib == Xpetra::UseEpetra) {  // Epetra specific test: run only once.

      RCP<TestProblem<SC,LO,GO,NO,LMO> > p = TestHelpers::getTestProblem<SC,LO,GO,NO,LMO>(lib);
      
      typedef Belos::MultiVec<double> MV;
      typedef Belos::Operator<double> OP;

      // Construct a Belos LinearProblem object
      RCP<Epetra_CrsMatrix> A = Utils::Op2NonConstEpetraCrs(p->GetA());
      RCP<OP> belosOp   = rcp(new Belos::EpetraOp(A));
      RCP<OP> belosPrec = rcp(new Belos::MueLuPrecOp<SC, LO, GO, NO, LMO>(p->GetH()));

      // X, B
      RCP<Epetra_MultiVector> eX = Utils::MV2NonConstEpetraMV(p->GetNewX0());
      RCP<Epetra_MultiVector> eB = Utils::MV2NonConstEpetraMV(p->GetRHS());
      RCP<MV> X = rcp(new Belos::EpetraMultiVec(*eX));
      RCP<MV> B = rcp(new Belos::EpetraMultiVec(*eB));

      // Run Belos      
      int numIters = MueLuTests::BelosAdaptersTest<SC, MV, OP>(belosOp, belosPrec, X, B, out, success);
      
      // Tests
      TEST_EQUALITY(MueLuTests::BelosAdaptersTestResults<Scalar>(numIters, X, out, success), true);

      // TODO: this do not work. Is it a bug?
      //  double norm;
      //  eX->Norm2(&norm);
    }
  }
#endif

#ifdef HAVE_MUELU_TPETRA
  // TEST:
  // - OP: Xpetra::Operator
  // - MV: Tpetra::MultiVector
  TEUCHOS_UNIT_TEST(BelosAdapters, XpetraOp_TpetraMV) {
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();
    if (lib == Xpetra::UseTpetra) {  // Tpetra specific test: run only once.

      RCP<TestProblem<SC,LO,GO,NO,LMO> > p = TestHelpers::getTestProblem<SC,LO,GO,NO,LMO>(lib);
      
      typedef Tpetra::MultiVector<SC> MV;
      typedef Belos::OperatorT<MV>    OP;
      
      // Construct a Belos LinearProblem object
      RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>    (p->GetA()));
      RCP<OP> belosPrec = rcp(new Belos::MueLuPrecOp<SC, LO, GO, NO, LMO>(p->GetH()));
      
      //X, B
      RCP<MV> X = Utils::MV2NonConstTpetraMV(p->GetNewX0());
      RCP<MV> B = Utils::MV2NonConstTpetraMV(p->GetRHS());

      // Run Belos      
      int numIters = MueLuTests::BelosAdaptersTest<SC, MV, OP>(belosOp, belosPrec, X, B, out, success);

      // Tests
      TEST_EQUALITY(MueLuTests::BelosAdaptersTestResults<Scalar>(numIters, X, out, success), true);
    }
  }
#endif

  // TODO : Tpetra and Epetra are not giving the same results on the test problem. Must be changed for this test.
#ifdef MUELU_DISABLED

  // TEST: Compare Epetra and Tpetra results for:
  // - OP: Xpetra::Operator
  // - MV: Xpetra::MultiVector
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_EPETRA)
  TEUCHOS_UNIT_TEST(BelosAdapters, XpetraOp_XpetraMV_EpetraVsTpetra) {
    // Unit test executable is called twice by ctest (for --linAlgebra=0 and 1). But there is no need to run this test twice. So it only runs for lib == Tpetra.
    if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
      
      // Test for Tpetra will be done by XpetraOp_XpetraMV.
      // We only need to run tests for Epetra to force the result comparisons.
      //
      // TODO
      //

    }
  }
#endif

#endif // MUELU_DISABLED

} // namespace MueLuTests

//TODO: norm test can be factorized, using Belos Adapter Norm function.
