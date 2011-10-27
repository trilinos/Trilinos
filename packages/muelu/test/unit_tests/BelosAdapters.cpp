#include "MueLu_TestHelpers2.hpp"

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"

// Belos / MueLu adapters
#include "BelosMueLuAdapter.hpp"

namespace MueLuTests {
   
  //
  // Test helpers
  //

  // Singleton for results comparisons across tests
  template <class Scalar>
  bool BelosAdaptersTestResults(typename Teuchos::ScalarTraits<Scalar>::magnitudeType r) {
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType ref = -1;
    if (ref == -1) {
      //std::cout << "BelosAdaptersTestResults(): Set reference results" << std::endl;
      ref = r;
      return true;
    }
    //std::cout << "BelosAdaptersTestResults(): Compare" << std::endl;
    return (r == ref);
  }

  // Test Belos adapters for the couple <MV,OP> 
  // TODO: add a bunch of 'const' on prototype
  template <class Scalar, class MV, class OP>
  void BelosAdaptersTest(RCP<OP> & belosOp, RCP<OP> & belosPrec, RCP<MV> & X, RCP<MV> & B, Teuchos::FancyOStream & out, bool & success) {
    RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, B));
    belosProblem->setLeftPrec(belosPrec);
    
    bool set = belosProblem->setProblem();
    TEST_EQUALITY(set, true);
    
    // Belos parameter list
    Teuchos::ParameterList belosList;
    belosList.set("Maximum Iterations",    10);   // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", 1e-7); // Relative convergence tolerance requested
    
    // Create an iterative solver manager.
    RCP<Belos::SolverManager<Scalar, MV, OP> > belosSolver = rcp( new Belos::BlockCGSolMgr<double,MV,OP>(belosProblem, rcp(&belosList,false)) );
    
    // Perform solve
    Belos::ReturnType ret = belosSolver->solve();
    TEST_EQUALITY(ret, Belos::Converged);
    
    // Get the number of iterations for this solve.
    int numIters = belosSolver->getNumIters();
    switch (TestHelpers::Parameters::getDefaultComm()->getSize()) { 
    case 0: TEST_EQUALITY(numIters, 5); break;
    case 4: 
      // Epetra TEST_EQUALITY(numIters, 6);
      // Tpetra TEST_EQUALITY(numIters, 7); 
      break;
    default:;
    }

    // Test norm equality across the unit tests
    Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(1);
    X->norm2(norms);
    TEST_EQUALITY(MueLuTests::BelosAdaptersTestResults<Scalar>(norms[0]), true);
  }

  //
  // Code factorization:
  //

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BelosAdaptersTest_XpetraOp_XpetraMV(Xpetra::UnderlyingLib lib, Teuchos::FancyOStream & out, bool & success) {
#include "MueLu_UseShortNames.hpp"
    
    RCP<TestProblem<SC,LO,GO,NO,LMO> > p   = TestHelpers::getTestProblem<SC,LO,GO,NO,LMO>(lib);

    typedef Xpetra::MultiVector<SC> MV;
    typedef Belos::OperatorT<MV>    OP;
    
    // Construct a Belos LinearProblem object
    RCP<OP> belosOp   = rcp(new Belos::MueLuOp<SC, LO, GO, NO, LMO>    (p->GetA()));
    RCP<OP> belosPrec = rcp(new Belos::MueLuPrecOp<SC, LO, GO, NO, LMO>(p->GetH()));
    
    // Test adapters
    RCP<MultiVector> X = p->GetNewX0();
    MueLuTests::BelosAdaptersTest<SC, MV, OP>(belosOp, belosPrec, X, p->GetRHS(), out, success);
  }

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

  // TEST:
  // - OP: Xpetra::Operator
  // - MV: Xpetra::MultiVector
  TEUCHOS_UNIT_TEST(BelosAdapters, XpetraOp_XpetraMV) {
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();
    BelosAdaptersTest_XpetraOp_XpetraMV<SC, LO, GO, NO, LMO>(lib, out, success);
  }
 
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
      MueLuTests::BelosAdaptersTest_XpetraOp_XpetraMV<SC, LO, GO, NO, LMO>(Xpetra::UseEpetra, out, success);
    }
  }
#endif

#endif // MUELU_DISABLED

} // namespace MueLuTests
