#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include "MueLu_IfpackSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

/*
  Comments about tests with hard coded results:
  1) Chebyshev smoothing must pass for any number of processors.
  2) Gauss-Seidel must pass for 1 and 4 processors.
  3) For any processor count except 1 and 4, the Gauss-Seidel test will
  report "passing", but this is only because the Teuchos test macro is skipped.
*/

namespace MueLuTests {

  using namespace Smoother;

  TEUCHOS_UNIT_TEST(IfpackSmoother, NotSetup)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      IfpackSmoother smoother("point relaxation stand-alone", Teuchos::ParameterList());
      testApplyNoSetup(smoother, out, success);

    }
  }

  // Tests interface to Ifpack's Gauss-Seidel preconditioner.
  TEUCHOS_UNIT_TEST(IfpackSmoother, HardCodedResult_GaussSeidel1)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();

      Teuchos::ParameterList paramList;
      paramList.set("relaxation: type", "Gauss-Seidel");
      paramList.set("relaxation: sweeps", (int) 1);
      paramList.set("relaxation: damping factor", (double) 1.0);
      paramList.set("relaxation: zero starting solution", false);
    
      IfpackSmoother smoother("point relaxation stand-alone", paramList);

      ST::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);
    
      switch (comm->getSize()) {
      case 1:
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms, 5.773502691896257e-01,1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;
      } // switch

    }
  }

  // Tests interface to Ifpack's Gauss-Seidel preconditioner.
  TEUCHOS_UNIT_TEST(IfpackSmoother, HardCodedResult_GaussSeidel2)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();

      Teuchos::ParameterList paramList;
      paramList.set("relaxation: type", "Gauss-Seidel");
      paramList.set("relaxation: sweeps", (int) 10);
      paramList.set("relaxation: damping factor", (double) 1.0);
      paramList.set("relaxation: zero starting solution", false);
    
      IfpackSmoother smoother("point relaxation stand-alone", paramList);

      ST::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);

      switch (comm->getSize()) {
      case 1:
        TEST_FLOATING_EQUALITY(residualNorms,8.326553652741774e-02,1e-12);
        break;
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms,8.326553653078517e-02,1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;

      } // switch

    }
  } // GaussSeidel

    // Tests interface to Ifpack's Chebyshev preconditioner.
  TEUCHOS_UNIT_TEST(IfpackSmoother, HardCodedResult_Chebyshev)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      // TODO This test should really calculate the reduction analytically.

      Teuchos::ParameterList paramList;
      paramList.set("chebyshev: degree", (int) 3);
      paramList.set("chebyshev: max eigenvalue", (double) 1.98476);
      paramList.set("chebyshev: min eigenvalue", (double) 1.0);
      paramList.set("chebyshev: ratio eigenvalue", (double) 20);
      paramList.set("chebyshev: zero starting solution", false);

      IfpackSmoother smoother("Chebyshev", paramList);

      ST::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);

      TEST_FLOATING_EQUALITY(residualNorms, 5.269156e-01, 1e-7); // Compare to residual reported by ML

    }
  } // Chebyshev

  // Tests interface to Ifpack's ILU(0) preconditioner.
  TEUCHOS_UNIT_TEST(IfpackSmoother, HardCodedResult_ILU)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      // FIXME this will probably fail in parallel b/c it becomes block Jacobi

      IfpackSmoother smoother("ILU", Teuchos::ParameterList());

      ST::magnitudeType residualNorms = testApply_A125_X0_RandomRHS(smoother, out, success);

      RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
      if (comm->getSize() == 1) {
        TEST_EQUALITY(residualNorms < 1e-10, true);
      } else {
        out << "Pass/Fail is only checked in serial." << std::endl;
      }

    }
  } // ILU

} // namespace MueLuTests

//TODO: test if 10 its of ifpack == 10 apply call
//TODO:     TEST_EQUALITY(smoother->GetNIts(),10);
//TODO: check failed if TPETRA
