#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include "MueLu_AmesosSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  using namespace Smoother;

  TEUCHOS_UNIT_TEST(AmesosSmoother, NotSetup)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
      {

        AmesosSmoother smoother("Amesos_Klu", Teuchos::ParameterList());
        testApplyNoSetup(smoother, out, success);

      }
  }

  TEUCHOS_UNIT_TEST(AmesosSmoother, Apply_Correctness_KLU)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
      {

        AmesosSmoother smoother("Amesos_Klu", Teuchos::ParameterList());
        ST::magnitudeType residualNorms = testApply_A125_X0_RandomRHS(smoother, out, success);
        TEST_EQUALITY(residualNorms < 1e-12, true);

      }
  }
  
} // namespace MueLuTests

//TODO
// Test with SuperLU
// Test with invalid std::string
// Test if paramList takes into account
