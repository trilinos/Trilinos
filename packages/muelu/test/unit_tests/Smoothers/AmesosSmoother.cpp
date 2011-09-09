#include <Teuchos_UnitTestHarness.hpp>
#include <Amesos_config.h>
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include "MueLu_AmesosSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  using namespace TestHelpers::Smoothers;

  TEUCHOS_UNIT_TEST(AmesosSmoother, NotSetup)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
      {
#ifdef HAVE_AMESOS_KLU
        AmesosSmoother smoother;
        testApplyNoSetup(smoother, out, success);
#endif
      }
  }

  TEUCHOS_UNIT_TEST(AmesosSmoother, Apply_Correctness)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
      {
#ifdef HAVE_AMESOS_KLU
        AmesosSmoother smoother("Klu");
        testDirectSolver(smoother, out, success);
#endif

#ifdef HAVE_AMESOS_SUPERLU
        AmesosSmoother smoother("Superlu");
        testDirectSolver(smoother, out, success);
#endif
      }
  }
  
} // namespace MueLuTests

//TODO
// Test with invalid std::string
// Test with invalid parameterList? (== a characterization test for Amesos)
// Test if paramList takes into account
// Check if all the defaults that are used by MueLu are tested
