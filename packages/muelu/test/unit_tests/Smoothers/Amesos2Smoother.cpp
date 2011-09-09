#include <Teuchos_UnitTestHarness.hpp>
#include <Amesos2_config.h>
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include "MueLu_Amesos2Smoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  using namespace Smoother;

  TEUCHOS_UNIT_TEST(Amesos2Smoother, NotSetup)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra)
      {
#ifdef HAVE_AMESOS2_KLU2
        Amesos2Smoother smoother("Klu", Teuchos::ParameterList());
        testApplyNoSetup(smoother, out, success);
#endif
      }
  }

  TEUCHOS_UNIT_TEST(Amesos2Smoother, Apply_Correctness)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra)
      {
#ifdef HAVE_AMESOS2_KLU2
        Amesos2Smoother smoother("Klu", Teuchos::ParameterList());
        testDirectSolver(smoother, out, success);
#endif

#ifdef HAVE_AMESOS2_SUPERLU
        Amesos2Smoother smoother("Superlu", Teuchos::ParameterList());
        testDirectSolver(smoother, out, success);
#endif
      }
  }
  
} // namespace MueLuTests
