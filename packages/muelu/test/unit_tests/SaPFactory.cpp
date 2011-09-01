#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_SaPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(SaPFactory, Test0)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
    TEUCHOS_TEST_EQUALITY(sapFactory != Teuchos::null, true, out, success);

    out << *sapFactory << std::endl;

  }

  TEUCHOS_UNIT_TEST(SaPFactory, GetSetMethods)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
    sapFactory->SetDampingFactor( (Scalar)4/3 );
    TEUCHOS_TEST_EQUALITY(((Scalar)4/3) == sapFactory->GetDampingFactor(), true, out, success);
    sapFactory->SetDiagonalView("roomWithAView");
    TEUCHOS_TEST_EQUALITY( sapFactory->GetDiagonalView(), "roomWithAView", out, success);

  } //GetSetMethods


}//namespace MueLuTests

