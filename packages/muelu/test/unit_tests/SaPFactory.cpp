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
    TEST_EQUALITY(sapFactory != Teuchos::null, true);

    out << *sapFactory << std::endl;

  }

  TEUCHOS_UNIT_TEST(SaPFactory, GetSetMethods)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
    sapFactory->SetDampingFactor( (Scalar)4/3 );
    TEST_EQUALITY(((Scalar)4/3) == sapFactory->GetDampingFactor(), true);
    sapFactory->SetDiagonalView("roomWithAView");
    TEST_EQUALITY( sapFactory->GetDiagonalView(), "roomWithAView");

  } //GetSetMethods


}//namespace MueLuTests

