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
    sapFactory->TentativeWithQR(true);
    TEUCHOS_TEST_EQUALITY( sapFactory->TentativeWithQR(), true, out, success);
    sapFactory->ReUseP(true);
    TEUCHOS_TEST_EQUALITY( sapFactory->ReUseP(), true, out, success);
    sapFactory->ReUsePtent(true);
    TEUCHOS_TEST_EQUALITY( sapFactory->ReUsePtent(), true, out, success);
    sapFactory->SetDiagonalView("roomWithAView");
    TEUCHOS_TEST_EQUALITY( sapFactory->GetDiagonalView(), "roomWithAView", out, success);
    TEST_THROW( sapFactory->SetUseAFiltered(true), MueLu::Exceptions::NotImplemented ); //FIXME

  } //GetSetMethods


}//namespace MueLuTests

