#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_SaPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

#include "MueLu_ExtendedHashtable.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(ExtendedHashtable, rcpTests)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
    TEUCHOS_TEST_EQUALITY(sapFactory != Teuchos::null, true, out, success);

    RCP<SaPFactory> sapFactory2 = rcp(new SaPFactory);
    TEUCHOS_TEST_EQUALITY(sapFactory2 != Teuchos::null, true, out, success);

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // build Operator
    Teuchos::ParameterList params;
    const RCP<const Map> map = MapFactory::Build(Xpetra::UseTpetra, 20, 0, comm);
    RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>("Laplace1D", map, params);

    // an extended hashtable
    RCP<MueLu::UTILS::ExtendedHashtable> exh = Teuchos::rcp(new MueLu::UTILS::ExtendedHashtable());

    exh->Set<RCP<Operator> >("op",Op,sapFactory);
    RCP<Operator> test = exh->Get<RCP<Operator> > ("op",sapFactory);
    TEUCHOS_TEST_EQUALITY_CONST( test, Op, out, success );


    exh->Set<RCP<Operator> >("op",Op,sapFactory2);
    test = exh->Get<RCP<Operator> > ("op",sapFactory2);
    TEUCHOS_TEST_EQUALITY_CONST( test, Op, out, success );


    exh->Set("op2",24,Teuchos::null);
    int test2 = exh->Get<int>("op2",Teuchos::null);
    TEUCHOS_TEST_EQUALITY_CONST( test2, 24, out, success );

    exh->Set("op2",12,Teuchos::null);
    test2 = exh->Get<int>("op2",Teuchos::null);
    TEUCHOS_TEST_EQUALITY_CONST( test2, 12, out, success );

    exh->Remove("op",sapFactory2);
    exh->Remove("op",sapFactory);

    exh->Set<std::string>("op","xxx",sapFactory2);
    std::string test3 = exh->Get<std::string> ("op",sapFactory2);
    TEUCHOS_TEST_EQUALITY_CONST( test3, "xxx", out, success );

  }

  /*TEUCHOS_UNIT_TEST(SaPFactory, GetSetMethods)
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

  } //GetSetMethods*/


}//namespace MueLuTests

