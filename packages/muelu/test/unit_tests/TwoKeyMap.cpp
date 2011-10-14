#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_TentativePFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

#include "MueLu_TwoKeyMap.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(TwoKeyMap, ptrTests)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<TentativePFactory> sapFactory = rcp(new TentativePFactory);
    TEST_EQUALITY(sapFactory != Teuchos::null, true);

    RCP<TentativePFactory> sapFactory2 = rcp(new TentativePFactory);
    TEST_EQUALITY(sapFactory2 != Teuchos::null, true);

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // build Operator
    Teuchos::ParameterList params;
    const RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), 20, 0, comm);
    RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>("Laplace1D", map, params);

    // an extended hashtable
    RCP<MueLu::UTILS::TwoKeyMap> exh = Teuchos::rcp(new MueLu::UTILS::TwoKeyMap());

    exh->Set<RCP<Operator> >("op", sapFactory.get(), Op);
    RCP<Operator> test = Teuchos::null;
    test = exh->Get<RCP<Operator> >("op", sapFactory.get());
    TEST_EQUALITY_CONST( test, Op );

    exh->Set("op", sapFactory.get(), 22);
    int test2 = exh->Get<int> ("op", sapFactory.get());
    TEST_EQUALITY_CONST( test2, 22 );

    exh->Set("op", sapFactory2.get(), Op);
    RCP<Operator> test3 = exh->Get<RCP<Operator> > ("op", sapFactory2.get());
    TEST_EQUALITY_CONST( test3, Op );
    TEST_EQUALITY_CONST( exh->GetType("op", sapFactory.get()), "int" );

    exh->Remove("op", sapFactory2.get());
    TEST_EQUALITY_CONST( exh->IsKey("op", sapFactory.get()),  true );
    TEST_EQUALITY_CONST( exh->IsKey("op", sapFactory2.get()), false );

    exh->Remove("op", sapFactory.get());
    TEST_EQUALITY_CONST( exh->IsKey("op", sapFactory.get()),  false );
    TEST_EQUALITY_CONST( exh->IsKey("op", sapFactory2.get()), false );

  }

}//namespace MueLuTests

