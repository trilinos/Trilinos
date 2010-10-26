#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_SaPFactory.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(SaPFactory, Test0)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;
  typedef MueLu::SaPFactory<Scalar,LO,GO,Node>    SaPFactory;

  using namespace Teuchos;
  using namespace MueLu;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  TEUCHOS_TEST_EQUALITY(sapFactory != Teuchos::null, true, out, success);

  sapFactory->SetOutputLevel(42);
  TEUCHOS_TEST_EQUALITY(sapFactory->GetOutputLevel() == 42, true, out, success);

  out << *sapFactory << std::endl;

}

TEUCHOS_UNIT_TEST(SaPFactory, GetSetMethods)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;
  typedef MueLu::SaPFactory<Scalar,LO,GO,Node>    SaPFactory;

  using namespace Teuchos;
  using namespace MueLu;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  sapFactory->SetMaxCoarseSize( (GO)55 );
  TEUCHOS_TEST_EQUALITY(((GO)55 ) == sapFactory->GetMaxCoarseSize(), true, out, success);
  sapFactory->SetDampingFactor( (Scalar)4/3 );
  TEUCHOS_TEST_EQUALITY(((Scalar)4/3) == sapFactory->GetDampingFactor(), true, out, success);
  sapFactory->TentativeWithQR(true);
  TEUCHOS_TEST_EQUALITY( sapFactory->TentativeWithQR(), true, out, success);
  sapFactory->ReUseP(true);
  TEUCHOS_TEST_EQUALITY( sapFactory->ReUseP(), true, out, success);
  sapFactory->ReUsePtent(true);
  TEUCHOS_TEST_EQUALITY( sapFactory->ReUsePtent(), true, out, success);
  sapFactory->ReUseAggregates(true);
  TEUCHOS_TEST_EQUALITY( sapFactory->ReUseAggregates(), true, out, success);
  sapFactory->ReUseGraph(true);
  TEUCHOS_TEST_EQUALITY( sapFactory->ReUseGraph(), true, out, success);

} //GetSetMethods


}//namespace <anonymous>

