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

}//namespace <anonymous>

