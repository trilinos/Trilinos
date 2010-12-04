#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_RAPFactory.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(RAPFactory, Test0)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;
  typedef Kokkos::DefaultKernels<Scalar,LO,Node>::SparseOps LMO;

  typedef MueLu::RAPFactory<Scalar,LO,GO,Node,LMO>    RAPFactory;

  using namespace Teuchos;
  using namespace MueLu;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<RAPFactory> rapFactory = rcp(new RAPFactory);
  TEUCHOS_TEST_EQUALITY(rapFactory != Teuchos::null, true, out, success);

  rapFactory->SetOutputLevel(Teuchos::VERB_MEDIUM);
  TEUCHOS_TEST_EQUALITY(rapFactory->GetOutputLevel() == Teuchos::VERB_MEDIUM, true, out, success);

  out << *rapFactory << std::endl;

}

}//namespace <anonymous>

